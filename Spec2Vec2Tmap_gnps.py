import os
import pandas as pd
from matchms.filtering import add_losses
from matchms.filtering import add_parent_mass
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms.filtering import reduce_to_number_of_peaks
from matchms.filtering import require_minimum_number_of_peaks
from matchms.filtering import select_by_mz
from matchms.filtering import derive_smiles_from_inchi

from matchms.importing import load_from_mgf
from spec2vec import SpectrumDocument

from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from mhfp.encoder import MHFPEncoder

import re
from collections import Counter
import numpy as np
import tmap as tm
from faerun import Faerun
import scipy.stats as ss

from matplotlib import pyplot as plt

from sklearn import preprocessing



### PART 1: IMPORT AND CLEAN SPECTRA + METADATA and CONVERT SPECTRA TO DOCUMENTS ###

# Spectra processing: https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html for help
def apply_my_filters(s):
    """This is how one would typically design a desired pre- and post-
    processing pipeline."""
    s = default_filters(s)
    s = add_parent_mass(s)
    s = normalize_intensities(s)
    s = reduce_to_number_of_peaks(s, n_required=10, ratio_desired=0.5)
    s = select_by_mz(s, mz_from=0, mz_to=1000)
    s = add_losses(s, loss_mz_from=10.0, loss_mz_to=200.0)
    s = require_minimum_number_of_peaks(s, n_required=10)
    return s

# Load data from MGF file and apply filters

spectrums = [apply_my_filters(s) for s in load_from_mgf("/Users/pma/Dropbox/Research_UNIGE/git_repos/pf_project/data/spectral/VGF_20plates_MN_spectra.mgf")]

# Omit spectrums that didn't qualify for analysis
spectrums = [s for s in spectrums if s is not None]

# Create spectrum documents
reference_documents = [SpectrumDocument(s) for s in spectrums]

# Create a df with the metdata from all spectra
metadata_df = pd.DataFrame(s.metadata for s in spectrums)

# At this stage we optionally import features metadata

### Loading the data 
# optional_meta_annot = pd.read_csv("/Users/pma/Dropbox/Research_UNIGE/git_repos/pf_project/data/chemspace/VGF_20plates_MN_results_top1_repond_classyfied_DW.txt", sep = '\t')
# optional_meta_pl = pd.read_csv("/Users/pma/tmp/Full_meta_mean_VGF_20_plates.tsv", sep = '\t')
optional_meta_bio = pd.read_csv("/Users/pma/tmp/VGF_20_plates_Bioactivity_scores.csv", sep = ',')
optional_meta_top5 = pd.read_csv("/Users/pma/VGF_20plates_MN_results_top5_repond_classyfired.out", sep = '\t')

optional_meta_bio.rename(columns={'Unnamed: 0': 'feature_id'}, inplace=True)

optional_meta_all = pd.merge(optional_meta_top5, optional_meta_bio, left_on = "ID_ISDB", right_on = "feature_id",  how="left")


### converting the type to int
metadata_df.scans = metadata_df.scans.astype(int)


metadata_df = pd.merge(metadata_df, optional_meta_all, left_on = "scans", right_on = "ID_ISDB",  how="left")


metadata_df['label'] = 'node_' + metadata_df.index.astype(str)

metadata_df.to_csv('metadata_pf.tsv', sep = '\t')


metadata_df.loc[metadata_df['Score_Tcruzi_10ug/ml_subtracted'].idxmin()]

metadata_df.rename(columns={'Score_Tcruzi_10ug/ml_subtracted': 'Score_Tcruzi_10ug_ml_subtracted'}, inplace=True)


metadata_df['Cruzi_cat'] = pd.cut(metadata_df.Score_Tcruzi_10ug_ml_subtracted, bins=[-10000,-7500,-5000,-2500,0,2500, 5000, 7500, 10000],labels=['H','G','F','E','D','C','B','A'])


metadata_df.info()


# Create a list of all our documents-spectra 
texts = []
for doc in reference_documents:
    texts.append(doc.words)

len(texts)


texts[0]

# Here we convert the text list as a df. 
# the objective is to add a biological value to the fragments and neutral losses in order to have a bioinformed spectral tmap

texts_df = pd.DataFrame(texts)

texts_df_informed = pd.merge(texts_df, metadata_df, left_index=True, right_index=True)


cols = texts_df_informed.columns[pd.to_numeric(texts_df_informed.columns, errors='coerce').to_series().notnull()]

cols_list = list(cols)

cols_list += ['Cruzi_cat']

for i in cols:
    texts_df_informed[i] = texts_df_informed[i] + texts_df_informed['Cruzi_cat'].astype(str)


# now we want to get this back as a list of list 

text_informed_sep = texts_df_informed[cols_list].values.tolist()


len(text_informed_sep)

### PART 2: SPECTRAL T-MAP ###


enc = tm.Minhash()
lf = tm.LSHForest()

# (Dirty) trick to pass from a counter object to a list
ctr = Counter()
for text in text_informed_sep:
    ctr.update(text)
n = 0
ctr = ctr.most_common()[: -(len(ctr) - n) - 1 : -1]

all_words = {}
for i, (key, _) in enumerate(ctr):
    all_words[key] = i

fingerprints = []
for text in text_informed_sep:
    fingerprint = []
    for t in text:
        if t in all_words:
            fingerprint.append(all_words[t])
    fingerprints.append(tm.VectorUint(fingerprint))

lf.batch_add(enc.batch_from_sparse_binary_array(fingerprints))


lf.index()
config = tm.LayoutConfiguration()
config.k = 100
x, y, s, t, _ = tm.layout_from_lsh_forest(lf, config=config)



### Export to cytoscape

coordinates_tuples = list(zip(x,y))
node_df = pd.DataFrame(coordinates_tuples, columns=['x','y'])
node_df['label'] = node_df.index
node_df['label'] = 'node_' + node_df['label'].astype(str)
node_df['cluster'] = 1
node_df = node_df[['label', 'cluster', 'x', 'y']]
node_df = node_df.round({'x': 6, 'y': 6})


edge_tuples = list(zip(s,t))
edge_df = pd.DataFrame(edge_tuples, columns=['s','t'])

node_df.to_csv('node_coordinates_bioinformed.csv', sep = ' ', index=True, header=False)
edge_df.to_csv('edges_bioinformed.csv', sep = ' ', index=False, header=False)




type_labels, type_data = Faerun.create_categories(metadata_df["Origin Type"])

genera_labels, genera_data = Faerun.create_categories(metadata_df["Family_Bio"])
top_genera = [i for i, _ in Counter(genera_data).most_common(9)]

top_genera_labels = [(7, "Other")]
genera_map = [7] * len(genera_data)
value = 1
for i, name in genera_labels:
    if i in top_genera:
        v = value
        if v == 7:
            v = 0
        top_genera_labels.append((v, name))
        genera_map[i] = v
        value += 1

genera_data = [genera_map[val] for _, val in enumerate(genera_data)]



tab_10 = plt.cm.get_cmap("tab10")
colors = [i for i in tab_10.colors]
colors[7] = (0.17, 0.24, 0.31)
tab_10.colors = tuple(colors)

# Compute molecular descriptors  
hac = []
c_frac = []
ring_atom_frac = []
largest_ring_size = []

for i, row in metadata_df.iterrows():
    mol = AllChem.MolFromSmiles(row["Smiles_DNP"])
    if mol: 
        atoms = mol.GetAtoms()
        size = mol.GetNumHeavyAtoms()
        n_c = 0
        n_ring_atoms = 0
        for atom in atoms:
            if atom.IsInRing():
                n_ring_atoms += 1
            if atom.GetSymbol().lower() == "c":
                n_c += 1
        c_frac.append(n_c / size)
        ring_atom_frac.append(n_ring_atoms / size)
        sssr = AllChem.GetSymmSSSR(mol)
        if len(sssr) > 0:
            largest_ring_size.append(max([len(s) for s in sssr]))
        else:
            largest_ring_size.append(0)
        hac.append(size)

c_frak_ranked = ss.rankdata(np.array(c_frac) / max(c_frac)) / len(c_frac)

# Creating the smiles label
metadata_df["label_smiles"] = (
    metadata_df["smiles"]
    + '__'
    + metadata_df.smiles.astype(str)
    + '__'
    + "</a>"
    )

# Colors 
tab_10 = plt.cm.get_cmap("tab10")
colors = [i for i in tab_10.colors]
colors[7] = (0.17, 0.24, 0.31)
tab_10.colors = tuple(colors)


tab_10 = plt.cm.get_cmap("viridis")
colors = [i for i in tab_10.colors]
colors[7] = (0.17, 0.24, 0.31)
tab_10.colors = tuple(colors)

metadata_df.info()

metadata_df.fillna('None',inplace=True)




# Plotting the spectral t-map
f = Faerun(view="front", coords=False)
f.add_scatter(
        "Spectral_Tmap",
        {
        "x": x,
        "y": y,
        "c": [
            metadata_df["parent_mass"],
            genera_data,
            hac,
            c_frak_ranked,
            ring_atom_frac,
            largest_ring_size,
        ],
        "labels": metadata_df["Smiles_DNP"],
    },
    shader="smoothCircle",
    point_scale= 5.0,
    max_point_size= 20,
    categorical= [False, True, False, False, False, False],
    colormap= ["viridis", "tab10", "rainbow", "rainbow", "rainbow", "Blues"],
    legend_labels=[None,top_genera_labels],
    series_title= ["Parent mass", 
    "Family_Bio","HAC","C Frac","Ring Atom Frac","Largest Ring Size",],
    has_legend= True,
)
f.add_tree("Spectral_Tmap_tree", {"from": s, "to": t}, point_helper="Spectral_Tmap")
f.plot("Spectral_Tmap", template="smiles")




# Plotting a plain spectral t-map
f = Faerun(view="front", coords=False)
f.add_scatter(
        "Spectral_Tmap",
        {
        "x": x,
        "y": y,
        "c": [
            metadata_df["parent_mass"],
            metadata_df["Score_Tcruzi_10ug_ml_subtracted"],
            # hac,
            # c_frak_ranked,
            # ring_atom_frac,
            # largest_ring_size,
        ],
        "labels": metadata_df["ID_ISDB"],
    },
    shader="smoothCircle",
    point_scale= 5.0,
    max_point_size= 20,
    categorical= [False, False],
    colormap= ["viridis", "viridis"],
    legend_labels=[None],
    series_title= ["Parent mass", "Cruzi_cat"],
    has_legend= True,
)
f.add_tree("Spectral_Tmap_tree", {"from": s, "to": t}, point_helper="Spectral_Tmap")
f.plot("Spectral_Tmap", template="default")


# Plotting a plain spectral t-map
f = Faerun(view="front", coords=False)
f.add_scatter(
        "Spectral_Tmap_bioinformed_sep",
        {
        "x": x,
        "y": y,
        "c": [
            metadata_df["parent_mass"],
            metadata_df["Score_Tcruzi_10ug_ml_subtracted"],
            # hac,
            # c_frak_ranked,
            # ring_atom_frac,
            # largest_ring_size,
        ],
        "labels": metadata_df["Consensus_ci_cf"],
    },
    shader="smoothCircle",
    point_scale= 5.0,
    max_point_size= 20,
    categorical= [False, False],
    colormap= ["viridis", "viridis"],
    legend_labels=[None],
    series_title= ["Parent mass", "Cruzi_cat"],
    has_legend= True,
)
f.add_tree("Spectral_Tmap_bioinformed_tree_sep", {"from": s, "to": t}, point_helper="Spectral_Tmap_bioinformed_sep")
f.plot("Spectral_Tmap_bioinformed_sep", template="default")



### PART 3: STRUCTURAL T-MAP ###


def main():
    """ Main funciton """

    enc = MHFPEncoder(1024)
    lf = tm.LSHForest(1024, 64)

    fps = []

    for i, row in metadata_df.iterrows():
        if i != 0 and i % 1000 == 0:
            print(100 * i / len(metadata_df))

        mol = AllChem.MolFromSmiles(row["smiles"])
        fps.append(tm.VectorUint(enc.encode_mol(mol)))

    lf.batch_add(fps)
    lf.index()

    c_frak_ranked = ss.rankdata(np.array(c_frac) / max(c_frac)) / len(c_frac)

    cfg = tm.LayoutConfiguration()
    cfg.node_size = 1 / 26
    cfg.mmm_repeats = 2
    cfg.sl_extra_scaling_steps = 5
    cfg.k = 20
    cfg.sl_scaling_type = tm.RelativeToAvgLength
    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)


    tab_10 = plt.cm.get_cmap("tab10")
    colors = [i for i in tab_10.colors]
    colors[7] = (0.17, 0.24, 0.31)
    tab_10.colors = tuple(colors)

    f = Faerun(view="front", coords=False)
    f.add_scatter(
        "Structural_tmap",
        {
            "x": x,
            "y": y,
            "c": [
                metadata_df["parent_mass"],
                hac,
                c_frak_ranked,
                ring_atom_frac,
                largest_ring_size,
            ],
            "labels": metadata_df["label_smiles"],
        },
        shader="smoothCircle",
        point_scale=8.0,
        max_point_size=20,
        categorical=[False, False, False, False, False],
        colormap=["rainbow", "rainbow", "rainbow", "rainbow", "Blues"],
        series_title=[
            "Parent mass",
            "HAC",
            "C Frac",
            "Ring Atom Frac",
            "Largest Ring Size",
        ],
        has_legend=True,
    )
    f.add_tree("Structural_tmap_tree", {"from": s, "to": t}, point_helper="Structural_tmap")
    f.plot('Structural_tmap', template="smiles")


if __name__ == "__main__":
    main()