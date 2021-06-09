import pickle
import sys
from collections import Counter
from time import sleep

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scipy.stats as ss
import tmap as tm
from faerun import Faerun
from map4 import MAP4Calculator
from rdkit.Chem import AllChem


# Define useful functions:

def Top_N_classes(N=0, input_data=[], input_labels=[]):
    """Keep only the top N classes for the input classes and replace the following by 'Other', default N = 0 return input"""
    if N == 0:
        return input_data, input_labels
    else:
        top_input = [i for i, _ in Counter(input_data).most_common(N)]
        output_labels = [(15, "Other")]
        input_map = [15] * len(input_data)
        value = 1
        for i, name in input_labels:
            if i in top_input:
                v = value
                if v == 15:
                    v = 0
                output_labels.append((v, name))
                input_map[i] = v
                value += 1
        output_data = [input_map[val] for _, val in enumerate(input_data)]
        return output_data, output_labels


def keep_only_given_class(to_keep=[], input_data=[], input_labels=[]):
    """Keep only the class specified in to_keep argument for mapping"""
    output_labels = []
    output_data = []
    i_to_keep = []
    values = [*range(1, len(to_keep) + 1, 1)]
    d = dict(zip(to_keep, values))
    print(d)
    d2 = {}

    for i, name in input_labels:
        if name in to_keep:
            output_labels.append(tuple([d[name], name]))
            i_to_keep.append(i)
            d2[i] = d[name]
        else:
            output_labels.append(tuple([0, 'Other']))
    output_labels = list(set(output_labels))

    print(d2)
    for val in input_data:
        if val in i_to_keep:
            output_data.append(val)
        else:
            output_data.append(0)
    output_data = [d2.get(item, item) for item in output_data]
    return output_data, output_labels


########################################################
### PART 1: LOAD Lotus and transform for plotting ###
########################################################

# Load lotus
df_meta = pd.read_csv(
    '210523_lotus_dnp_metadata.csv',
    usecols=['structure_smiles_2D', 'structure_taxonomy_npclassifier_01pathway',
    'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class', 
    'structure_taxonomy_classyfire_01kingdom', 'structure_taxonomy_classyfire_02superclass', 
    'structure_taxonomy_classyfire_03class'],
    sep=",")

# Keep only 1 occurence of structure - biosource pair
df_meta = df_meta.drop_duplicates(
    subset=['structure_smiles_2D'])

########################################################
###           PART 2: LOAD annotations               ###
########################################################

# Load annotations
df_annotations = pd.read_csv('Table_PF_pos_190905.csv', usecols=['CRC_Number_DNP', 'Smiles_DNP'], low_memory=False, sep=",")

df_annotations['Smiles'] = df_annotations['Smiles_DNP'].str.split('|').str[0]
df_annotations['CRC_Number_DNP'] = df_annotations['CRC_Number_DNP'].str.split('|').str[0]
df_annotations.dropna(subset=['Smiles'], inplace=True)

anisomeric_smiles = []
inchikey = []
for i, row in df_annotations.iterrows():
    mol = AllChem.MolFromSmiles(row["Smiles"])
    if mol != None:
        aniso_smiles = AllChem.MolToSmiles(mol, isomericSmiles = False)
        ik = AllChem.MolToInchiKey(mol)
        anisomeric_smiles.append(aniso_smiles)
        inchikey.append(ik)
    else:
        anisomeric_smiles.append(np.nan)
        inchikey.append(aniso_smiles)

df_annotations['anisomeric_smiles'] = anisomeric_smiles
df_annotations['inchikey'] = inchikey
df_annotations['inchikey_2D'] = df_annotations['inchikey'].str[:14]

df_annotations = df_annotations.drop_duplicates(
    subset=['inchikey_2D'])
    
df_annotations.dropna(subset=['anisomeric_smiles'], inplace=True)

df_annotations.drop(['Smiles', 'Smiles_DNP'], axis=1, inplace=True)

print('Number of unique annotations (2D structures): ' + str(len(df_annotations)))

# Load classyfire classes for DNP annotations:

df_classes_dnp = pd.read_csv('190602_DNP_TAXcof_CF.tsv', sep='\t', usecols= ['CRC_Number_DNP', 'Kingdom_cf_DNP', 'Superclass_cf_DNP', 'Class_cf_DNP'])
df_classes_dnp.drop_duplicates(subset='CRC_Number_DNP', inplace=True)

df_annotations = df_annotations.merge(df_classes_dnp, on = 'CRC_Number_DNP', how ='left')
df_annotations.drop(['CRC_Number_DNP'], axis=1, inplace=True)


# Merge both tables:

anisomeric_smiles = []
inchikey = []
for i, row in df_meta.iterrows():
    mol = AllChem.MolFromSmiles(row["structure_smiles_2D"])
    if mol != None:
        ik = AllChem.MolToInchiKey(mol)
        inchikey.append(ik)
    else:
        anisomeric_smiles.append(np.nan)
        inchikey.append(aniso_smiles)

df_meta['inchikey'] = inchikey
df_meta['inchikey_2D'] = df_meta['inchikey'].str[:14]
df_meta = df_meta.drop_duplicates(
    subset=['inchikey_2D'])


df_merged = pd.merge(df_meta, df_annotations, left_on = 'inchikey_2D', right_on='inchikey_2D', how='outer')
df_merged['is_annotated'] = np.where(df_merged['anisomeric_smiles'].isnull(), 'no', 'yes')

df_merged['structure_smiles_2D'].fillna(df_merged['anisomeric_smiles'], inplace=True)
df_merged['structure_taxonomy_classyfire_01kingdom'].fillna(df_merged['Kingdom_cf_DNP'], inplace=True)
df_merged['structure_taxonomy_classyfire_02superclass'].fillna(df_merged['Superclass_cf_DNP'], inplace=True)
df_merged['structure_taxonomy_classyfire_03class'].fillna(df_merged['Class_cf_DNP'], inplace=True)

df_merged.drop(['anisomeric_smiles', 'Kingdom_cf_DNP', 'Superclass_cf_DNP', 'Class_cf_DNP'], axis=1, inplace=True)

values = {
    # 'structure_taxonomy_npclassifier_01pathway_first': 'Unknown',
    # 'structure_taxonomy_npclassifier_02superclass_first': 'Unknown',
    # 'structure_taxonomy_npclassifier_03class_first': 'Unknown',
    'structure_taxonomy_classyfire_01kingdom': 'Unknown',
    'structure_taxonomy_classyfire_02superclass': 'Unknown',
    'structure_taxonomy_classyfire_03class': 'Unknown',
          }
df_merged.fillna(value=values, inplace=True)

# Generating class for plotting
dic_categories = {
    "structure_taxonomy_classyfire_01kingdom": {'Ncat': 0},
    "structure_taxonomy_classyfire_02superclass": {'Ncat': 19},
    "structure_taxonomy_classyfire_03class": {'Ncat': 19},
    "is_annotated": {'Ncat': 0}
}

for dic in dic_categories:
    labels, data = Faerun.create_categories(df_merged[str(dic)])
    dic_categories[dic]['data'], dic_categories[dic]['labels'] = Top_N_classes(dic_categories[dic]['Ncat'], data, labels)

# # Publication examples, but you can use it as you want!
# labels, data = Faerun.create_categories(df_gb['organism_taxonomy_06family_<lambda_0>'])
# simaroubaceae_data, simaroubaceae_labels = keep_only_given_class(['Simaroubaceae'], data, labels)

# labels, data = Faerun.create_categories(df_gb['structure_taxonomy_npclassifier_03class_first'])
# NPclass_data, NPclass_labels = keep_only_given_class([
#     'Oleanane triterpenoids', 'Germacrane sesquiterpenoids', 'Carotenoids (C40, β-β)',
#     'Flavonols', 'Lanostane, Tirucallane and Euphane triterpenoids', 'Cyclic peptides',
#     'Guaiane sesquiterpenoids', 'Flavones', 'Colensane and Clerodane diterpenoids',
#     'Quassinoids'], data, labels)

# count_simaroubaceae = df_gb.organism_taxonomy_06family_join.str.count("Simaroubaceae")
# simaroubaceae_specificity = count_simaroubaceae / df_gb['biosource_count']

# # Generating colormaps for plotting
# cmap = mcolors.ListedColormap(["gainsboro", "peachpuff", "salmon", "tomato"])
# cmap2 = mcolors.ListedColormap(["gainsboro", "tomato"])
# cmap3 = mcolors.ListedColormap(["lightgray", "firebrick"])
# cmap4 = mcolors.ListedColormap(
#     ["gainsboro", "#83C644", "#04AC6E", "#00ff99", "#008ECA", "#008C83", "#D61E1E", "#006A7A", "#00B2C4", "#2F4858",
#      "#EAD400"])

# Generate a labels column
df_merged["labels"] = (
        df_merged["structure_smiles_2D"]
        + '__'
        + df_merged["structure_taxonomy_classyfire_03class"]
        + "</a>"
)

########################################################
# OPTION 1: START FROM SCRATCH AND GENERATE LSH FOREST AND TMAP FROM SMILES, PLUS CHEMICAL DESCRIPTORS
########################################################

MAP4 = MAP4Calculator(dimensions=1024)
# ENC = tm.Minhash(1024)

# enc = MHFPEncoder(1024)
lf = tm.LSHForest(1024, 64)

fps = []
hac = []
c_frac = []
ring_atom_frac = []
largest_ring_size = []

n_iter = len(df_merged)
for i, row in df_merged.iterrows():
    j = (i + 1) / n_iter
    sys.stdout.write('\r')
    sys.stdout.write(f"[{'=' * int(50 * j):{50}s}] {round((100 * j), 2)}%")
    sys.stdout.flush()
    sleep(0.05)

    mol = AllChem.MolFromSmiles(row["structure_smiles_2D"])
    atoms = mol.GetAtoms()
    size = mol.GetNumHeavyAtoms()
    n_c = 0
    n_ring_atoms = 0
    for atom in atoms:
        if atom.IsInRing():
            n_ring_atoms += 1
        if atom.GetSymbol().lower() == "c":
            n_c += 1

    c_frac.append(n_c / size if size != 0 else 0)

    ring_atom_frac.append(n_ring_atoms / size if size != 0 else 0)

    sssr = AllChem.GetSymmSSSR(mol)
    if len(sssr) > 0:
        largest_ring_size.append(max([len(s) for s in sssr]))
    else:
        largest_ring_size.append(0)
    hac.append(size)
    fps.append(mol)

fps = MAP4.calculate_many(fps)
lf.batch_add(fps)
lf.index()

# Store lsh forest and structure metadata
lf.store("210608_pf_set_annotation_vs_lotus.dat")
with open("210608_pf_set_annotation_vs_lotus_attributes.pickle", "wb+") as f:
    pickle.dump(
        (hac, c_frac, ring_atom_frac, largest_ring_size),
        f,
        protocol=pickle.HIGHEST_PROTOCOL,
    )

# tmap configuration
cfg = tm.LayoutConfiguration()
cfg.k = 20
cfg.sl_extra_scaling_steps = 10
cfg.node_size = 1 / 50
cfg.mmm_repeats = 2
cfg.sl_repeats = 2

# tmap generation
x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)

# To store coordinates
x = list(x)
y = list(y)
s = list(s)
t = list(t)
pickle.dump(
    (x, y, s, t), open("210608_pf_set_annotation_vs_lotus_coords.dat", "wb+"), protocol=pickle.HIGHEST_PROTOCOL
)

del (lf)

########################################################
# OPTION 2: LOAD PRE_COMPUTED LSH FOREST AND CHEMICAL DESCRIPTORS
########################################################

lf = tm.LSHForest(1024, 64)
lf.restore(
    "../../data/interim/tmap/210506_lotus_2D_map4.dat")  # Version "210312_lotus.dat" contains 270'336 resulting from 210223_frozen_metadata groupby(structure_wikidata)

hac, c_frac, ring_atom_frac, largest_ring_size = pickle.load(
    open("../../data/interim/tmap/210506_lotus_2D_map4.pickle", "rb")
)

# tmap configuration

cfg = tm.LayoutConfiguration()
cfg.k = 20
cfg.sl_extra_scaling_steps = 10
cfg.node_size = 1 / 50
cfg.mmm_repeats = 2
cfg.sl_repeats = 2

# tmap generation
x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)

# To store coordinates
# x = list(x)
# y = list(y)
# s = list(s)
# t = list(t)
# pickle.dump(
#     (x, y, s, t), open("../data/interim/tmap/210505_coords_lotus_2D_map4.dat", "wb+"), protocol=pickle.HIGHEST_PROTOCOL
# )

del (lf)

########################################################
# OPTION 3: LOAD PRE-COMPUTED COORDINATES, SOURCES AND TARGETS
# AND CHEMICAL DESCRIPTORS
########################################################

x, y, s, t = pickle.load(open("210608_pf_set_annotation_vs_lotus_coords.dat",
                              "rb"))

hac, c_frac, ring_atom_frac, largest_ring_size = pickle.load(
    open("210608_pf_set_annotation_vs_lotus_attributes.pickle", "rb")
)

########################################################
# COMMON PART: PLOT THE TMAP
########################################################
c_frak_ranked = ss.rankdata(np.array(c_frac) / max(c_frac)) / len(c_frac)

# Plotting function
f = Faerun(view="front", coords=False)#clear_color='#ffffff'
f.add_scatter(
    "annotations",
    {
        "x": x,
        "y": y,
        "c": [
            dic_categories['structure_taxonomy_classyfire_01kingdom']['data'],
            dic_categories['structure_taxonomy_classyfire_02superclass']['data'],
            dic_categories['structure_taxonomy_classyfire_03class']['data'],
            dic_categories['is_annotated']['data']
        ],
        "labels": df_merged["labels"],
    },
    shader="smoothCircle",
    point_scale=2.0,
    max_point_size=10,
    legend_labels=[
        dic_categories['structure_taxonomy_classyfire_01kingdom']['labels'],
        dic_categories['structure_taxonomy_classyfire_02superclass']['labels'],
        dic_categories['structure_taxonomy_classyfire_03class']['labels'],
        dic_categories['is_annotated']['labels']
        ],
    categorical=[True, True, True, True],
    colormap=["tab20", "tab20", "tab20", "tab20"],
    series_title=[ "chemo_cf_kingdom", "chemo_cf_superclass", "chemo_cf_class", "annotation_status"],
    has_legend=True
)

f.add_tree("annotation_tree", {"from": s, "to": t}, point_helper="annotations", color='#e6e6e6')
f.plot('210609_annotation_vs_lotusdnp_tmap', template="smiles")


            # hac,
            # c_frak_ranked,
            # ring_atom_frac,
            # largest_ring_size,
        # "HAC",
        # "C Frac",
        # "Ring Atom Frac",
        # "Largest Ring Size",
        # #, "rainbow", "rainbow", "rainbow", "Blues"],