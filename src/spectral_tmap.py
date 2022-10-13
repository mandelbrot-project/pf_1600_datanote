import pandas as pd
import pickle
import os
import requests
import subprocess
import shlex
import zipfile
from matchms.filtering import add_losses
from matchms.filtering import normalize_intensities
from matchms.filtering import require_minimum_number_of_peaks
from matchms.filtering import select_by_relative_intensity

from matchms.importing import load_from_mgf
from spec2vec import SpectrumDocument

from collections import Counter
import numpy as np
import tmap as tm
from faerun import Faerun

from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors

#########################################################
###  PART 0: Defining paths and download input files  ###
#########################################################

# These lines allows to make sure that we are placed at the repo directory level 
from pathlib import Path

p = Path(__file__).parents[1]
os.chdir(p)

tmp_path = 'tmp/'

# input files 

dataset_annotations_filename = 'PF_full_datanote_spectral_match_results_repond_flat.tsv'
dataset_annotations_path = tmp_path + dataset_annotations_filename

metadata_filename = 'plant_extract_dataset_metadata.tsv'
metadata_path = tmp_path + metadata_filename

feature_table_filename = '210302_feature_list_filtered_strong_aligned_cosine03_quant_without_QCs_and_Blanks.csv'
feature_table_path = tmp_path + feature_table_filename

spectra_filename = '210302_feature_list_filtered_strong_aligned_cosine03_without_QCs_and_Blanks.mgf'
spectra_path = tmp_path + spectra_filename

component_index_filename = 'gnps_job/clusterinfosummarygroup_attributes_withIDs_withcomponentID/5babd364facf4f46aa8a2ae8e3b9c6bb.clustersummary'
component_index_path = tmp_path + component_index_filename

# tmap outputs
set_coordinates_path = tmp_path + "pf_set_spectral_coord.dat"

tmap_filename = 'pf1600_spectral_tmap_pos'

# outputs folder

output_folder = 'docs/'
tmap_output_folder = 'docs/data/outputs/spectral/tmap'

# do you want do download all the files (true) or do you have them locally (false)

download_option = False


if download_option == True:

    # downloading required inputs
    print('download data')

    print('Downloading ISDB annotations ...')
    # ISDB annotations
    URL = "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000087728/updates/2022-04-28_pmallard_b0e0f70c/other/PF_full_datanote/PF_full_datanote_spectral_match_results_repond_flat.tsv"
    response = requests.get(URL)
    open(dataset_annotations_path, "wb").write(response.content)

    print('Downloading sample Metadata ...')
    # Metadata
    URL = "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000087728/updates/2021-12-21_pmallard_cde57c76/metadata/211221_metadata_vgf.tsv"
    response = requests.get(URL)
    open(metadata_path, "wb").write(response.content)

    print('Downloading sample Spectra ...')
    # Spectra
    URL = "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000087728/updates/2022-09-14_pmallard_c2fb482f/other/210302_feature_list_filtered_strong_aligned_cosine03_without_QCs_and_Blanks.mgf"
    response = requests.get(URL)
    open(spectra_path, "wb").write(response.content)

    print('Downloading sample Feature table ...')
    # Feature table
    URL = "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000087728/updates/2022-09-14_pmallard_c2fb482f/other/210302_feature_list_filtered_strong_aligned_cosine03_quant_without_QCs_and_Blanks.csv"
    response = requests.get(URL)
    open(feature_table_path, "wb").write(response.content)


    print('Downloading Dataset Molecular Network ...')
    # Componentindex
    path_to_file = tmp_path + 'gnps_job.zip'
    path_to_folder = tmp_path + 'gnps_job/'
    job_url_zip = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResult?task="+"3197f70bed224f9ba6f59f62906839e9"+"&view=download_cytoscape_data"
    cmd = 'curl -d "" ' + job_url_zip + ' -o '+ path_to_file + ' --create-dirs'
    subprocess.call(shlex.split(cmd))
    with zipfile.ZipFile(path_to_file, 'r') as zip_ref:
        zip_ref.extractall(path_to_folder)
    os.remove(path_to_file)


#########################################################
###               PART 1: IMPORT DATA                 ###
#########################################################
print('import data')

# Spectra processing function
def apply_my_filters(s):
    """This is how one would typically design a desired pre- and post-
    processing pipeline."""
    s = normalize_intensities(s)
    s = select_by_relative_intensity(s, intensity_from = 0.01, intensity_to = 1.0)
    s = add_losses(s, loss_mz_from=10.0, loss_mz_to=400.0)
    s = require_minimum_number_of_peaks(s, n_required=5)
    return s

# Load samples metadata
optional_meta = pd.read_csv(metadata_path, sep = '\t', usecols=[
    'sample_id', 'organism_phylum', 'organism_class', 'organism_order', 'organism_family',\
    'organism_genus', 'organism_species'
])

optional_meta.dropna(inplace=True)

# Load feature table
quant_table = pd.read_csv(feature_table_path)
quant_table.columns = quant_table.columns.str.replace(".mzXML", "")

# save RT and m/z for later use
quant_table_rt_mz_index = quant_table[['row ID', 'cluster_index', 'row retention time', 'row m/z']].set_index('cluster_index', drop=True)

''' Get the taxon with the maximal MEAN intensity for each feature '''

# quant_table.set_index('cluster_index', drop=True, inplace=True)
# samples_col = [col for col in quant_table.columns if 'VGF' in col]
# quant_table = quant_table[samples_col].transpose()
# taxonomy = ['organism_phylum', 'organism_class', 'organism_order', 'organism_family', 'organism_genus', 'organism_species']
# max_taxon = {}
# for col in taxonomy:
#     merged = quant_table.merge(optional_meta[['sample_id', col]], left_index=True, right_on= 'sample_id', how='left')#.drop('sample_id', inplace=1, axis = 1)
#     max_taxon[col] = merged.groupby(col).mean().idxmax()
# max_taxon_df = pd.DataFrame(max_taxon).merge(quant_table_rt_mz_index, left_index=True, right_index=True)

''' Or get the taxon with the maximal intensity for each feature '''

samples_col = [col for col in quant_table.columns if 'VGF' in col]
quant_table['max'] = quant_table[samples_col].idxmax(axis=1)
max_taxon_df = quant_table[['cluster_index', 'row ID','row retention time', 'row m/z', 'max']].merge(optional_meta, left_on = 'max', right_on = 'sample_id', how='left')

''' Load component index from GNPS job '''
component_index = pd.read_csv(
    component_index_path,
    sep = '\t', usecols=['cluster index', 'componentindex', 'Smiles']
    )

max_taxon_df = max_taxon_df.merge(component_index, left_on='cluster_index', right_on='cluster index')
max_taxon_df = max_taxon_df.rename(columns={'Smiles': 'gnps_smiles'}).drop('cluster index', axis=1)
max_taxon_df['has_gnps_smiles'] = np.where(pd.isnull(max_taxon_df['gnps_smiles']), 'No', 'Yes')

''' Load ISDB annotations '''
df_annotations = pd.read_csv(dataset_annotations_path, usecols=['feature_id', 'structure_smiles', 'score_taxo', 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'], low_memory=False, sep="\t")
df_annotations = df_annotations[df_annotations['score_taxo'] > 0] # Keep only reweighted annotations
df_annotations.drop_duplicates(subset = 'feature_id', inplace=True)


max_taxon_df = max_taxon_df.merge(df_annotations, left_on = 'cluster_index', right_on='feature_id', how = 'left').drop('feature_id', axis=1)

# Load Spectra
spectrums = [apply_my_filters(s) for s in load_from_mgf(spectra_path)]
spectrums = [s for s in spectrums if s is not None]
reference_documents = [SpectrumDocument(s, n_decimals=2) for s in spectrums]
metadata_df = pd.DataFrame(s.metadata for s in spectrums)

### converting the type to int
metadata_df.feature_id = metadata_df.scans.astype(int)

### merging with the additional metadata
metadata_df = pd.merge(metadata_df[['feature_id']], max_taxon_df, left_on = "feature_id", right_on='row ID', how="left")
metadata_df['has_isdb_smiles'] = np.where(pd.isnull(metadata_df['structure_smiles']), 'No', 'Yes')
metadata_df.structure_smiles = metadata_df.structure_smiles.fillna('O1CCOC1c1c(C#CC(C)(C)C)cc(c(C#CC(C)(C)C)c1)C#Cc1cc(C#CCCC)cc(C#CCCC)c1')
metadata_df.fillna('None',inplace=True)
metadata_df.componentindex = metadata_df.componentindex.astype(str)
metadata_df.score_taxo = metadata_df.score_taxo.astype(str)

del(max_taxon_df, df_annotations, quant_table)

# Create a list of all our documents-spectra 
texts = []
for doc in reference_documents:
    texts.append(doc.words)

#########################################################
###             PART 2: SPECTRAL T-MAP                ###
#########################################################
print('build tmap')

# (Dirty) trick to pass from a counter object to a list
ctr = Counter()
for text in texts:
    ctr.update(text)
ctr = list(set(ctr))

all_words = {}
for i, key in enumerate(ctr):
    all_words[key] = i

############ Un-weighted TMAP generation ###############

enc = tm.Minhash()
lf = tm.LSHForest()

fingerprints = []
for text in texts:
    fingerprint = []
    for t in text:
        if t in all_words:
            fingerprint.append(all_words[t])
    fingerprints.append(tm.VectorUint(fingerprint))

lf.batch_add(enc.batch_from_sparse_binary_array(fingerprints))
lf.index()

config = tm.LayoutConfiguration()
# config.k = 50 #50
# config.k = 10 #50
config.node_size = 1 / 40 # 1/26
# config.mmm_repeats = 2 #2
# config.sl_extra_scaling_steps = 5 #5
# config.sl_scaling_type = tm.RelativeToDesiredLength #RelativeToAvgLength

x, y, s, t, _ = tm.layout_from_lsh_forest(lf, config=config)


# Save for further use
cord_path = tmp_path + 'pf_set_spectral_coord.dat'
x = list(x)
y = list(y)
s = list(s)
t = list(t)
pickle.dump(
    (x, y, s, t), open(cord_path, "wb+"), protocol=pickle.HIGHEST_PROTOCOL
)

del (lf)

# Load precomputed coordinates
# x, y, s, t = pickle.load(open("pf_set_spectral_coord.dat",
#                               "rb"))


#### Tmap categories to map

###### Fuction to restrain the classes to the top N most frequent for better visualisation in tmap
def Top_N_classes(N, input_data, input_labels):
    """Keep only the top N classes for the input classes and replace the following by 'Other'"""
    top_input = [i for i, _ in Counter(input_data).most_common(N)]

    output_labels = [(7, "Other")]
    input_map = [7] * len(input_data)
    value = 1
    for i, name in input_labels:     
        if i in top_input:
            v = value
            if v == 7:
                v = 0
            output_labels.append((v, name))
            input_map[i] = v
            value += 1
    output_data = [input_map[val] for _, val in enumerate(input_data)]

    return output_data, output_labels
 ######

species_labels, species_data = Faerun.create_categories(metadata_df["organism_species"])
species_data, species_labels = Top_N_classes(9, species_data, species_labels)

genus_labels, genus_data = Faerun.create_categories(metadata_df["organism_genus"])
genus_data, genus_labels = Top_N_classes(9, genus_data, genus_labels)

family_labels, family_data = Faerun.create_categories(metadata_df["organism_family"])
family_data, family_labels = Top_N_classes(9, family_data, family_labels)

order_labels, order_data = Faerun.create_categories(metadata_df["organism_order"])
order_data, order_labels = Top_N_classes(9, order_data, order_labels)

class_labels, class_data = Faerun.create_categories(metadata_df["organism_class"])
#class_data, class_labels = Top_N_classes(19, class_data, class_labels)

phylum_labels, phylum_data = Faerun.create_categories(metadata_df["organism_phylum"])
#phylum_data, phylum_labels = Top_N_classes(19, phylum_data, phylum_labels)

ci_labels, ci_data = Faerun.create_categories(metadata_df["componentindex"])
ci_data, ci_labels = Top_N_classes(9, ci_data, ci_labels)

npc_pathway_labels, npc_pathway_data = Faerun.create_categories(metadata_df["structure_taxonomy_npclassifier_01pathway"])
npc_pathway_data, npc_pathway_labels = Top_N_classes(9, npc_pathway_data, npc_pathway_labels)

npc_superclass_labels, npc_superclass_data = Faerun.create_categories(metadata_df["structure_taxonomy_npclassifier_02superclass"])
npc_superclass_data, npc_superclass_labels = Top_N_classes(9, npc_superclass_data, npc_superclass_labels)

npc_class_labels, npc_class_data = Faerun.create_categories(metadata_df["structure_taxonomy_npclassifier_03class"])
npc_class_data, npc_class_labels = Top_N_classes(9, npc_class_data, npc_class_labels)

score_taxo_labels, score_taxo_data = Faerun.create_categories(metadata_df["score_taxo"])

has_smiles_labels, has_smiles_data = Faerun.create_categories(metadata_df["has_isdb_smiles"])

# Selected classes to plot
def selected_classes(list_selection, input_data, input_labels):
    """Keep only the top N classes for the input classes and replace the following by 'Other'"""
    output_labels = []
    selected_i = []
    key_dict = {}
    j = 1
    for i, name in input_labels:     
        if name in list_selection:
            output_labels.append((j, name))
            selected_i.append(i)
            key_dict[i] = j
            j += 1
        else:
            output_labels.append((0, "Other"))
    output_labels = list(set(output_labels))        
    output_data = [key_dict[item] if item in selected_i else 0 for item in input_data]
    return output_data, output_labels

temp_labels, temp_data = Faerun.create_categories(metadata_df["organism_family"])
sel_fam_data, sel_fam_labels = selected_classes(['Apocynaceae', 'Annonaceae', 'Meliaceae', 'Asteraceae', 'Fabaceae'], temp_data, temp_labels)

temp_labels, temp_data = Faerun.create_categories(metadata_df["structure_taxonomy_npclassifier_02superclass"])
sel_sc_data, sel_sc_labels = selected_classes(['Tryptophan alkaloids', 'Linear polyketides', 'Triterpenoids', 'Sesquiterpenoids', 'Flavonoids'], temp_data, temp_labels)

temp_labels, temp_data = Faerun.create_categories(metadata_df["structure_taxonomy_npclassifier_03class"])
sel_c_data, sel_c_labels = selected_classes(['Pregnane steroids', 'Acetogenins', 'Limonoids'], temp_data, temp_labels)


size_superclass = []
for value in npc_superclass_data:
    if value in [7, 5]:
        size_superclass.append(0.7)
    else:
        size_superclass.append(1)
    
size_class = []
for value in npc_class_data:
    if value in [7, 6]:
        size_class.append(0.7)
    else:
        size_class.append(1)
    
size_ci = []
for value in ci_data:
    if value in [7, 1]:
        size_ci.append(0.5)
    else:
        size_ci.append(1)
        
size_sel_sc = []
for value in sel_sc_data:
    if value == 0:
        size_sel_sc.append(0.5)
    else:
        size_sel_sc.append(1)
 
size_sel_c = []
for value in sel_c_data:
    if value == 0:
        size_sel_c.append(0.5)
    else:
        size_sel_c.append(1)
                
size = [1]*len(size_ci)

# Colors stuff
tab_10 = cm.get_cmap("tab10")
colors = [tab_10(i) for i in range(tab_10.N)]
colors[7] = '#e6e6e6' #(0.230,0.230,0.230)
#colors[9] = '#e6e6e6' #(0.230,0.230,0.230)
cmap2 = mcolors.ListedColormap(colors)

tab_10 = cm.get_cmap("tab10")
colors = [tab_10(i) for i in range(tab_10.N)]
colors[7] = '#e6e6e6' #(0.230,0.230,0.230)
colors[1] = '#e6e6e6' #(0.230,0.230,0.230)
cmap3 = mcolors.ListedColormap(colors)

tab_10 = cm.get_cmap("tab10")
colors = [tab_10(i) for i in range(tab_10.N)]
colors[7] = '#e6e6e6' #(0.230,0.230,0.230)
colors[5] = '#e6e6e6' #(0.230,0.230,0.230)
cmap4 = mcolors.ListedColormap(colors)

tab_10 = cm.get_cmap("tab10")
colors = [tab_10(i) for i in range(tab_10.N)]
colors[7] = '#e6e6e6' #(0.230,0.230,0.230)
colors[6] = '#e6e6e6' #(0.230,0.230,0.230)
cmap5= mcolors.ListedColormap(colors)

tab_10 = cm.get_cmap("tab10")
colors = [tab_10(i) for i in range(tab_10.N)]
colors[7] = '#e6e6e6' #(0.230,0.230,0.230)
colors[5] = '#e6e6e6' #(0.230,0.230,0.230)
cmap6= mcolors.ListedColormap(colors)

cmap7 = mcolors.ListedColormap(['#e6e6e6', "#023047"])

cmap_sel_fam = mcolors.ListedColormap(['#e6e6e6', '#ef476f', '#ffd166', '#06d6a0', '#118ab2', '#073b4c'])

cmap_sel_sc = mcolors.ListedColormap(['#e6e6e6', '#457b9d', '#43aa8b', '#786150', '#e8ac65', '#ffd60a'])

cmap_sel_c = mcolors.ListedColormap(['#e6e6e6', '#2a9d8f', '#FF8C61', '#bc6c25'])

# Creating the smiles label
metadata_df["label_smiles"] = (
    metadata_df["structure_smiles"]
    + '__'
    + 'Feature ID is: '
    + '__'
    + metadata_df.feature_id.astype(str)
    + '__'
    + 'componentindex is: '
    + '__'
    + metadata_df.componentindex.astype(str)
    + '__'
    + "</a>"
    )

# Plotting the spectral t-map
print('plot tmap')

os.chdir(os.path.join(os.getcwd(), tmap_output_folder))

f = Faerun(view="front", coords=False, clear_color='#FFFFFF')
f.add_scatter(
        "attributes",
        {
        "x": x,
        "y": y,
         "c": [
            species_data,
            genus_data,
            family_data,
            order_data,
            class_data,
            phylum_data,
            ci_data,
            npc_pathway_data,
            npc_superclass_data,
            npc_class_data,
            score_taxo_data,
            has_smiles_data,
            sel_fam_data,
            sel_sc_data,
            sel_c_data,
            metadata_df['row m/z'],
            metadata_df['row retention time']
        ],
         "s": [
             size,
             size,
             size,
             size,
             size,
             size,
             size_ci,
             size,
             size_superclass,
             size_class,
             size,
             size,
             size,
             size_sel_sc, 
             size_sel_c,
         ],
        "labels": metadata_df["label_smiles"],
    },
    shader="smoothCircle",
    point_scale= 4.0,
    max_point_size= 10,
    categorical= [
        True,
        True,
        True,
        True,
        True,
        True,
        True, 
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        False,
        False
        ],
    colormap= [
        cmap2,
        cmap2,
        cmap2,
        cmap2,
        cmap2,
        cmap2,
        cmap3,
        cmap4,
        cmap6,
        cmap5,
        cmap2,
        cmap7,
        cmap_sel_fam,
        cmap_sel_sc,
        cmap_sel_c,
        "rainbow",
        "rainbow"
        ],
    legend_labels=[
        species_labels,
        genus_labels,
        family_labels,
        order_labels,
        class_labels,
        phylum_labels,
        ci_labels,
        npc_pathway_labels,
        npc_superclass_labels,
        npc_class_labels,
        score_taxo_labels,
        has_smiles_labels,
        sel_fam_labels,
        sel_sc_labels,
        sel_c_labels
        ],
    series_title= [
        "Species max", 
        "Genus max",
        "Family max",
        "Order max",
        "Class max",
        "Phylum max",
        "Componentindex",
        "NPC pathway",
        'NPC superclass',
        'NPC class',
        'Score taxo',
        'Annotation status',
        "Selected families",
        "Selected NPC superclasses",
        "Selected NPC classes",        
        "Parent mass",
        "RT"
        ],
    has_legend= True,
)
f.add_tree("tree",
           {"from": s, "to": t},
           point_helper="attributes",
           color='#a89e9e'
           )
f.plot(tmap_filename,
       template="smiles"
)
