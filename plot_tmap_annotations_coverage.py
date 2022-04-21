import pickle
import sys
from collections import Counter
from time import sleep

import matplotlib.colors as mcolors
from matplotlib import cm
from matplotlib import pyplot as plt
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
    usecols=['organism_taxonomy_02kingdom', 'structure_smiles', 'structure_smiles_2D', 'structure_taxonomy_npclassifier_01pathway',
    'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'],
    sep=",")

plant_smiles = []
for i, row in df_meta.iterrows():
    if row['organism_taxonomy_02kingdom'] == 'Archaeplastida':
        plant_smiles.append(row['structure_smiles_2D'])
plant_smiles = set(plant_smiles)

# Keep only 1 occurence of each 2d structure
df_meta = df_meta.drop_duplicates(subset=['structure_smiles_2D'])

########################################################
###           PART 2: LOAD annotations               ###
########################################################

# Load annotations
df_annotations = pd.read_csv('PF_full_datanote_spectral_match_results_repond_flat.tsv', usecols=['structure_smiles', 'score_taxo', 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'], low_memory=False, sep="\t")

df_annotations = df_annotations[df_annotations['score_taxo'] > 0] # Keep only reweighted annotations

anisomeric_smiles = []
inchikey = []
for i, row in df_annotations.iterrows():
    mol = AllChem.MolFromSmiles(row["structure_smiles"])
    if mol != None:
        aniso_smiles = AllChem.MolToSmiles(mol, isomericSmiles = False)
        ik = AllChem.MolToInchiKey(mol)
        anisomeric_smiles.append(aniso_smiles)
        inchikey.append(ik)
    else:
        anisomeric_smiles.append(np.nan)
        inchikey.append(np.nan)

df_annotations['anisomeric_smiles'] = anisomeric_smiles
df_annotations['inchikey'] = inchikey
df_annotations['inchikey_2D'] = df_annotations['inchikey'].str[:14]

df_annotations = df_annotations.drop_duplicates(subset=['inchikey_2D'])
df_annotations.drop(['structure_smiles', 'inchikey', 'score_taxo'], axis=1, inplace=True)

df_annotations.rename({'structure_taxonomy_npclassifier_01pathway':'pathway_annotation', 'structure_taxonomy_npclassifier_02superclass':'superclass_annotation', 'structure_taxonomy_npclassifier_03class':'class_annotation'}, axis=1, inplace=True)

print('Number of unique annotations (2D structures): ' + str(len(df_annotations)))

# Merge both tables:

inchikey = []
for i, row in df_meta.iterrows():
    mol = AllChem.MolFromSmiles(row["structure_smiles_2D"])
    if mol != None:
        ik = AllChem.MolToInchiKey(mol)
        inchikey.append(ik)
    else:
        inchikey.append(np.nan)

df_meta['inchikey'] = inchikey
df_meta.dropna(subset=['inchikey'], inplace=True)
df_meta['inchikey_2D'] = df_meta['inchikey'].str[:14]
df_meta = df_meta.drop_duplicates(subset=['inchikey_2D'])

# Check is structure reported in plant
df_meta['is_in_plant'] = df_meta['structure_smiles_2D'].isin(plant_smiles)

df_merged = pd.merge(df_meta, df_annotations, left_on = 'inchikey_2D', right_on='inchikey_2D', how='outer')

df_merged['is_annotated'] = np.where(df_merged['anisomeric_smiles'].isnull(), 'no', 'yes')
df_merged['structure_smiles_2D'].fillna(df_merged['anisomeric_smiles'], inplace=True)
df_merged['structure_taxonomy_npclassifier_01pathway'].fillna(df_merged['pathway_annotation'], inplace=True)
df_merged['structure_taxonomy_npclassifier_02superclass'].fillna(df_merged['superclass_annotation'], inplace=True)
df_merged['structure_taxonomy_npclassifier_03class'].fillna(df_merged['class_annotation'], inplace=True)
df_merged['is_in_plant'].fillna(True, inplace=True)

df_merged.drop(['anisomeric_smiles', 'pathway_annotation', 'superclass_annotation', 'class_annotation'], axis=1, inplace=True)

values = {
    'structure_taxonomy_npclassifier_01pathway': 'Unknown',
    'structure_taxonomy_npclassifier_02superclass': 'Unknown',
    'structure_taxonomy_npclassifier_03class': 'Unknown',
          }
df_merged.fillna(value=values, inplace=True)

# Special tricks to plot selected classes
chemical_classes_to_plot = [
    'Cholestane steroids', 'Cyclic peptides', 'Flavones', 'Flavonols',
    'Germacrane sesquiterpenoids', 'Guaiane sesquiterpenoids', 'Limonoids', 'Oleanane triterpenoids'
    ]

df_merged['chemical_classes_to_plot'] = np.where(
    df_merged['structure_taxonomy_npclassifier_03class'].isin(chemical_classes_to_plot),
    df_merged['structure_taxonomy_npclassifier_03class'], 'Others')

# Generating class for plotting
dic_categories = {
    "structure_taxonomy_npclassifier_01pathway": {'Ncat': 9},
    "structure_taxonomy_npclassifier_02superclass": {'Ncat': 9},
    "structure_taxonomy_npclassifier_03class": {'Ncat': 9},
    "chemical_classes_to_plot": {'Ncat': 9},
    "is_annotated": {'Ncat': 0},
    "is_in_plant": {'Ncat': 0}
}

for dic in dic_categories:
    labels, data = Faerun.create_categories(df_merged[str(dic)])
    dic_categories[dic]['data'], dic_categories[dic]['labels'] = Top_N_classes(dic_categories[dic]['Ncat'], data, labels)

df_merged.reset_index(inplace=True)
df_merged['index'] = df_merged['index'].astype(str)

# Generate a labels column
df_merged["labels"] = (
        df_merged["structure_smiles_2D"]
        + '__'
        + df_merged["structure_taxonomy_npclassifier_03class"]
        + '__'
        + df_merged["index"]
        + '__'
        + "</a>"
)

cmap = mcolors.ListedColormap(['#e6e6e6', "#023047"])

tab_10 = cm.get_cmap("tab10")
colors = [tab_10(i) for i in range(tab_10.N)]
colors[7] = '#e6e6e6' #(0.230,0.230,0.230)
colors[9] = '#e6e6e6' #(0.230,0.230,0.230)
cmap2 = mcolors.ListedColormap(colors)

cmap3 = mcolors.ListedColormap(['#EC9F05', "#04a777"])

list_hex = ['#FF8C61', '#CE6A85', '#04a777', '#5fa8d3', '#1b4965', '#985277', '#8f2600', '#e6e6e6', '#5C374C', '#e6e6e6']
cmap4 = mcolors.ListedColormap(list_hex)


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
lf.store("210803_pf_set_annotation_vs_lotus.dat")
with open("210803_pf_set_annotation_vs_lotus_attributes.pickle", "wb+") as f:
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
    (x, y, s, t), open("210803_pf_set_annotation_vs_lotus_coords.dat", "wb+"), protocol=pickle.HIGHEST_PROTOCOL
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

x, y, s, t = pickle.load(open("210803_pf_set_annotation_vs_lotus_coords.dat",
                              "rb"))

hac, c_frac, ring_atom_frac, largest_ring_size = pickle.load(
    open("210803_pf_set_annotation_vs_lotus_attributes.pickle", "rb")
)

########################################################
# COMMON PART: PLOT THE TMAP
########################################################
c_frak_ranked = ss.rankdata(np.array(c_frac) / max(c_frac)) / len(c_frac)

# Plotting function
f = Faerun(view="front", coords=False, clear_color='#ffffff')
f.add_scatter(
    "attributes",
    {
        "x": x,
        "y": y,
        "c": [
            dic_categories['structure_taxonomy_npclassifier_01pathway']['data'],
            dic_categories['structure_taxonomy_npclassifier_02superclass']['data'],
            dic_categories['structure_taxonomy_npclassifier_03class']['data'],
            dic_categories['chemical_classes_to_plot']['data'],
            dic_categories['is_annotated']['data'],
            dic_categories['is_in_plant']['data'],
            hac,
            c_frak_ranked,
            ring_atom_frac,
            largest_ring_size,
        ],
        "labels": df_merged["labels"],
    },
    shader="smoothCircle",
    point_scale=4.0,
    max_point_size=10,
    legend_labels=[
        dic_categories['structure_taxonomy_npclassifier_01pathway']['labels'],
        dic_categories['structure_taxonomy_npclassifier_02superclass']['labels'],
        dic_categories['structure_taxonomy_npclassifier_03class']['labels'],
        dic_categories['chemical_classes_to_plot']['labels'],
        dic_categories['is_annotated']['labels'],
        dic_categories['is_in_plant']['labels']
        ],
    categorical=[True, True, True, True, True, True, False, False, False, False],
    colormap=[cmap2, cmap2, cmap2,  cmap4, cmap, cmap3, "rainbow", "rainbow", "rainbow", "Blues"],
    series_title=["chemo_np_kingdom", "chemo_np_superclass", "chemo_np_class", "Selected classes", "Is annotated", "Found in plants", "HAC", "C Frac", "Ring Atom Frac", "Largest Ring Size",],
    has_legend=True
)
f.add_tree("annotation_tree", {"from": s, "to": t}, point_helper="attributes", color='#e6e6e6')
f.plot('210803_annotation_vs_lotusdnp_tmap', template="smiles")


### Bar Chart

import plotly.express as px
import plotly.graph_objects as go

chemical_classes_to_plot = [
    'Cholestane steroids', 'Cyclic peptides', 'Flavones',
    'Flavonols', 'Germacrane sesquiterpenoids', 'Guaiane sesquiterpenoids',
    'Limonoids',
    'Oleanane triterpenoids']

colors_chemical_classes =[
    '#FF8C61', '#985277', '#8f2600', '#5C374C', '#5fa8d3', '#1b4965', '#CE6A85', '#04a777'
    ]

df_for_plot = df_merged[df_merged['structure_taxonomy_npclassifier_03class'].isin(chemical_classes_to_plot)]

df_for_plot = df_for_plot[['structure_taxonomy_npclassifier_03class', 'is_annotated', 'structure_smiles_2D']]

df_for_plot = df_for_plot[df_for_plot['structure_taxonomy_npclassifier_03class'] != 'Unknown']
df_for_plot['is_annotated'] = np.where(df_for_plot['is_annotated'] == 'yes', 1, 0)

df_for_plot = df_for_plot.groupby('structure_taxonomy_npclassifier_03class').agg(
    {'structure_smiles_2D': 'count',
     'is_annotated': 'sum',
    })

df_for_plot.sort_values(['is_annotated'], ascending=False, inplace=True)
df_for_plot['structure_smiles_2D'] = df_for_plot['structure_smiles_2D'] - df_for_plot['is_annotated']

# Version 1
fig = go.Figure(go.Bar(x=df_for_plot.index, y=df_for_plot.is_annotated, name='Annotated', text=df_for_plot.is_annotated, marker_color='#023047'))
fig.add_trace(go.Bar(x=df_for_plot.index, y=df_for_plot.structure_smiles_2D, name='Reported', text=df_for_plot.structure_smiles_2D, marker_color='#e6e6e6'))

fig.update_traces(marker_line_color=colors_chemical_classes,
                  marker_line_width=10)

fig.update_layout(
    barmode='stack', template='simple_white', width=1500, height=1000,
    font=dict(
        size=40,
    ),
    xaxis_title="NPClassifier class",
    yaxis_title="Count",)

fig.update_xaxes(tickfont = dict(size=40))
fig.show()
fig.write_image("barchart_2.jpg")


# Version 2

fig = go.Figure(data=[
    go.Bar(name='Total', x=df_for_plot.index, y=df_for_plot.is_annotated,marker_color=colors_chemical_classes),
    go.Bar(name='Annotation', x=df_for_plot.index, y=df_for_plot.structure_smiles_2D,marker_color=colors_chemical_classes, opacity =0.4)
])
# Change the bar mode
fig.update_layout(
    barmode='stack', template='simple_white', width=1500, height=1000,
    xaxis_title="Selected chemical classes",
    yaxis_title="Count",
    font=dict(
        size=40,
    ))
fig.show()
fig.write_image("barchart_3.jpg")
