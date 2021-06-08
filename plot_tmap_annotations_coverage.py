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
df_meta = pd.read_csv('../../data/processed/210523_frozen_metadata.csv.gz', sep=",")

# Fill NaN with 'Unknown'
values = {'organism_taxonomy_02kingdom': 'Not attributed (Bacteria and Algae)', 'organism_taxonomy_03phylum': 'Unknown',
          'organism_taxonomy_04class': 'Unknown',
          'organism_taxonomy_05order': 'Unknown', 'organism_taxonomy_06family': 'Unknown',
          'organism_taxonomy_08genus': 'Unknown', 'organism_taxonomy_09species': 'Unknown'}
df_meta.fillna(value=values, inplace=True)

# Keep only 1 occurence of structure - biosource pair
df_meta = df_meta.drop_duplicates(
    subset=['structure_wikidata', 'organism_taxonomy_02kingdom', 'organism_taxonomy_03phylum',
            'organism_taxonomy_04class',
            'organism_taxonomy_05order', 'organism_taxonomy_06family', 'organism_taxonomy_08genus',
            'organism_taxonomy_09species'])

# Group by structure and return for each taxonomical level the most frequent taxa, its occurence and the number of unique taxa
df_gb = df_meta.groupby('structure_wikidata').agg(
    {'structure_inchi': 'first',
     'structure_inchikey': 'first',
     'structure_smiles': 'first',
     'structure_smiles_2D': 'first',
     'structure_taxonomy_npclassifier_01pathway': 'first',
     'structure_taxonomy_npclassifier_02superclass': 'first',
     'structure_taxonomy_npclassifier_03class': 'first',
     'structure_taxonomy_classyfire_01kingdom': 'first',
     'structure_taxonomy_classyfire_02superclass': 'first',
     'structure_taxonomy_classyfire_03class': 'first',
     'structure_taxonomy_classyfire_04directparent': 'first',
    #  'organism_taxonomy_02kingdom': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1), 'count',
    #                                  'nunique', lambda x: Counter(x)],
    #  'organism_taxonomy_03phylum': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1),
    #                                 'nunique', lambda x: Counter(x)],
    #  'organism_taxonomy_04class': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1),
    #                                'nunique', lambda x: Counter(x)],
    #  'organism_taxonomy_05order': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1),
    #                                'nunique', lambda x: Counter(x)],
    #  'organism_taxonomy_06family': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1),
    #                                 'nunique', lambda x: Counter(x)],
    #  'organism_taxonomy_08genus': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1),
    #                                'nunique', lambda x: Counter(x)],
    #  'organism_taxonomy_09species': [lambda x: pd.Series.mode(x), lambda x: x.value_counts().head(1),
    #                                  'nunique', lambda x: Counter(x)],
     })

del (df_meta)

# df_gb.columns = ['_'.join(col).strip() for col in df_gb.columns.values]

# df_gb.rename(columns={"organism_taxonomy_02kingdom_count": "biosource_count"}, inplace=True)

# taxo_cols = ['organism_taxonomy_02kingdom', 'organism_taxonomy_03phylum', 'organism_taxonomy_04class',
#              'organism_taxonomy_05order', 'organism_taxonomy_06family', 'organism_taxonomy_08genus',
#              'organism_taxonomy_09species']

# bins = [0, 2, 6, 21, np.inf]
# names = ['1 unique biosource', '2-5 unique biosources', '6-20 unique biosources', '>20 unique biosources']

# for col in taxo_cols:
#     df_gb[col + '_<lambda_0>'] = df_gb[col + '_<lambda_0>'].astype(str)
#     df_gb[col + '_specificity'] = df_gb[col + '_<lambda_1>'] / df_gb['biosource_count']
#     df_gb.drop([col + '_<lambda_1>'], axis=1, inplace=True)
#     df_gb[col + '_nunique_cat'] = pd.cut(df_gb[col + '_nunique'], bins, labels=names, right=False)

# names = ['1 biosource', '2-5 biosources', '6-20 biosources', '>20 biosources']
# df_gb['biosource_count_cat'] = pd.cut(df_gb.biosource_count, bins, labels=names, right=False)

# df_gb.reset_index(inplace=True)

# df_gb["id_wikidata"] = df_gb['structure_wikidata'].str[31:]

values = {'structure_taxonomy_npclassifier_01pathway_first': 'Unknown',
          'structure_taxonomy_npclassifier_02superclass_first': 'Unknown',
          'structure_taxonomy_npclassifier_03class_first': 'Unknown',
          'structure_taxonomy_classyfire_01kingdom_first': 'Unknown',
          'structure_taxonomy_classyfire_02superclass_first': 'Unknown',
          'structure_taxonomy_classyfire_03class_first': 'Unknown',
          'structure_taxonomy_classyfire_04directparent_first': 'Unknown'}
df_gb.fillna(value=values, inplace=True)

# # Generating mean specificity score by taxonomical level and class, with a minimal score for plotting, set it to 0 not to use it. 

# dic_mean_specificity = {
#     'organism_taxonomy_02kingdom': {'min': 0},
#     'organism_taxonomy_06family': {'min': 0},
#     'organism_taxonomy_08genus': {'min': 0},
# }
# structure_taxo = ['structure_taxonomy_npclassifier_01pathway_first',
#                   'structure_taxonomy_npclassifier_02superclass_first',
#                   'structure_taxonomy_npclassifier_03class_first']

# for dic in dic_mean_specificity:
#     for level in structure_taxo:
#         print(str(dic)+'_<lambda_2>')
#         df_gb['counter_sum'] = df_gb.groupby(by = [level])[str(dic)+'_<lambda_2>'].transform('sum')
#         df_gb['counter_sum_values'] = df_gb['counter_sum'].apply(lambda x: sum(x.values()))
#         df_gb['counter_sum'] = df_gb['counter_sum'].apply(lambda x: x.most_common(1)[0])
#         df_gb['specificity_taxon'] = df_gb['counter_sum'].apply(lambda x: x[0])
#         df_gb['specificity_score'] = df_gb['counter_sum'].apply(lambda x: x[1])
#         df_gb['specificity_score'] = df_gb['specificity_score']/df_gb['counter_sum_values']

#         dic_mean_specificity[dic][level+'_score'] = df_gb['specificity_score']
#         dic_mean_specificity[dic][level+'_taxon'] = df_gb['specificity_taxon']

#         df_gb.drop(['counter_sum', 'counter_sum_values','specificity_taxon', 'specificity_score'], axis=1, inplace=True)

########################################################
###           PART 2: LOAD anotations                ###
########################################################

# Load lotus
df_annotations = pd.read_csv('../../data/processed/210523_frozen_metadata.csv.gz', sep=",")



# Generating class for plotting
dic_categories = {
    'structure_taxonomy_npclassifier_01pathway_first': {'Ncat': 19},
    "structure_taxonomy_npclassifier_02superclass_first": {'Ncat': 19},
    "structure_taxonomy_npclassifier_03class_first": {'Ncat': 19},
    "structure_taxonomy_classyfire_01kingdom_first": {'Ncat': 0},
    "structure_taxonomy_classyfire_02superclass_first": {'Ncat': 19},
    "structure_taxonomy_classyfire_03class_first": {'Ncat': 19},
    "structure_taxonomy_classyfire_04directparent_first": {'Ncat': 19},
    "organism_taxonomy_02kingdom_<lambda_0>": {'Ncat': 19},
    "organism_taxonomy_03phylum_<lambda_0>": {'Ncat': 19},
    "organism_taxonomy_04class_<lambda_0>": {'Ncat': 19},
    "organism_taxonomy_05order_<lambda_0>": {'Ncat': 19},
    "organism_taxonomy_06family_<lambda_0>": {'Ncat': 19},
    "organism_taxonomy_08genus_<lambda_0>": {'Ncat': 19},
    "organism_taxonomy_09species_<lambda_0>": {'Ncat': 19},
    "biosource_count_cat": {'Ncat': 0},
    "organism_taxonomy_02kingdom_nunique_cat": {'Ncat': 0},
    "organism_taxonomy_03phylum_nunique_cat": {'Ncat': 0},
    "organism_taxonomy_04class_nunique_cat": {'Ncat': 0},
    "organism_taxonomy_05order_nunique_cat": {'Ncat': 0},
    "organism_taxonomy_06family_nunique_cat": {'Ncat': 0},
    "organism_taxonomy_08genus_nunique_cat": {'Ncat': 0},
    "organism_taxonomy_09species_nunique_cat": {'Ncat': 0}
}

for dic in dic_categories:
    labels, data = Faerun.create_categories(df_gb[str(dic)])
    dic_categories[dic]['data'], dic_categories[dic]['labels'] = Top_N_classes(dic_categories[dic]['Ncat'], data,
                                                                               labels)

# Publication examples, but you can use it as you want!
labels, data = Faerun.create_categories(df_gb['organism_taxonomy_06family_<lambda_0>'])
simaroubaceae_data, simaroubaceae_labels = keep_only_given_class(['Simaroubaceae'], data, labels)

labels, data = Faerun.create_categories(df_gb['structure_taxonomy_npclassifier_03class_first'])
NPclass_data, NPclass_labels = keep_only_given_class([
    'Oleanane triterpenoids', 'Germacrane sesquiterpenoids', 'Carotenoids (C40, β-β)',
    'Flavonols', 'Lanostane, Tirucallane and Euphane triterpenoids', 'Cyclic peptides',
    'Guaiane sesquiterpenoids', 'Flavones', 'Colensane and Clerodane diterpenoids',
    'Quassinoids'], data, labels)

# count_simaroubaceae = df_gb.organism_taxonomy_06family_join.str.count("Simaroubaceae")
# simaroubaceae_specificity = count_simaroubaceae / df_gb['biosource_count']

# Generating colormaps for plotting
cmap = mcolors.ListedColormap(["gainsboro", "peachpuff", "salmon", "tomato"])
cmap2 = mcolors.ListedColormap(["gainsboro", "tomato"])
cmap3 = mcolors.ListedColormap(["lightgray", "firebrick"])
cmap4 = mcolors.ListedColormap(
    ["gainsboro", "#83C644", "#04AC6E", "#00ff99", "#008ECA", "#008C83", "#D61E1E", "#006A7A", "#00B2C4", "#2F4858",
     "#EAD400"])

# Generate a labels column
df_gb["labels"] = (
        df_gb["structure_smiles_first"]
        + '__<a target="_blank" href="'
        + df_gb["structure_wikidata"]
        + '">'
        + df_gb["id_wikidata"]
        + '__'
        + df_gb['structure_inchikey_first']
        + '__'
        + df_gb["structure_taxonomy_npclassifier_03class_first"]
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

n_iter = len(df_gb)
for i, row in df_gb.iterrows():
    j = (i + 1) / n_iter
    sys.stdout.write('\r')
    sys.stdout.write(f"[{'=' * int(50 * j):{50}s}] {round((100 * j), 2)}%")
    sys.stdout.flush()
    sleep(0.05)

    mol = AllChem.MolFromSmiles(row["structure_smiles_2D_first"])
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
    # fps.append(tm.VectorUint(enc.encode_mol(mol)))
    # fps.append(tm.VectorUint(MAP4.calculate(mol)))

fps = MAP4.calculate_many(fps)
lf.batch_add(fps)
lf.index()

# Store lsh forest and structure metadata
# lf.store("../../data/interim/tmap/210506_lotus_2D_map4.dat")
# with open("../../data/interim/tmap/210506_lotus_2D_map4.pickle", "wb+") as f:
#     pickle.dump(
#         (hac, c_frac, ring_atom_frac, largest_ring_size),
#         f,
#         protocol=pickle.HIGHEST_PROTOCOL,
#     )

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
#     (x, y, s, t), open("../../data/interim/tmap/210506_coords_lotus_2D_map4.dat", "wb+"), protocol=pickle.HIGHEST_PROTOCOL
# )

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

x, y, s, t = pickle.load(open("../../data/interim/tmap/210506_coords_lotus_2D_map4.dat",
                              "rb"))  # Version "coords_210312_lotus.dat" contains 270'336 resulting from 210223_frozen_metadata groupby(structure_wikidata)

hac, c_frac, ring_atom_frac, largest_ring_size = pickle.load(
    open("../../data/interim/tmap/210506_lotus_2D_map4.pickle", "rb")
)

########################################################
# COMMON PART: PLOT THE TMAP
########################################################
c_frak_ranked = ss.rankdata(np.array(c_frac) / max(c_frac)) / len(c_frac)

# Plotting function
f = Faerun(view="front", coords=False, clear_color='#ffffff')
f.add_scatter(
    "lotus",
    {
        "x": x,
        "y": y,
        "c": [
            dic_categories['structure_taxonomy_npclassifier_01pathway_first']['data'],
            dic_categories['structure_taxonomy_npclassifier_02superclass_first']['data'],
            dic_categories['structure_taxonomy_npclassifier_03class_first']['data'],
            dic_categories['structure_taxonomy_classyfire_01kingdom_first']['data'],
            dic_categories['structure_taxonomy_classyfire_02superclass_first']['data'],
            dic_categories['structure_taxonomy_classyfire_03class_first']['data'],
            dic_categories['structure_taxonomy_classyfire_04directparent_first']['data'],
            dic_categories['organism_taxonomy_02kingdom_<lambda_0>']['data'],
            dic_categories['organism_taxonomy_03phylum_<lambda_0>']['data'],
            dic_categories['organism_taxonomy_04class_<lambda_0>']['data'],
            dic_categories['organism_taxonomy_05order_<lambda_0>']['data'],
            dic_categories['organism_taxonomy_06family_<lambda_0>']['data'],
            dic_categories['organism_taxonomy_08genus_<lambda_0>']['data'],
            dic_categories['organism_taxonomy_09species_<lambda_0>']['data'],
            dic_categories['biosource_count_cat']['data'],
            dic_categories['organism_taxonomy_02kingdom_nunique_cat']['data'],
            dic_categories['organism_taxonomy_03phylum_nunique_cat']['data'],
            dic_categories['organism_taxonomy_04class_nunique_cat']['data'],
            dic_categories['organism_taxonomy_05order_nunique_cat']['data'],
            dic_categories['organism_taxonomy_06family_nunique_cat']['data'],
            dic_categories['organism_taxonomy_08genus_nunique_cat']['data'],
            dic_categories['organism_taxonomy_09species_nunique_cat']['data'],
            simaroubaceae_data,
            NPclass_data,
            dic_mean_specificity['organism_taxonomy_02kingdom'][
                'structure_taxonomy_npclassifier_01pathway_first_score'],
            dic_mean_specificity['organism_taxonomy_02kingdom'][
                'structure_taxonomy_npclassifier_02superclass_first_score'],
            dic_mean_specificity['organism_taxonomy_02kingdom'][
                'structure_taxonomy_npclassifier_03class_first_score'],
            dic_mean_specificity['organism_taxonomy_06family'][
                'structure_taxonomy_npclassifier_01pathway_first_score'],
            dic_mean_specificity['organism_taxonomy_06family'][
                'structure_taxonomy_npclassifier_02superclass_first_score'],
            dic_mean_specificity['organism_taxonomy_06family'][
                'structure_taxonomy_npclassifier_03class_first_score'],
            dic_mean_specificity['organism_taxonomy_08genus'][
                'structure_taxonomy_npclassifier_01pathway_first_score'],
            dic_mean_specificity['organism_taxonomy_08genus'][
                'structure_taxonomy_npclassifier_02superclass_first_score'],
            dic_mean_specificity['organism_taxonomy_08genus'][
                'structure_taxonomy_npclassifier_03class_first_score'],
            hac,
            c_frak_ranked,
            ring_atom_frac,
            largest_ring_size,
        ],
        "labels": df_gb["labels"],
    },
    shader="smoothCircle",
    point_scale=2.0,
    max_point_size=10,
    legend_labels=[
        dic_categories['structure_taxonomy_npclassifier_01pathway_first']['labels'],
        dic_categories['structure_taxonomy_npclassifier_02superclass_first']['labels'],
        dic_categories['structure_taxonomy_npclassifier_03class_first']['labels'],
        dic_categories['structure_taxonomy_classyfire_01kingdom_first']['labels'],
        dic_categories['structure_taxonomy_classyfire_02superclass_first']['labels'],
        dic_categories['structure_taxonomy_classyfire_03class_first']['labels'],
        dic_categories['structure_taxonomy_classyfire_04directparent_first']['labels'],
        dic_categories['organism_taxonomy_02kingdom_<lambda_0>']['labels'],
        dic_categories['organism_taxonomy_03phylum_<lambda_0>']['labels'],
        dic_categories['organism_taxonomy_04class_<lambda_0>']['labels'],
        dic_categories['organism_taxonomy_05order_<lambda_0>']['labels'],
        dic_categories['organism_taxonomy_06family_<lambda_0>']['labels'],
        dic_categories['organism_taxonomy_08genus_<lambda_0>']['labels'],
        dic_categories['organism_taxonomy_09species_<lambda_0>']['labels'],
        dic_categories['biosource_count_cat']['labels'],
        dic_categories['organism_taxonomy_02kingdom_nunique_cat']['labels'],
        dic_categories['organism_taxonomy_03phylum_nunique_cat']['labels'],
        dic_categories['organism_taxonomy_04class_nunique_cat']['labels'],
        dic_categories['organism_taxonomy_05order_nunique_cat']['labels'],
        dic_categories['organism_taxonomy_06family_nunique_cat']['labels'],
        dic_categories['organism_taxonomy_08genus_nunique_cat']['labels'],
        dic_categories['organism_taxonomy_09species_nunique_cat']['labels'],
        simaroubaceae_labels,
        NPclass_labels
    ],
    categorical=[True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                 True, True, True, True, True, True, True,
                 False, False, False, False, False, False, False, False, False, False, False, False],
    colormap=["tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20", "tab20",
              "tab20", "tab20", "tab20", cmap, cmap2, cmap, cmap, cmap, cmap, cmap, cmap, cmap3, cmap4,
              "inferno", "inferno", "inferno", "inferno", "inferno", "inferno", "inferno", "inferno",
              "inferno", "rainbow", "rainbow", "rainbow", "Blues"],
    series_title=[
        "chemo_np_classifier_pathway",
        "chemo_np_classifier_superclass",
        "chemo_np_classifier_class",
        "chemo_cf_kingdom",
        "chemo_cf_superclass",
        "chemo_cf_class",
        "chemo_cf_parent",
        "kingdom_bio",
        "phylum_bio",
        "class_bio",
        "order_bio",
        "family_bio",
        "genus_bio",
        "species_bio",
        "biosource_count",
        "kingdom_nunique",
        "phylum_nunique",
        "class_nunique",
        "order_nunique",
        "family_nunique",
        "genus_nunique",
        "species_nunique",
        "Simaroubaceae vs others",
        "Top8 classes, quassinoids and carotenoids",
        "Pathway kingdom specificity",
        "Superclass kingdom specificity",
        "Class kingdom specificity",
        "Pathway family specificity",
        "Superclass family specificity",
        "Class family specificity",
        "Pathway genus specificity",
        "Superclass genus specificity",
        "Class genus specificity",
        "HAC",
        "C Frac",
        "Ring Atom Frac",
        "Largest Ring Size",
    ],
    has_legend=True,
)
f.add_tree("lotus_tree", {"from": s, "to": t}, point_helper="lotus", color='#e6e6e6')
f.plot('../../res/html/210506_lotus_map4_2D', template="smiles")
