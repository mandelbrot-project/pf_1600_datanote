# Scripts for the PF1600 project

This repository gathers scripts used for the paper _"Open and re-usable annotated mass spectrometry dataset of a chemodiverse collection of 1,600 plant extracts."_  

#### -- Project Status: [Active]

## Intro 

The purpose of the project is to explore the chemodiversity of a large collection of plant extracts (>15,000 extracts) in the frame of a collaboration between Pierre Fabre Laboratoires and the University of Geneva. To establish and benchmark the data acquisition and data analysis methods, a subset of 1,600 plants was selected.
Their exploration makes the object of the current datanote.

### Partners

#### Pierre-Fabre Laboratories

- https://www.pierre-fabre.com/en

#### Phytochemistry & Bioactive Natural Products 

- School of Pharmaceutical Sciences, University of Geneva, CMU, Rue Michel-Servet 1, CH-1211 Geneva 4, Switzerland
- https://ispso.unige.ch/phytochimie/index.php

#### COMMONS Lab

- Departement of Biology, University of Fribourg, 1700 Fribourg, Switzerland
- https://www.unifr.ch/bio/en/groups/allard/


### Methods Used
* High-Performance Liquid Chromatography
* High-Resolution Mass Spectrometry
* Computational Mass Spectrometry (molecular networking, metabolite annotation etc.)
* Data Visualization


### Technologies
* R 
* Python

## Project Abstract 
> 
> As privileged structures, natural products often display potent biological activities. However, the discovery of novel bioactive scaffolds is often hampered by the chemical complexity of the biological matrices they are found in. Large natural extracts collections are thus extremely valuable for their chemical novelty potential but also complicated to exploit in the frame of drug-discovery projects. In the end, it is the pure chemical substances that are desired for structural determination purposes and bioactivity evaluation. Researchers interested in the exploration of large and chemodiverse extracts collections should thus establish strategies aiming to efficiently tackle such chemical complexity and access these structures. Establishing carefully crafted digital layers documenting the spectral and chemical complexity as well as bioactivity results of natural products extracts collections can help to prioritize time-consuming but mandatory isolation efforts. In this note, we report the results of our initial exploration of a collection of 1,600 plant extracts in the frame of a drug discovery effort. After describing the taxonomic coverage of this collection, we present the results of its liquid chromatography high-resolution mass spectrometric profiling and the exploitation of these profiles using computational solutions. The resulting annotated mass spectral dataset and associated chemical and taxonomic metadata are made available to the community and data reuse cases are proposed. We are currently continuing our exploration of this plant extracts collection for drug-discovery purposes (notably looking for novel anti-trypanosomatids, anti-infective and prometabolic compounds) and eco-metabolomics insights. We believe that such a dataset can be exploited and reused by researchers interested in computational natural products exploration.
> 



## Data availability  

Raw Data are available in the the MassIVE repository [MSV000087728](https://doi.org/doi:10.25345/C59J97).

## Getting Started

1. Clone this repo (for help see this [tutorial](https://help.github.com/articles/cloning-a-repository/)).

```
git clone https://github.com/mandelbrot-project/pf_1600_datanote.git
cd pf_1600_datanote
```

2. Install and activate the Conda environment from the environment.yml file. The TMAP package is only available on MacOS and Linux. If you are using Windows, you can use WSL.

```
conda env create --file environment.yml
conda activate pf_datanote
```

## Phylogenetic coverage calculation

Here are the instruction to calculate the phylogenetic coverage of the PF1600 dataset.
The resulting data is used to generate Figure 1. and the associated barplots

1. Set parameters in the [phylo_config.yaml](https://github.com/mandelbrot-project/pf_1600_datanote/blob/e23e573e011c498eeb9664337cd1e74229c133d5/config/phylo_config.yaml)

2. Run the [phylogenetic_coverage.R](https://github.com/mandelbrot-project/pf_1600_datanote/blob/e23e573e011c498eeb9664337cd1e74229c133d5/src/phylogenetic_coverage.R) script. Make sure to run it from the repository directory level.

```
Rscript src/phylogenetic_coverage.R
```

## Spectral TMAP establishment    
Here are the instruction to calculate the TMAP clustering of the features' MS/MS spectra of the PF1600 dataset.
The resulting data is used to generate Figure 2.
  
Run the [spectral_tmap.py](https://github.com/mandelbrot-project/pf_1600_datanote/blob/e23e573e011c498eeb9664337cd1e74229c133d5/src/spectral_tmap.py) script.
```
python src/spectral_tmap.py
```

## Structural TMAP establishment    
Here are the instruction to calculate the TMAP clustering of the annotated compounds from the PF1600 dataset.
The resulting data is used to generate Figure 3.
  
Run the [structural_tmap.py](https://github.com/mandelbrot-project/pf_1600_datanote/blob/e23e573e011c498eeb9664337cd1e74229c133d5/src/structural_tmap.py) script.
```
python src/structural_tmap.py
```

## Featured Notebooks/Analysis/Deliverables
* [Notebook/Markdown/Slide Deck Title](link)
* [Notebook/Markdown/Slide DeckTitle](link)
* [Blog Post](link)


## Contributing Members

|Name     |  Slack Handle   | 
|---------|-----------------|
|[Pierre-Marie Allard](https://github.com/oolonek)| @johnDoe        |
|[Arnaud Gaudry](https://github.com/ArnaudGaudry) |     @janeDoe    |
|[Adriano Rutz](https://github.com/Adafede) |     @janeDoe    |
|[Emmanuel Defossez](https://github.com/Edefossez) |     @janeDoe    |


## Contact

* If you have a problem to run the scripts please file an issue here https://github.com/mandelbrot-project/pf_1600_datanote/issues
* Feel free to directly contribute to the code !
