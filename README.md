# MIGRENE toolbox,
## Description 
MIGRENE toolbox is an integrated pipeline for Microbial and personalized GEM (Genome-scale metabolic model), REactiobiome, and community NEtwork modeling. It enables the generation of species and community-level models from any reference gene catalogs and metagenome species to be applied to personalized microbiome studies. Using the toolbox, GEMs could be generated based on the gut microbial gene catalogs and metagenomic species pan-genomes (MSPs). This toolbox also contains functions for performing community modelling using GEMs, determining reaction abundance and richness and reaction set enrichment (RSE), and reactobiome that describes an aggregate of the metabolic repertoires of an individual gut microbiome, or the biochemical state of the microbiome and profiling the reactobiome.

# Download and installation
1. Download this repository. You can clone the repository using:
```
git clone https://github.com/sysbiomelab/MIGRENE.git
```
Or you can download this repository as a <a href="https://codeload.github.com/sysbiomelab/MIGRENE/zip/refs/heads/master">compressed archive</a>.

2. Change to the folder MIGRENE/ and run from Matlab
```
addpath(genpath("MIGRENE"))
```
Or you can use the <a href="https://uk.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html">link</a> to learn how to set path in MATLAB to the directory.

# three tutorials are provided in the toolbox:
Three tutorials shows the steps that MIGRENE Toolbox automatically generate and simulate MAGMA
(MSP Associated Genome scale MetAbolic) models and personalized metabolic microbiome data
using Bacterial gene catalog, metagenome species (MSP) and metagenomic data integration.
note: if you have your GEMS and you need to create the personalized metabolic microbiome data
i.e. reaction richness, Microbiome, reaction abundance, community models and pRSE (personalized
reaction set enrichment), go to Box-c.

IntegrationCatalogToModel.m: Box-a|
integration of bacterial gene catalog into metabolic model to generate a microbiome reference genome
scale metabolic model (GEM).

MAGMAgeneration.m: Box-b|
calculation of reactionScore, constraining the model based on diet, species specific GEMs or MAGMA generation

PersonalizedMicrobiomeMetabolism.m: Box-c| generating personalized metabolic microbiome

# Contact
gholamreza.bidhkori@kcl.ac.uk,
gbidkhori@gmail.com,
# Citation
in preparation, 2023
# Acknowledgments

