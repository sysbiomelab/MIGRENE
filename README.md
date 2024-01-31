# MIGRENE toolbox,
## Description 
MIGRENE toolbox is an integrated pipeline for Microbial and personalized GEM (Genome-scale metabolic model), REactiobiome, and community NEtwork modeling. It enables the generation of species and community-level models from any reference gene catalogs and metagenome species to be applied to personalized microbiome studies. Using the toolbox, GEMs could be generated based on the gut microbial gene catalogs and metagenomic species pan-genomes (MSPs). This toolbox also contains functions for performing community modelling using GEMs, determining reaction abundance and richness and reaction set enrichment (RSE), and reactobiome that describes an aggregate of the metabolic repertoires of an individual gut microbiome, or the biochemical state of the microbiome.

## Download and installation
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
## Tutorials
<a href="https://github.com/sysbiomelab/MIGRENE/wiki/generation-of-microbiome-GEM"> Generation of microbiome GEM using a generic metabolic model and a microbiome catalog</a>: This tutorial shows how to integrate a bacterial gene catalog 
into the metabolic model to generate a microbiome reference genome-scale metabolic model (GEM).

<a href="https://github.com/sysbiomelab/MIGRENE/wiki/generation-of-microbiome-GEM"> Generation of Bacterial (species-specific) GEM </a>: This tutorial shows how to calculate the reaction score and threshold for bacteria, to constrain the model and to generate species-specific bacterial GEMs. 

<a href="https://github.com/sysbiomelab/MIGRENE/wiki/generation-of-microbiome-GEM"> Generation of Personalized Microbiome Metabolism </a>: It shows how to calculate reaction richness, reactobiome, reaction abundance, community models and iRSE (individualized reaction set enrichment)

## Integration of a gene catalog into a metabolic model.
### Data usage
* `<catalog>`: is a txt file containing gene names and KO (KEGG orthology) such as [SubSet_hs_10_4_igc2_annot.txt](data/SubSet_hs_10_4_igc2_annot.txt)
* `<mapping file>`: (optional) a txt file contains the mapping information for KO to KEGG reaction ID.
* `<Metabolic model>`: a mat file containig a metabolic models whether COBRA or RAVEN format such as [RefMetabolicModel.mat](mat/RefMetabolicModel.mat) 
* `<database_type>`: has to be `n` if the user is using a nucleotide database or `a` if the user is using an amino acids database
### functions
* [checkCatalog](Functions/checkCatalog.m): check the `<catalog>` to make sure it is ready for integration.
* [checkCatalog](Functions/convertCatalogAnnotation.m): Convert KO annotations to KEGG reaction IDs in the  `<catalog>`. If no mapping file `<mapping file>` is provided, the latest information is automatically downloaded from the KEGG API. the output is a `<converted catalog>`.

microbiomeGEMgeneration
* 
*  
### output
## GEM generation

## Reactobiome generation

## Community Model Generation


## Contact
gholamreza.bidhkori@kcl.ac.uk,
gbidkhori@gmail.com,
# Citation
in preparation, 2024

