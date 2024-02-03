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

<a href="https://github.com/sysbiomelab/LiverCirrhosis_MS"> Reactobiome and reaction richness for Liver Cirrhosis gut microbiome samples </a>


## Integration of a gene catalog into a metabolic model.
### Data usage
* `<catalog>`: is a txt file containing gene names and KO (KEGG orthology) such as [SubSet_hs_10_4_igc2_annot.txt](data/SubSet_hs_10_4_igc2_annot.txt)
* `<mapping file>`: (optional) a txt file contains the mapping information for KO to KEGG reaction ID.
* `<Metabolic_model>`: a mat file containig a metabolic models whether COBRA or RAVEN format such as [RefMetabolicModel.mat](mat/RefMetabolicModel.mat) 
### functions
* [checkCatalog](Functions/checkCatalog.m): check the `<catalog>` to make sure it is ready for integration.
* [convertCatalogAnnotation](Functions/convertCatalogAnnotation.m): Convert KO annotations to KEGG reaction IDs in the  `<catalog>`. If no mapping file `<mapping file>` is provided, the latest information is automatically downloaded from the KEGG API. the output is a `<converted catalog>`.
* [microbiomeGEMgeneration](Functions/microbiomeGEMgeneration.m): integrate the `<converted catalog>` into the `<Metabolic_model>` to generate a microbiome `<reference_GEM>`


## Generation of Bacterial (species‐specific) GEM
### Data usage
* `<reference_GEM>`: the genome scale metabolic Model with COBRA or RAVEN format such as the `<reference_GEM>` produced above by [microbiomeGEMgeneration](Functions/microbiomeGEMgeneration.m)
* `<bacterial_info>`: an structure from a binary matrix containing gene-level data for bacterial species (exanmple [MSPgeneProfile.txt](data/MSPgeneProfile.txt)).
    
  ```
  T = readtable('MSPgeneProfile.txt');
  bacterial_info = struct();
  bacterial_info.genes=table2cell(T(:,1));
  bacterial_info.msp=T.Properties.VariableNames;
  bacterial_info.msp=transpose(MSPinfo.msp(1,2:end));
  bacterial_info.expression=table2array(T(1:end ,2:end));
  ```
* `<taxonomy>`: an Excel file contains taxonomy classification info i.e. Kingdom, Phylum, Class, Order, Family, Genus, Species (exanmple [MSPgeneProfile.txt](data/MSPgeneProfile.txt)). 
* `<Bibliome_Data>`: (optional, example [here](mat/bibliome.mat)) any bibliome data about phenotypic features of the bacteria can be provided as a structure with four fields: "bacteria" is a cell array listing the names of the bacteria. "rxn" lists the name of the reactions having bibliome. "value" is a matrix of numbers: zero means no information, 1 means consumed, 2 means produced, -1 means not consumed and -2 means not produced by the corresponding bacteria. "aerobeInfo" a cell array provides the info that the bacteria require oxygen for growth or not, specifying with "aerobe", "anaerobe" or "facultative".
  
### functions
* [DietConstrain](Functions/DietConstrain.m): (optional) this function constrains `<reference_GEM>` based on the provided diet `<diet_number>` (1 to 5). Five diets have been provided by the toolbox: 1: high Fibre Plant Based, 2: high Fibre omnivore, 3: high Protein Plant based, 4: high protein omnivore, 5:UK average. Set the number of the diet for constraining the model.

* [MetagenomeToReactions](Functions/MetagenomeToReactions.m): This function needs `<reference_GEM>`and `<bacterial_info>` as inputs and must be seperately run for each bacterial species in the `<bacterial_info.msp>` using a loop, as below:
```
for h=1:length(MSPinfo.msp) 	
            metagenomeData=struct();
            metagenomeData.gene=bacterial_info.genes;
            metagenomeData.value=bacterial_info.expression(:,h);
            [Reaction_State, bacterial_model] = MetagenomeToReactions(microbiomeGEM, metagenomeData);
            save(['save\to\directory\' bacterial_info.msp{h} '.mat'],'Reaction_State','bacterial_model');
end
```
`<bacterial_model>` and `<Reaction_State>` for each species are the output that must be saved in the output directry in the same `mat` file entitled the bacterial name.
`<bacterial_model>` is a model with the bacterial genes and gene rules and `<Reaction_State>` is a vector showing the state of the reaction (zero or one) for the bacterial species.
* [GenerateMSPInformation](Functions/GenerateMSPInformation.m): this function generates  `<bacterial_Information>`, a structure that includes the following fields: taxoLevel, the taxonomy names. taxoInfo, taxonomy information for each species. taxoGroup: taxonomy group for bacteria. rxns, the reaction name in the reference model. bacteria, list of MSP IDs. BacteriaNames, list of species names. RxnStateAll, the reaction state (absent/present) for each bacteria. the input is the address to the directory including saved `Reaction_State>` and `bacterial_model>` for each bacterium (MSP), `<reference_GEM>` and address to `taxonomy>` file.
* [MetaGenomicsReactionScore](Functions/MetaGenomicsReactionScore.m): This function utilize `<bacterial_Information>` to converts reaction states to reaction scores (`<reaction_Score>`) and calculate a threshold (`<threshold>`) for each bacterial species. `<reaction_Score>`, `<threshold>` must be added to the `mat` file including `<bacterial_model>` and `<Reaction_State>`
* [contextSpecificModelTune](Functions/contextSpecificModelTune.m): `<bacterial_model>`,`<reaction_Score>`, `<threshold>` and `<Bibliome_Data>`to genrate context specefic species genome scale metabolic model (`<bacterial_GEM>`) as the output. [contextSpecificModelTune](Functions/contextSpecificModelTune.m) function tunes `<bacterial_GEM>` and also provides the level and the details of gap filling.

## Reactobiome and reaction richness Generation

### Data usage

* `<modelList>`: an array to provide the names of bacterial GEM, such as row names of this [file]
* `<sampleName>`: an array to provide the names of samples or subject, such as column names of this [file]
* `<bacterial_abundance>`: the bacterial abundance matrix showing the abundance of bacteria in `<modelList>` in each subject in `<sampleName>`  note: for example go to [link1](https://github.com/sysbiomelab/LiverCirrhosis_MS) or [link2](https://github.com/sysbiomelab/MIGRENE/wiki/Generation-of-Personalized-Microbiome-Metabolism). 
* `<PathToModels>`: Provides the path where the [bacterial GEMs](https://github.com/sysbiomelab/LiverCirrhosis_MS/blob/main/GEMmodels.zip) are saved in `PathToModels.path` field and the name of the model assigned in the .mat files in `PathToModels.name`. note: if you use the provided [bacterial GEMs](https://github.com/sysbiomelab/LiverCirrhosis_MS/blob/main/GEMmodels.zip), unzip it and then ```PathToModels.name='model'```.

### Functuions
* [RxnRichnessGenerator](Functions/RxnRichnessGenerator.m): uses `<modelList>`,`<PathToModels>`,`<abundance>` and `<sampleName>` and generates reaction richness for all the subjects regarding the provided GEMs.
* [ReactobiomeGenerator](Functions/ReactobiomeGenerator.m): utelizes `<modelList>`,`<PathToModels>`,`<abundance>` and `<sampleName>` and generates reactobiome for all the subjects regarding the provided GEMs.

## Contact
gholamreza.bidhkori@kcl.ac.uk,
gbidkhori@gmail.com,
# Citation
in preparation, 2024

