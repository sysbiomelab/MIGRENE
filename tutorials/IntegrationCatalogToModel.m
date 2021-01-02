%% Box-A| generation of microbiome GEM using a generic metabolic model and a microbime catalog 
%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk
%% start
% get path to where the MIGRENE Toolbox is located
MIGDIR = fileparts(which('MIGRENE_pipeline'));
% provide the path to microbiome catalog .
CATDIR=[MIGDIR filesep 'data'];
% provide the path to reference metabolic model.
MATDIR=[MIGDIR filesep 'mat'];
% define a directory to save microbiomeGEM, here it will be saved in mat
% directory
SAVEDIR=[MIGDIR filesep 'mat'];
% number of cores specified for parallelization. it can be a positive
% integer or a range specified as a 2-element vector of integers
numWorkers=4

%% integration
% integration of bacterial gene catalog into metabolic model to generate a
% generic genome scale metabolic.
%First the cataloge data is read from the text file.
% here, a subset of updated gut catalog is used to run the pipeline and
% generate MAGMA. this small annotated cataloge includes genes
%of 10 Bacteroides. their taxonomy information is also provided for MAGMA
%generation
catalog=[CATDIR filesep 'SubSet_hs_10_4_igc2_annot.txt'];
T = readtable(catalog);
catalogData=table2cell(T) ;
%before using the catalog, the following function rearranges the catalog
%for the genes with more that KO linked, checks the structure, and provides
%the format compatible for downstream functions
[catalogFileChecked]= checkCatalog(catalogData,numWorkers);

%convert KO in the catalog to KEGG reaction IDs. if there is no ID mapping file,
%assign an empty cell array. then, it will automatically download the last
%updated information from KEGG API and saves in directory "data" where 
%the MIGRENE Toolbox is located.(make sure, you are connected to internet) 
mapping={};
inputFile=catalogFileChecked;
[catalogConverted]= convertCatalogAnnotation(inputFile,mapping,numWorkers);

% gene-protein-reaction (GPR) association is assigned by integrating the
% the catalog genes into metabolic model and generate generic genome scale
% metabolic model. it will be compatible with both COBRA and RAVEN
% toolboxes so you can use any functions provided by both toolboxes.

% first, the annotated metabolic model model is loaded. you can load your
% generic model.
load([MATDIR filesep 'RefMetabolicModel.mat']);
genericModel=model;
% if you have annotation file, please load or import it. otherwise, leave
% the annotationFile empty. then, the function automatically find all type
% of annotations in your model and download the corresponding info from
% KEGG API, (make sure, you are connected to internet). here a generic
% metabolic model is used. the reactions in the model are annotated by KO
% and kegg RN ID.
annotationFile={};

[microbiomeGEM]=microbiomeGEMgeneration(genericModel,catalogConverted,...
    annotationFile,numWorkers);

%save microbiomeGEM to a MAT-file
save([SAVEDIR filesep 'microbiomeGEM.mat'],'microbiomeGEM')
% done, congrats. go to MAGMAgeneration.m in tutorials directory