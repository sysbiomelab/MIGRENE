%% Box-B| calculation of reactionScore, constraining the model based on diet, MAGMA generation
%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk
%% start
% get path to where the MIGRENE Toolbox is located
MIGDIR = fileparts(which('MIGRENE_pipeline'));
% provide the path to bacterial species (MSP) gene info and taxonomy
% information.
CATDIR=[MIGDIR filesep 'data'];
taxo=[CATDIR filesep 'Taxonomy.xlsx'];
speciesGeneInfo=[CATDIR filesep 'MSPgeneProfile.txt'];
% provide the path to microbiomeGEM (generated in
% IntegrationCatalogToModel.m) and bibliome data.
MATDIR=[MIGDIR filesep 'mat'];
MODEL=[MATDIR filesep 'microbiomeGEM.mat'];
BIBLIOME=[MATDIR filesep 'bibliome.mat'];
% define a directory to save microbiomeGEM
SAVEDIR=[MIGDIR filesep 'saveDir'];
% number of cores specified for parallelization. it can be a positive
% integer or a range specified as a 2-element vector of integers
numWorkers=4

% for some functions such as FBA simulation, you need to install
% cobra toolbox.
%please visit https://opencobra.github.io/cobratoolbox/latest/installation.html
% for installing the cobra toolbox. after installation, initiate COBRA 
initCobraToolbox()

%% 
%load microbiomeGEM model as a reference GEM.
load(MODEL)

% you can constrain the model before generating your bacterial GEMs
% for example here, the reference model is constraned based on the high fiber 
% animal based model.
% we already generated 5 diets. the following function constrain the model.

dietNumber=2; % 1:high Fibre Plant Based, 2:high Fibre omnivore, 3:high Protein Plant based
              % 4:high protein omnivore, 5:UK average.
[microbiomeGEM]=DietConstrain(microbiomeGEM,dietNumber);
% acetate and lactate are added as popular carbon sources for bacteria
microbiomeGEM.lb(find(strcmp(microbiomeGEM.rxns, 'Ex_Acetate')))= -2.597426442;
microbiomeGEM.lb(find(strcmp(microbiomeGEM.rxns, 'Ex_L-Lactate')))= -0.074862638;

% generating new diet based on USDA dataset
newDiet=false
% if you have a diet and want to apply for modeling, you can use the
% USDAcreatingDiet function as below. your compounds also could be matched
% with USDA IDs by using searchfood function. here we provided an example
% for high protein plant based diet. the diet will be converted to mmol/gDW
% for modeling
if newDiet 
    %Example: high protein plant based diet
    load([MATDIR filesep 'DietInput.mat']);
    %DietInput.mat includes: food_id_item: the iD of the food items that are
    %present in a high protein plant based diet, food_grams_item: the amount
    %(grams) of each food ID item
    [micronutrients_diet_mmol, macronutrients_diet]= USDAcreatingDiet(food_id_item,food_grams_item)
end

%filling a structure for species gene info.
MSPinfo = struct();
% First, the species(MSP) gene info is read from the text file. it provides
% the genes of 10 Bacteroides in a matrix.
T = readtable(speciesGeneInfo);
MSPinfo.genes=table2cell(T(:,1));
MSPinfo.msp=T.Properties.VariableNames;
MSPinfo.msp=transpose(MSPinfo.msp(1,2:end));
MSPinfo.expression=table2array(T(1:end ,2:end));

% second, generate reaction state (absent/present reaction) for bacteria
% and prune the genes based on bacterial (MSP) profile for each species.
for h=1:length(MSPinfo.msp)
    if ~exist([SAVEDIR filesep MSPinfo.msp{h} '.mat']) 	
            metagenomeData=struct();
            metagenomeData.gene=MSPinfo.genes;
            metagenomeData.value=MSPinfo.expression(:,h);
            [RxnState, modelforMSP] = MetagenomeToReactions(microbiomeGEM, metagenomeData);
            % save RxnState and modelforMSP in a mat file entitled the
            % corresponding bacterium (MSP)
            save([SAVEDIR filesep MSPinfo.msp{h} '.mat'],'RxnState','modelforMSP');
    else
        disp(['the reaction state and modelforMSP for ' MSPinfo.msp{h} ' are already generated. see your directory'])
    end 
end

% collect MSP (bacterial) Information
% get path to where RxnState and modelforMSP for each bacterium (MSP) were
% saved
RXNDIR=SAVEDIR;
% the following function generates MSPInformation, a structure includes the
% following fields:
% taxoLevel, the taxonomy names. taxoInfo, taxonomy information for each
% species. taxoGroup: taxonomy group for bacteria. rxns, the reaction name
% in reference model. bacteria, list of MSP IDs. BacteriaNames, list of
% species names. RxnStateAll, the reaction state (absent/present) for each
% bacteria

[MSPInformation]= GenerateMSPInformation(taxo,RXNDIR,microbiomeGEM)

% convet reaction state to reaction score and calculate the threshold for
% gapfilling
for h=1:length(MSPInformation.bacteria)
    % add the msp name into MSPInformation 
    MSPInformation.species=MSPInformation.bacteria{h};
    [reactionScore, threshold] = MetaGenomicsReactionScore(MSPInformation);
%   adds new variables to the corresponding MAT-file
    save([RXNDIR filesep MSPInformation.bacteria{h} '.mat'],'reactionScore','threshold','-append');
end 
%% collecting the bibliome data and the constraining the model

% if you have any bibliome data about phenotypic features of the bacteria
% that you are making GEM, provide it as a structure with four fields:
% "bacteria" a cell array listing the name of the bacteria; "rxn" list the
% name of the reactions having bibliome; "value" a matrix of numbers: zero
% means no information, 1 means consumed, 2 means produced, -1 not-consumed
% and -2 means not-produced by the corresponding bacteria. "aerobeInfo" a
% cell array provides the info that the bacteria require oxygen for growth
% or not so specefiy with "aerobe" or "anaerobe" or "facultative".
% if you do not provide the information, then the models are just generated
% and tuned based on the reactionScore and threshold.
load(BIBLIOME)

%collect all exchange and transport reactions
indexEx=strfind(microbiomeGEM.rxns,'Ex');
IndexEx = find(not(cellfun('isempty',indexEx)));
indexTr=strfind(microbiomeGEM.rxns,'t_');
IndexTr = find(not(cellfun('isempty',indexTr)));
indexOfTrEx=union(IndexEx,IndexTr);
modelseed='true' %if the reference model is based on KBase or modelSEED  
for h=1:length(MSPInformation.bacteria)
    load([RXNDIR filesep MSPInformation.bacteria{h} '.mat'])
    %to keep all the exchange reactions and transport reactions in the models
    %and prune it later, before generating the species-specific GEMs (MAGMA)
    %the score of the reaction changed as 1
    reactionScore(indexOfTrEx,:)=1;
    % add a field, named species, to the biblome for finding the
    % corresponding info in the structure including all the bacteria
    if exist('bibliome')
        bibliome.species=MSPInformation.bacteria{h};
        [contextSpecificModel] = contextSpecificModelGenertion(modelforMSP,reactionScore,threshold,bibliome);
    else
        [contextSpecificModel] = contextSpecificModelGenertion(modelforMSP,reactionScore,threshold);
    end
     MSPInformation.species=MSPInformation.bacteria{h};
    [contextSpecificModel, modelInfo] = contextSpecificModelTune(contextSpecificModel,MSPInformation,reactionScore,threshold,modelseed);
    %adds new variables to the corresponding MAT-file
    save([RXNDIR filesep MSPInformation.bacteria{h} '.mat'],'contextSpecificModel','modelInfo','-append');
end

% collect modelInfo Tables and save as excel file. regarding the level
% and the pecentage of gapfilling, you can decide which bacterial models
% were generated based on enough gene info and discard the poor models 
modelInfoFinal=table()
for h=1:length(MSPInformation.bacteria)
    load([RXNDIR filesep MSPInformation.bacteria{h} '.mat'],'modelInfo')
    modelInfo=[table(MSPInformation.bacteria(h)) modelInfo];
    modelInfoFinal=[modelInfoFinal;modelInfo];
end

% Write Data to Excel Spreadsheets
filename=[RXNDIR filesep 'modelInfoFinal.xlsx']
writetable(modelInfoFinal,filename,'Sheet',1,'Range','A1')

% done, congrats. go to PersonalizedMicrobiomeMetabolism.m in tutorials
% directory for integration of bacterial profile into the models to
% investigate the metabolism of microbiome at personalized level.