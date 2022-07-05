%% Box-C| generating personalized metabolic microbiome
%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk
%% start
% get path to where the MIGRENE Toolbox is located
MIGDIR = fileparts(which('MIGRENE_pipeline'));
% provide the path to bacterial species (MSP) gene info and bacterial
% abundance obtained from metagenomics analysis
CATDIR=[MIGDIR filesep 'data'];
ABUNDANCE=[CATDIR filesep 'BacterialAbundance.xlsx'];
PATHWAY=[CATDIR filesep 'pathways.xlsx'];
% provide the path to microbiomeGEM (generated in
% IntegrationCatalogToModel.m) and bibliome data.
MATDIR=[MIGDIR filesep 'mat'];
MODEL=[MATDIR filesep 'microbiomeGEM.mat'];
BIBLIOME=[MATDIR filesep 'bibliome.mat'];
% define a directory to save microbiomeGEM
SAVEDIR=[MIGDIR filesep 'saveDir'];
% number of cores specified for parallelization. it can be a positive
% integer or a range specified as a 2-element vector of integers
numWorkers=4;
% for some functions such as FBA simulation, you need to install
% cobra toolbox
initCobraToolbox()

%% 
%load microbiome (MSP) abundance profile. it could be metagenomics or 16s based 
[abundance,infoFile,~]=xlsread(ABUNDANCE);
%name of models
modelList = infoFile(2:end,1);
%name of samples
sampleName=infoFile(1,2:end);
%check the samples,
%remove the MSP name if the abundance of bacteria in all samples are zero
abundance=abundance(sum(abundance,2)~=0,:);
modelList=modelList(sum(abundance,2)~=0,:);
%remove the samples if the there is no bacterial abundance
sampleName=sampleName(sum(abundance,1)~=0);
abundance=abundance(:,sum(abundance,1)~=0);
% get the number of bacteria (bacterial richness) in each sample 
temp1=abundance;
temp1(find(temp1>0))=1;
BactrialRichness=table(sampleName',sum(temp1,1)');
BactrialRichness.Properties.VariableNames = {'sampleName','BactrialRichness'};

% give the path where the models are available and the name of model assgined in the .mat files
PathToModels.path=SAVEDIR;
PathToModels.name='contextSpecificModel';
% generate gut microbiome reaction composition (reaction richness) of all individuals 
richness= RxnRichnessGenerator(modelList,PathToModels,abundance,sampleName);

% generate reaction abundance for all individuals; the function generates both reaction abundance
% and relative reaction abundance
[reactionRelativeAbun, rxnAbunPerSample]= ReactionAbundanceGenerator(modelList,PathToModels,abundance,sampleName);
% generate reactobiome for all individuals
countPerFiveBacteria= ReactobiomeGenerator(modelList,PathToModels,abundance,sampleName);

% for the enrichment analysis, you need to prepare two files:
%1) a file includes pathway terms and the IDs that could be KO,EC,kegg
%reaction ID and etc. here we use the KEGG pathway terms with Kegg reaction
%ID
%provide the pathway profiles
[~,terms,~]=xlsread(PATHWAY);
terms=terms(2:end,:);
%2) a file for ID mapping between reaction names in the models and and 
% the IDs in the pathway file. here, we use the info in the reference model. 
load(MODEL)
IDmap=[microbiomeGEM.rxns microbiomeGEM.rxnRN];
Index = find(not(cellfun('isempty',IDmap(:,2))));
IDmap=IDmap(Index,:);
% if you assign the p-value, then coverage of non-significance terms is set
% as zero. if you dont define the p-value, it returns all. 
p_value=0.05;
[coverageRSE,pRSE]= pRSEGenerator(modelList,PathToModels,abundance,sampleName,IDmap,terms,p_value);
% coverage is a table shows the coverage of each pathway in samples
% pRSE is a table shows the p value of the pathways in samples
%%
%community modeling
% define the number of top abundant bacteria for community modeling 
%here we generate communities for top 10 bacteria 
top=5;
thre=[];
for i=1:size(abundance,2)
t1=sort(abundance(:,i),1,'descend');
thre(i,1)=t1(top,1);
abundance(find(abundance(:,i) < thre(i,1)),i)=0;
end
boxplot(thre)
median(thre)

%specify the metabolite ID and exchange reaction for biomass (optional)
biomass.EXrxn='Ex_Biomass';
biomass.mets='cpd11416ee[lu]';
% make a directory to save generated community models
if ~exist([SAVEDIR filesep 'community'],'dir')
mkdir([SAVEDIR filesep 'community']);
end
PathToSave=[SAVEDIR filesep 'community'];
% in report, "one" next to sample name shows that community model has been
% generated for the individuals in PathToSave directory
[report]= MakeCommunity(modelList,PathToModels,abundance,sampleName,PathToSave,biomass);
