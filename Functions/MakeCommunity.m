function [report]= MakeCommunity(modelList,PathToModels,abundance,sampleName,PathToSave,biomass)
%inputs:
%   modelList: 				list of model names.
%   PathToModels:			a structure includes the path where the models are available
%							and the name of model assigned in the .mat files
%	abundance:				matrix of microbiome (MSP) abundance profile
%	sampleName:				list of sample names
%	PathToSave:				a string showing the directory to save generated community models
%OPTIONAL INPUT
%	biomass:				a structure that specify the metabolite ID and exchange reaction name for biomass
%outputs:					
%   report:					shows generated community models for samples 

%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidhkori@kcl.ac.uk
if nargin<6
    biomass={};
end
index=[];
for h1 = 1:size(modelList,1)
	if exist([PathToModels.path filesep modelList{h1} '.mat'])
		index=[index;h1];
	end
end
modelList=modelList(index,:);
abundance=abundance(index,:);
%keep the bacteria with at leaset one nonzero observation
abundance=abundance(sum(abundance,2)~=0,:);
modelList=modelList(sum(abundance,2)~=0,:);
%keep the sample with at leaset one nonzero observation
abundance=abundance(:,sum(abundance,1)~=0);
sampleName=sampleName(sum(abundance,1)~=0);
exchangeMetabolites={};
report={};
for h1 = 1:size(modelList,1)
	load([PathToModels.path filesep modelList{h1}],PathToModels.name)
	model=eval(PathToModels.name);
    models{h1,1}=model;
    exchangeMets=GetExchangeMetabolite(model);
    exchangeMetabolites=vertcat(exchangeMetabolites,exchangeMets);
end
exchangeMetabolites=unique(exchangeMetabolites);
mets_art=exchangeMetabolites;
mets_art(:,2) = strrep(mets_art(:,1), '[e]', '[lu]');
mets_art(:,3) = strrep(mets_art(:,1), '[e]', '[fo]');
mets_art(:,4) = strrep(mets_art(:,1), '[e]', '[fe]');
reactions={};
for h1=1:length(models)
    temp=models{h1};
    indexEx=strfind(temp.rxns,'Ex_');
    IndexEx = find(not(cellfun('isempty',indexEx)));
    rxn=temp.rxns(IndexEx,1);
    rxn(:,2)=printRxnFormula(temp,rxn);
    reactions=vertcat(reactions,rxn);
end
[~,idx]=unique(strcat(reactions(:,1), 'rows'));
reactions1=reactions(idx,:)

Ex_art=[]
for h1=1:size(mets_art,1)
IndexC = strfind(reactions1(:,2),mets_art{h1,1});
Index = find(not(cellfun('isempty',IndexC)));
Ex_art{h1,1}=reactions1{Index,1}
end

Ex_art(:,2) = strrep(Ex_art(:,1), 'Ex_', 'FoEx_');
Ex_art(:,3) = strrep(Ex_art(:,1), 'Ex_', 'Fo_');
Ex_art(:,4) = strrep(Ex_art(:,1), 'Ex_', 'Fe_');
Ex_art(:,5) = strrep(Ex_art(:,1), 'Ex_', 'FeEx_');
[fakemodel, fakemodel1]=fakeModelGenerator(mets_art,Ex_art);
relativeAbundance=abundance;
if any(sum(abundance,1)~=1)
    SumOfAbun=sum(abundance,1)
    for i=1:length(sampleName)
       relativeAbundance(:,i)=abundance(:,i)/SumOfAbun(i);
    end
end
abundance=relativeAbundance;

for h1=1:numel(sampleName)
    abundanceS=abundance(:,h1)
    modelListN=modelList(find(abundanceS ~= 0));
    abundanceN=abundanceS(find(abundanceS ~= 0));
models={};
for h11 = 1:size(modelListN,1)
	load([PathToModels.path filesep modelListN{h11}],PathToModels.name)
	model=eval(PathToModels.name);
    minLB=min(model.lb);
    %remove the constraint from the model
    model.lb(find(model.lb > minLB & model.lb < 0))=-1000;
    models{h11,1}=model
end    
    
modelList1=strcat(modelListN, '_');
modelsM={};
for h2 = 1:size(modelListN,1)
    model=models{h2}
    exchangeCom=intersect(model.rxns, fakemodel1.rxns);
    fakemodelS = removeRxns(fakemodel1,setdiff(fakemodel1.rxns,exchangeCom))
    model = removeRxns(model,Ex_art(:,1));
    model= mergeTwoModels(model, fakemodelS, 1, false);
    if ~isempty(biomass)
        model.lb(find(strcmp(model.rxns,biomass.EXrxn)))=0;
    else
        model.lb(find(strcmp(model.rxns,'Ex_Biomass')))=0;
    end
    model.rxns = strcat(modelList1{h2, 1}, model.rxns);
    Lumen=model.mets(find(cellfun('isempty',strfind(model.mets,'ee[lu]'))))
    if ~isempty(biomass)
        Lumen=union(Lumen,biomass.mets)
    else
        Lumen=union(Lumen,'cpd11416ee[lu]')
    end
    model.mets(find(ismember(model.mets,Lumen)))=strcat(modelList1{h2, 1}, model.mets(find(ismember(model.mets,Lumen))))
    modelsM{h2,1}=model;
end
   
merged=mergeTwoModels(modelsM{1}, modelsM{2}, 1, false)
for i = 3:size(modelsM,1)
    merged= mergeTwoModels(merged, modelsM{i}, 1, false);
end
mergedModelS=mergeTwoModels(merged, fakemodel, 1, false)   
      
% make a global biomass including the biomasses of the bacteria in the community  
if ~isempty(biomass)
     BiomassAll=mergedModelS.mets(find(~(cellfun('isempty',strfind(mergedModelS.mets,biomass.mets)))));
else
     BiomassAll=mergedModelS.mets(find(~(cellfun('isempty',strfind(mergedModelS.mets,'cpd11416ee[lu]')))));
end
biomassmodel.mets=BiomassAll;
biomassmodel.rxns={'BiomassAll'};
biomassmodel.lb=0.1;
biomassmodel.ub=1;
biomassmodel.S=zeros(numel(biomassmodel.mets),numel(biomassmodel.rxns));
if ~isempty(biomass)
    biomassmodel.S(find(strcmp(biomassmodel.mets,biomass.mets)))=1;
else
    biomassmodel.S(find(strcmp(biomassmodel.mets,'cpd11416ee[lu]')))==1;
end
% add bacterial abundance as Stoichiometric Coefficients into the global biomass
for w12 =1:numel(biomassmodel.mets)
    mgs1 = strrep(biomassmodel.mets{w12}, '_cpd11416ee[lu]', '');
    value=abundanceN(find(strcmp(modelListN,mgs1)));
    if ~isempty(value)
      biomassmodel.S(w12,1)=-(value);
    end
end
% remove the FoEx_ and Fo_ global biomass    
if ~isempty(biomass)
    temp3=strrep(biomass.EXrxn, 'Ex_', '');
else
    temp3=strrep('Ex_Biomass', 'Ex_', '');
end
Toremove={strcat('FoEx_',temp3) strcat('Fo_',temp3)}
mergedModelS = removeRxns(mergedModelS,Toremove);
%add the global biomass to the community model
PmergedModel=mergeTwoModels(mergedModelS,biomassmodel, 1, false);
PmergedModel.c(:,1)=0;
PmergedModel.c(find(strcmp(PmergedModel.rxns,'BiomassAll')))=1;
PmergedModel1=PmergedModel;

for jjj=1:numel(modelListN)   
         RXNs=PmergedModel1.rxns(find(~(cellfun('isempty',strfind(PmergedModel1.rxns,modelListN{jjj})))));
         ExRXN=PmergedModel1.rxns(find(~(cellfun('isempty',strfind(PmergedModel1.rxns,'Biomass_Bacteria')))));
         ExRXN=intersect(ExRXN,RXNs);
         PmergedModel1=coupleRxnList2Rxn(PmergedModel1,RXNs,ExRXN); 
         jjj
end
save([PathToSave filesep sampleName{h1} '.mat'],'PmergedModel1','PmergedModel')
report{h1,1}=sampleName{h1}
report{h1,2}=true
end
end

function exchangeMets=GetExchangeMetabolite(model)
indexEx=strfind(model.rxns,'Ex_');
IndexEx = find(not(cellfun('isempty',indexEx)));
S=model.S(:,IndexEx);
exchangeMets=model.mets(find(any(S,2)));
end
