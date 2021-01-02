function [coverage,enrichment]= pRSEGenerator(modelList,PathToModels,abundance,sampleName,IDmap,GSE,significance)
%inputs:
%   modelList: 				list of model names.
%   PathToModels:			a structure includes the path where the models are available
%							and the name of model assigned in the .mat files
%	abundance:				matrix of microbiome (MSP) abundance profile
%	sampleName:				list of sample names
%	IDmap:					ID mapping between reaction names in the models and and 
% 							the IDs in the pathway file
%	GSE:					the pathway profile incuding reaction sets
%OPTIONAL INPUT
%	significance:			p value
%outputs:					
%   coverage:				a table shows the coverage of each pathway in samples
%   enrichment:				a table shows the p value of the pathways in samples
%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk
if nargin<7
    significance=1;
end

index=[];
for h1 = 1:size(modelList,1)
	if exist([PathToModels.path filesep modelList{h1} '.mat'])
		index=[index;h1];
	end
end
modelList=modelList(index,:);
abundance=abundance(index,:);
binary=abundance;
binary(find(binary>0))=1;
rxnTemp={};
for h1 = 1:size(modelList,1)
	load([PathToModels.path filesep modelList{h1}],PathToModels.name)
	model=eval(PathToModels.name);
	rxnTemp=vertcat(rxnTemp,model.rxns);
end
rxns=unique(rxnTemp);
temporary=zeros(numel(rxns),numel(modelList));
for w1=1:numel(modelList)
    load([PathToModels.path filesep modelList{h1}],PathToModels.name)
	model=eval(PathToModels.name);
    temporary(find(ismember(rxns,model.rxns)),w1)=1;
end
binary1=binary';
input=[];
for i= 1:numel(sampleName)
    abun=binary1(i,:);
    temporary2=temporary;
    for j=1:numel(abun)
        temporary2(:,j)=temporary2(:,j)*abun(:,j);
    end
    input(:,i)=any(temporary2 ==1,2);
end
% convert rxns ID to provided ID for enrichment analysis
rxnsTemp=rxns;
for ii=1:size(rxnsTemp,1)
    index=find(strcmp(IDmap(:,1), rxnsTemp{ii,1}));
    if ~isempty(index)
    rxnsTemp{ii,2}=IDmap{index,2};
    end
end
Index = find(not(cellfun('isempty',rxnsTemp(:,2))));
input=input(Index,:);
rxns=rxnsTemp(Index,2);
% calculate Total number of unique IDs in provided GSE.
IDs={};
for i=1:size(GSE,1)
    temp1=strsplit(GSE{i,2},',');
    IDs=vertcat(IDs,temp1');
end
IDs=unique(IDs);
N=size(IDs,1);
enrichment=zeros(size(GSE,1),size(input,2));
coverage=zeros(size(GSE,1),size(input,2));
for i=1:size(GSE,1)
    temp2=transpose(strsplit(GSE{i,2},','));
    m=size(temp2,1); %Number of IDs associated to the term 
    n = N - m; %Number of IDs not associated to the term
    for j=1:size(input,2)
        input1=input(:,j);
        tempOne=find(input1==1);
        g = length(tempOne); % Number of submitted reactions
        k = length(intersect(IDs,rxns(tempOne))); % number of submitted 
        x = length(intersect(temp2,rxns(tempOne))); % reactions with at least one annotation in IDmap
        %number of IDs in the term present in the sample
        enrichment(i,j)=hygecdf(x-1,N,m,k,'upper');
        coverage(i,j)=x/m;
    end
end
if significance ~=1
    coverage(find(enrichment>significance))=0;
end
% convert coverage matrix to table and add sample name and the terms to the table
coverage = array2table(coverage);
coverage.Properties.VariableNames = sampleName;
coverage=[array2table(GSE(:,1)) coverage];
% convert enrichment matrix to table and add sample name and the terms to the table
enrichment = array2table(enrichment);
enrichment.Properties.VariableNames = sampleName;
enrichment=[array2table(GSE(:,1)) enrichment];

end