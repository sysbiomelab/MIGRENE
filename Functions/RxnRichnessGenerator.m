function [richness]= RxnRichnessGenerator(modelList,PathToModels,abundance,sampleName)
%inputs:
%   modelList: 			list of model names.
%   PathToModels:		a structure includes the path where the models are available
%						and the name of model assigned in the .mat files
%	abundance:			matrix of microbiome (MSP) abundance profile
%	sampleName:			list of sample names
%outputs:
%   richness: 			gut microbiome reaction composition (reaction richness) of all individuals

%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk

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
    load([PathToModels.path filesep modelList{w1}],PathToModels.name)
	model=eval(PathToModels.name);
    temporary(find(ismember(rxns,model.rxns)),w1)=1;
end
binary1=binary';
rxnsBinary=[];
for i= 1:numel(sampleName)
    abun=binary1(i,:);
    temporary2=temporary;
    for j=1:numel(abun)
        temporary2(:,j)=temporary2(:,j)*abun(:,j);
    end
    rxnsBinary(:,i)=any(temporary2 ==1,2);
end

richness = sum(rxnsBinary,1);
richness=table(sampleName',richness');
richness.Properties.VariableNames = {'sampleName','rxn_richness'};

end
