function [countPerFive]= CPFGenerator(modelList,PathToModels,abundance,sampleName)
%inputs:
%   modelList: 				list of model names.
%   PathToModels:			a structure includes the path where the models are available
%							and the name of model assigned in the .mat files
%	abundance:				matrix of microbiome (MSP) abundance profile
%	sampleName:				list of sample names
%outputs:
%   countPerFive:			reactobiome profile

%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk
index=[];
for h1 = 1:size(modelList,1)
	if exist([PathToModels.path filesep modelList{h1} '.mat'])
		index=[index;h1];
	end
end
modelList=modelList(index,:);
abundance=abundance(index,:);
c={};
for h1 = 1:size(modelList,1)
	load([PathToModels.path filesep modelList{h1}],PathToModels.name)
	model=eval(PathToModels.name);
	models{h1,1}=model;
	c=vertcat(c,model.rxns);
end
rxns=unique(c);
compare=zeros(numel(rxns),numel(modelList));
for w1=1:numel(modelList)
	model=models{w1,1};
    compare(find(ismember(rxns,model.rxns)),w1)=1;
end
binary=abundance;
binary(find(binary>0))=1;
binary1=binary';
count_rxnstions=[];
for i= 1:numel(sampleName)
    abun=binary1(i,:);
    compare2=compare;
    for j=1:numel(abun)
        compare2(:,j)=compare2(:,j)*abun(:,j);
    end
    count_rxnstions(:,i)=sum(compare2,2);
end


count_rxnstionsNor=zeros(size(count_rxnstions));
biomassCount=count_rxnstions(find(strcmp(rxns, 'Biomass_Bacteria')),:);
if ~isempty(biomassCount)
	for i=1:size(count_rxnstionsNor,2)
		count_rxnstionsNor(:,i)=count_rxnstions(:,i)*500/biomassCount(1,i);
	end
else
	for i=1:size(count_rxnstionsNor,2)
		count_rxnstionsNor(:,i)=count_rxnstions(:,i)*500/max(count_rxnstions(:,i));
	end
end
count_rxnstionsNor = array2table(count_rxnstionsNor);
count_rxnstionsNor.Properties.VariableNames = sampleName;
countPerFive=[array2table(rxns) count_rxnstionsNor];

end