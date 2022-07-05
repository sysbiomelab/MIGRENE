function [reactionRelativeAbun, rxnAbunPerSample]= ReactionAbundanceGenerator(modelList,PathToModels,abundance,sampleName)
%inputs:
%   modelList: 				list of model names.
%   PathToModels:			a structure includes the path where the models are available
%							and the name of model assigned in the .mat files
%	abundance:				matrix of microbiome (MSP) abundance profile
%	sampleName:				list of sample names
%outputs:
%   reactionRelativeAbun	relative reaction abundance
%	rxnAbunPerSample		reaction abundance
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
abundance1=abundance';
FinalAbundance=[];
for i= 1:numel(sampleName)
    abun1=abundance1(i,:);
    compare1=compare;
     for j=1:numel(abun1)
        compare1(:,j)=compare1(:,j)*abun1(:,j);
     end
    FinalAbundance(:,i)=sum(compare1,2);
end

SumAbundancy=sum(FinalAbundance);
reactionRelativeAbun=[];
for j=1:numel(SumAbundancy)
   reactionRelativeAbun(:,j)=FinalAbundance(:,j)/SumAbundancy(:,j);
end
reactionRelativeAbun = array2table(reactionRelativeAbun);
reactionRelativeAbun.Properties.VariableNames = sampleName;
reactionRelativeAbun=[array2table(rxns) reactionRelativeAbun];
 
rxnAbunPerSample=zeros(size(FinalAbundance));
for i=1:size(FinalAbundance,2)
	x = FinalAbundance(:,i);
	minVal = min(x);
	maxVal = max(x);
	if minVal==maxVal
		rxnAbunPerSample(:,i)=0;
	else
		rxnAbunPerSample(:,i) = (x - minVal) / ( maxVal - minVal);
	end
end
rxnAbunPerSample = array2table(rxnAbunPerSample);
rxnAbunPerSample.Properties.VariableNames = sampleName;
rxnAbunPerSample=[array2table(rxns) rxnAbunPerSample];
end