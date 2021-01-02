function [RxnState, MSPmodel]= MetagemenomeToReactions(model,metagenomeData)
% creates a reaction state for each species based on absent/present genes in MSP into gut
% reference model and filter genes and gene rules in reference model for each species 
%inputs:
%   model: 				reference metabolic Model with COBRA or RAVEN format.
%   metagenomeData:		a structure contains two fields "gene" and "value", includes  gene name and 
%                       value, respectively.
%outputs:
%   MSPmodel: 			reference model with the bacterial genes and gene rules
%   RxnState: 			a vector showing the state of the reaction (zero or one) for the MSP

%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk


if isfield(model,'rules') && ~isfield(model,'grRules')
    model.grRules=cell([numel(model.rxns) 1]);
    for h=1:length(model.rules)
    if  ~isempty(model.rules{h,1})
        matchStr=regexp(model.rules{h,1},'\d*','match');
        collect={};
        for k=1:length(matchStr)
            index=str2num(matchStr{k});
            converted=model.genes{index,1};
            collect=vertcat(collect,converted);
        end
        if length(matchStr) >1
            model.grRules{h,1}=['(' strjoin(unique(collect),' or ') ')']; % generating model.grRules
        else
             model.grRules{h,1}=strjoin(unique(collect),' or ');
        end
    else
        model.grRules(h,1)={''}; 
    end
    end
elseif ~isfield(model,'rules') && ~isfield(model,'grRules')
    error('Either model.rules or model.grRules would be defined in the model. please provide one of them')
end

presentGenes=metagenomeData.gene(find(metagenomeData.value == 1));

RxnState=zeros(length(model.rxns),1);

RxnState(find(cellfun('isempty',model.grRules)),1)=-1;
% tic
% for i=1:length(presentGenes)
%     matchStr = regexp(model.grRules,presentGenes{i},'match');
%     indexx=find(not(cellfun('isempty',matchStr)));
%     RxnState(indexx,1)=1;
% end
% toc
model.grRules=strrep(model.grRules,'( ','');
model.grRules=strrep(model.grRules,'(','');
model.grRules=strrep(model.grRules,') ','');
model.grRules=strrep(model.grRules,')','');

totalgenes={};
for i=1:length(model.grRules)
    temp1=model.grRules{i,1};
    if ~isempty(temp1)
        genes=strsplit(temp1,' or ');
        inter=intersect(genes,presentGenes);
        if ~isempty(inter)
            RxnState(i,1)=1;
        end
        if length(inter) ==1
           model.grRules(i,1)= inter(1,1);
        elseif length(inter)>1
            model.grRules{i,1}=['(' strjoin(unique(inter),' or ') ')'];
        else
            model.grRules{i,1}='';
        end
        totalgenes=vertcat(totalgenes,inter);
    end
end

model.genes=unique(totalgenes);
model.geneNames=model.genes;

temp1=model.grRules;
temp1=strrep(temp1,'( ','');
temp1=strrep(temp1,'(','');
temp1=strrep(temp1,') ','');
temp1=strrep(temp1,')','');
for h=1:length(temp1)
    if  ~isempty(temp1{h,1})
        tra2=strsplit(temp1{h,1},' or ')';
        collect={};
        for k=1:length(tra2)
            index=find(strcmp(model.genes,tra2{k})) ;
            converted=['x(' num2str(index) ')'];
            collect=vertcat(collect,converted);
        end
        if length(tra2) >1
            model.rules{h,1}=['(' strjoin(unique(collect),' | ') ')']; % generating genericModel.rules
        else
            model.rules{h,1}= collect{1,1};
        end
    else
        model.rules(h,1)={''}; 
    end
end


if isfield(model,'rxnGeneMat')
    model=rmfield(model,'rxnGeneMat');
end

MSPmodel=model;
end

