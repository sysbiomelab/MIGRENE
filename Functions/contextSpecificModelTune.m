function [modelTuned,modelInfo]= contextSpecificModelTune(model,MSPInformation,reactionScore,threshold,modelseed)
% this function tune the species-specefic model and gather the model information.
%inputs:
%	model:				species-specefic model
%   MSPInformation:     a structure includes:
%						taxoLevel, the taxonomy names. taxoInfo, taxonomy information for each
% 						species. taxoGroup: taxonomy group fot bacteria. rxns, the reaction
% 						name in reference model. bacteria, a list of MSP IDs. BacteriaNames,
%						list of species names. species, names of the species.
% 						RxnStateAll, the reaction state (absent/present) for each bacteria
%   reactionScore:      a matrix includes 3 different scores for each reaction.
%	threshold:			includes a threshold for each specified taxonomy level.	
%	modelseed			true or false
%                       
% output:
%   modelTuned:     species-specefic model
%	modelInfo: 		provides the gap filling percentage at different levels: i.e. taxonomy 
%					proximity based, further taxonomy level, not annotated (without
%					gene-protein-reaction relation) and the total gap filling percentage

%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk


%find metabolites of exchange reactions
exchangeMets=GetExchangeMetabolite(model)
%remove dead end exchange reactions 
if ~isempty(exchangeMets)
    for hh1=1:size(exchangeMets,1)
        expression = exchangeMets{hh1};
        matchStr = regexp(model.mets,expression,'match');
        Index = find(~(cellfun('isempty',matchStr)));
        len=find(any(model.S(Index,:),1))';
        if length(len) <= 2
            model=removeRxns(model,model.rxns(len));
        end
    end
end

model.description='by MIGRENE toolbox';
if ~isempty(MSPInformation.BacteriaNames)
    index2=find(ismember(MSPInformation.bacteria,MSPInformation.species))
    model.modelName=MSPInformation.BacteriaNames{index2,1};  
end
model.modelID=MSPInformation.species; 
model.compNames={'Extracellular';'Cytosol';'ExtracellularForElectronTransportChain';'boundary'};
model.comps={'C_e';'C_c';'C_pe';'e'};

if modelseed  
    MIGDIR = fileparts(which('MIGRENE_pipeline'));
    METPATH=[MIGDIR filesep 'mat' filesep 'MetInformation.mat'];
    load(METPATH)
    for w1 =1:size(MetInformation,1)
        if ismember(MetInformation{w1,1},model.mets)
            model.metNames(find(strcmp(model.mets, MetInformation{w1,1})),1)=MetInformation(w1,3);
            model.metKEGG(find(strcmp(model.mets, MetInformation{w1,1})),1)=MetInformation(w1,4);
            model.metFormulas(find(strcmp(model.mets, MetInformation{w1,1})),1)=MetInformation(w1,5); 
        end
    end
end
modelTuned=model
% Assign an empty table to gather all the info of the model as below
modelInfo = array2table(zeros(1,8))
modelInfo.Properties.VariableNames = {'number_of_rxns',...
    'number_of_rxns_without_transport_rxns','number_of_rxns_with_genes','level_of_gapfilling','percentage_of_gapfillig'...
    'gapfillig_at_the_level','gapfillig_at_other_level',...
    'gapfillig_at_nonAnnotatedRxns'};
modelInfo.number_of_rxns=length(model.rxns);
indexEx=strfind(model.rxns,'Ex');
IndexEx = find(cellfun('isempty',indexEx));
indexTr=strfind(model.rxns,'t_');
IndexTr = find(cellfun('isempty',indexTr));
indexOfTrExRxns=intersect(IndexEx,IndexTr);
modelInfo.number_of_rxns_without_transport_rxns=length(indexOfTrExRxns);
%get the level of gapfilling
UsedThreshold=threshold(find(threshold,1,'first'));
if isempty(UsedThreshold)
UsedThreshold1=0;
UsedThreshold=0;
else
UsedThreshold1=find(threshold,1,'first');
end
%
if UsedThreshold1 == 0
    level='not classified';
else
    level=MSPInformation.taxoLevel{UsedThreshold1};
end

modelInfo.level_of_gapfilling=level;
rxns=model.rxns(indexOfTrExRxns)
g=find(ismember(MSPInformation.rxns,rxns));
metagenomesetMSP=reactionScore(g,3);
% number of reactions with corresponding genes in the species
modelInfo.number_of_rxns_with_genes=length(find(metagenomesetMSP>=1)); 
% gapfilling info
% percentage of gapfilling
modelInfo.percentage_of_gapfillig=(1-length(find(metagenomesetMSP>=1))/...
    length(indexOfTrExRxns))*100;
% percentage of number of reactions added by gap filling at the lowest classified level
modelInfo.gapfillig_at_the_level=length(find(metagenomesetMSP>=UsedThreshold...
    & metagenomesetMSP <1))/length(indexOfTrExRxns)*100;
% percentage of number of reactions added by gap filling at the other level 
modelInfo.gapfillig_at_other_level=length(find(metagenomesetMSP>=0 ...
    & metagenomesetMSP < UsedThreshold))/length(indexOfTrExRxns)*100; 
% percentage of number of reactions added by gap filling without any info in the catalog 
modelInfo.gapfillig_at_nonAnnotatedRxns=length(find(metagenomesetMSP==-1))...
    /length(indexOfTrExRxns)*100;
end

function exchangeMets=GetExchangeMetabolite(model)
indexEx=strfind(model.rxns,'Ex_');
IndexEx = find(not(cellfun('isempty',indexEx)));
S=model.S(:,IndexEx);
exchangeMets=model.mets(find(any(S,2)));
exchangeMets=strrep(exchangeMets,'ee[e]','');
exchangeMets=strrep(exchangeMets,'[e]','');
exchangeMets=strrep(exchangeMets,'e[e]','');
end