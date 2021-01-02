function [MSPInformation]= GenerateMSPInformation(taxo,RXNDIR,model)
% this function gather the taxonomy information and reaction state for bacteria.
%inputs:
%	taxo:				taxonomy profile
%   RXNDIR:             path to where RxnState and modelforMSP for each bacterium (MSP) was saved 
%   model:        		reference metabolic model with COBRA or RAVEN format, 

%                       
% output:
%   MSPInformation:     Structure includes:
%						taxoLevel, the taxonomy name. taxoInfo, taxonomy information for each
% 						bacteria. taxoGroup: taxonomy group fot each bacteria. rxns, the reaction
% 						name in reference model. bacteria, MSP IDs. BacteriaNames, species name.
% 						RxnStateAll, the reaction state (absent/present) for each bacteria

%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk

% read taxo info
T = readtable(taxo);
% get the bacterial species name
mspNames=table2cell(T(:,1));
% provide the taxo info from bottom to top taxonomical levels  (genus to
% phylum level)
TaxoAll=T.Properties.VariableNames;
levels={'genus' 'family' 'order' 'class' 'phylum'};
Index=[];
for i=1:length(levels)
    index=find(ismember(TaxoAll,levels{i}));
    if ~isempty(index)
        Index=[Index index];
    end
end
Taxo=TaxoAll(Index); % name of sorted levels
 if isempty(Index)
     error('there is no taxonomy info in the provided excel file. make sure the first row provide the taxonomy name i.e. genus to phylum')
 else
    infoFile=table2cell(T(1:end ,Index)); % taxonomy info for each species
    %dedicate the groups in each taxonomy level 
    infoFile1=zeros(size(infoFile,1),size(infoFile,2));
    for i=1:size(Taxo,2)
	[~,~,ic]=unique(infoFile(:,i));
	infoFile1(:,i)=ic;
    end
 end

% collect all the info in a structure
MSPInformation.taxoLevel=Taxo;
MSPInformation.taxoInfo=infoFile;
MSPInformation.taxoGroup=infoFile1;
MSPInformation.rxns=model.rxns;
MSPInformation.bacteria=mspNames;

index1=find(ismember(TaxoAll,'species'));
if ~isempty(index1)
    MSPInformation.BacteriaNames=table2cell(T(1:end ,index1));
else
    MSPInformation.BacteriaNames={};
end
    
% all the RxnStates were collected from reactionProfile directory
MSPInformation.RxnStateAll=[];
for i =1:numel(MSPInformation.bacteria)
    if exist([RXNDIR filesep MSPInformation.bacteria{i} '.mat'])
        load ([RXNDIR filesep MSPInformation.bacteria{i} '.mat'],'RxnState');
        if ~isempty(ismember(model.rxns, 'Biomass_Bacteria'))
            RxnState(find(strcmp(model.rxns, 'Biomass_Bacteria')))=1;
        else 
            matchStr = regexp(lower(model.rxns),'biomass','match');
            RxnState(find(not(cellfun('isempty',matchStr))))=1;
        end
       MSPInformation.RxnStateAll=[MSPInformation.RxnStateAll RxnState];
    else
        error(['there is no information for ' MSPInformation.bacteria{i} '. please check the dedicated directory. Besides, you might not generate it by MetagemenomeToReactions function'])
    end
end
