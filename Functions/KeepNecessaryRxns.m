function reducedModelTemp = KeepNecessaryRxns(model, score, threshold, min)
%inputs:
%   model: 				metabolic model with COBRA or RAVEN format.
%   score:				a numeric vector that shows the score of each reaction in the model
%	min:				minimum fraction of objective
%outputs:
%   reducedModelTemp: 			reference model with the bacterial genes and gene rules

%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk
[tempModel,~,IndexRev2irrev,IndexIrrev2rev] = convertToIrreversible(model);
expressionRxnsIrrev = zeros(length(tempModel.rxns),1);
for i1=1:length(tempModel.rxns)
    expressionRxnsIrrev(i1,1) = score(IndexIrrev2rev(i1,1),1);
end
expressionRxnsIrrev(find(expressionRxnsIrrev==-1))=0;
cc=optimizeCbModel(model);
tempModel.lb(find(tempModel.c ==1),1)=cc.f*min; % minimum fraction of objective 
tempModel.c(:,1)=0;
for i1=1:length(tempModel.rxns)
    if (expressionRxnsIrrev(i1,1) < threshold)
        tempModel.c(i1,1) = threshold-expressionRxnsIrrev(i1,1); %FIX: use expression level as weight
    end
end
gimmeSolution = optimizeCbModel(tempModel,'min');
reactionScoreTransition=zeros(length(expressionRxnsIrrev),1);
if (gimmeSolution.stat ~= 1)
    reactionScoreTransition(:,1) = 0;
end
reactionScoreTransition(find(gimmeSolution.x>0),1)=1;
reactionScoreTransition(find(expressionRxnsIrrev>threshold))=1;
%Translate reactionActivity to reversible model
reactionScoreRev = zeros(length(model.rxns),1);
for i=1:length(model.rxns)
    temp1=IndexRev2irrev{i,1}';
    for j=1:length(temp1)
        if reactionScoreTransition(temp1(j)) > 0
            reactionScoreRev(i,1) = reactionScoreTransition(temp1(j));
        end
    end
end
rxn2remove = model.rxns(reactionScoreRev == 0);
reducedModelTemp = removeRxns(model,rxn2remove); 
end
