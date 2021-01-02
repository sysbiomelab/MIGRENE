function [reactionScore, threshold] = MetaGenomicsReactionScore(BacteriaInformation)
% this function calculate reaction score and threshold for gap filling.
%inputs:
%	BacteriaInformation:	a structure includes:
%							taxoLevel, the taxonomy name. taxoInfo, taxonomy information for each
% 							bacteria. taxoGroup: taxonomy group fot each bacteria. rxns, the reaction
% 							name in reference model. bacteria, MSP IDs. BacteriaNames, species name.
% 							RxnStateAll, the reaction state (absent/present) for each bacteria. 
%						    species, names of species				
%                       
% output:
%   reactionScore:     a matrix includes 3 different scores for each reaction.
%	threshold:			includes a threshold for each specified taxonomy level.	

%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk

index=find(strcmp(BacteriaInformation.bacteria, BacteriaInformation.species));
BacteriaInformation.value=BacteriaInformation.RxnStateAll(:,index);

BacteriaInformation.RxnStateAll(find(BacteriaInformation.RxnStateAll <0))=0;
ScoreMatrix=[];
for j=1:size(BacteriaInformation.taxoLevel,2)
	s=BacteriaInformation.taxoGroup(index,j);
    s1=char(BacteriaInformation.taxoInfo(index,j));
    IndexC = isempty(strfind(s1,'unclassified'));
	if IndexC
		group=find(BacteriaInformation.taxoGroup(:,j)== s);
        expression_Group= BacteriaInformation.RxnStateAll(:,group);
        % calculate the freq of each reaction for the taxonomy level
        ScoreMatrix(:,j) = (sum(expression_Group,2))/size(expression_Group,2);
    else
        ScoreMatrix(:,j)=zeros(size(BacteriaInformation.expressionset,1),1);
    end
end

ScoreMatrix(:,j+1)=(sum(ScoreMatrix,2))/size(ScoreMatrix,2);
t1=sort(ScoreMatrix,1,'descend');
for j=1:size(ScoreMatrix,2)-1
    if sum(t1(:,j))>0
        threshold(1,j)=t1(find(t1(:,j)==0,1, 'first')-1,end);
    else
        threshold(1,j)=0;
    end
end

reactionScore(:,1)=ScoreMatrix(:,end);
reactionScore(:,2)=BacteriaInformation.value;
reactionScore(:,end+1)=sum(reactionScore,2);
reactionScore(find(reactionScore(:,2)==-1),end)=-1;
reactionScore(find(reactionScore(:,end)> 1),end)=1;

end
