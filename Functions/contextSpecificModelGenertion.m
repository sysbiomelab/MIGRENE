function [contextSpecificModel] = contextSpecificModelGenertion(model,metagenomeset,threshold,bibliome)
% generates species-specific model based on reaction score, threshold and mind the gap.
%inputs:
%   model:				reference metabolic model in COBRA or RAVEN format.
%	metagenomeset:		a matrix includes 3 different scores for each reaction.
%	threshold:			includes a threshold for each specified taxonomy level.
%OPTIONAL INPUTS:
%	bibliome:			any bibliome data on phenotypic features of the species.
% 						as structure with four fields:
% 						"bacteria" a cell array listing the name of the bacteria; "rxn" list the
% 						name of the reactions having bibliome; "value" a matrix of numbers: zero
% 						means no information, 1 means consumed, 2 means produced, -1 not-consumed
% 						and -2 means not-produced by the corresponding bacteria. "aerobIenfo" a
% 						cell array provides the info that the bacteria require oxygen for growth
% 						or not so specefiy with "aerobe" or "anaerobe" or "facultative".
%outputs: 
%	contextSpecificModel species-specific model	
%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk
if nargin<4
    bibliome={};
end
% Integrating the bibliome data into the model
if ~isempty(bibliome)
    if isfield(bibliome, 'value') && isfield(bibliome, 'rxn')
        %check all necessaries are provided
        if size(bibliome.value,2)>1 && ~isfield(bibliome, 'species')
            error('there are more than one column in bibliome.value and you didnt specify the name of species. please assign it by adding species field or do not dedicate bibliome as input of function')
        elseif size(bibliome.value,2)>1 && isfield(bibliome, 'species')
            % find the corresponding bibliome data for the species
            index=find(strcmp(bibliome.bacteria, bibliome.species));
            value=bibliome.value(:,index);
        elseif size(bibliome.value,2)==1
            value=bibliome.value;
        end
         model1=model;
        % find the exchange reactions and constrain the model for consumption 
        consumed=bibliome.rxn(find(value==1));
        model1.lb(find(ismember(model1.rxns,consumed)))=-1;
        model1.ub(find(ismember(model1.rxns,consumed)))=0;
        % find the exchange reactions and constrain the model for
        % production
        produced=bibliome.rxn(find(value==2));
        model1.lb(find(ismember(model1.rxns,produced)))=0.1;
        model1.ub(find(ismember(model1.rxns,produced)))=1000;
        % find the exchange reactions and constrain the model for
        % metabolites regarding the bibliome data
        notCon=bibliome.rxn(find(value==-1));
        model1.lb(find(ismember(model1.rxns,notCon)))=0;
        notPro=bibliome.rxn(find(value==-1));
        model1.ub(find(ismember(model1.rxns,notPro)))=0;
       % check the changes doesnt affect the functionality of the model
        g2=optimizeCbModel(model1);
        if g2.f>0
            disp('the bibliome info was added to the model')
            model=model1;
        else
            disp('the bibliome info could not be added to the model. some changes made the model nonfunctional')
        end
    end
end
%check the species requires oxygen for growth
if ~isempty(bibliome) && isfield(bibliome, 'aerobeInfo')
    if size(bibliome.aerobeInfo,1)>1 && ~isfield(bibliome, 'species')
         error('there are more than one row in bibliome.aerobeInfo and you didnt specify the name of species. please assign it by adding species field or do not dedicate bibliome as input of function.')
    elseif size(bibliome.aerobeInfo,1)>1 && isfield(bibliome, 'species')
        index=find(strcmp(bibliome.bacteria, bibliome.species));
        obic=bibliome.aerobeInfo(index,1);
    else
        obic=bibliome.aerobeInfo;
    end
else
    obic={'none'};
end

%generate the specefic model
if strcmp(obic,'anaerobe')
model.ub(find(strcmp(model.rxns, 'Ex_O2')))=0;
model.lb(find(strcmp(model.rxns, 'Ex_O2')))=0;
model.ub(find(strcmp(model.rxns, 'rxn08173')))=1000;
model.lb(find(strcmp(model.rxns, 'rxn08173')))=0;
end

if ~ismember(obic,{'aerobe';'anaerobe';'facultative'})
model.ub(find(strcmp(model.rxns, 'Ex_O2')))=0;
model.lb(find(strcmp(model.rxns, 'Ex_O2')))=-10;
model.ub(find(strcmp(model.rxns, 'rxn08173')))=1000;
model.lb(find(strcmp(model.rxns, 'rxn08173')))=0;
end

if strcmp(obic,'aerobe')
model.ub(find(strcmp(model.rxns, 'Ex_O2')))=-1;
model.lb(find(strcmp(model.rxns, 'Ex_O2')))=-19.2;
model.ub(find(strcmp(model.rxns, 'rxn08173')))=1000;
model.lb(find(strcmp(model.rxns, 'rxn08173')))=0.5;
if metagenomeset(find(strcmp(reference_model.rxns,'rxn08173')),3) < 1
	metagenomeset(find(strcmp(reference_model.rxns,'rxn08173')),3)=1
end
end
%get the closest taxonomy level
s=threshold(find(threshold,1,'first'));
if isempty(s)
threshold=0;
else
threshold=s;
end

if ~strcmp(obic,'facultative')
    % collect the reactions with score from other levels of taxonomy and check
    % the essentiality of the reactions.
    rxn2remove= model.rxns(find(metagenomeset(:,3)<threshold & metagenomeset(:,3)>=0),1);
    tissueModel= removeRxns(model,rxn2remove);
    cc=optimizeCbModel(tissueModel);
    matrix_General(1,1)=cc.f;
    exp=metagenomeset(:,3);
    TempModel = KeepNecessaryRxns(model, exp, 0.99, 0.1);
    ddd=setdiff(TempModel.rxns,tissueModel.rxns);
    ddd_general=ddd;
    matrix_General(1,2)=length(setdiff(TempModel.rxns,tissueModel.rxns));
    metagenomeset1=metagenomeset;
    for t=1:length(ddd)
        metagenomeset1(find(strcmp(model.rxns, ddd{t})),3)=1;
    end
    % remove the reactions with score from other levels of taxonomy
    rxn2remove=model.rxns(find(metagenomeset1(:,3)<threshold & metagenomeset1(:,3)>=0),1) ;
    tissueModel= removeRxns(model,rxn2remove);

    exp3=metagenomeset(:,3);
    tissueModel1=tissueModel;
    %remove deadEnd reactions that are not supported by metagenomics. 
    [~,~, removedRxns] = removeDeadEnds(tissueModel1);
    % exp3(:,2)=1;
    % exp3(find(ismember(tissueModel1.rxns,removedRxns)),2)=0;
    % a1=tissueModel1.rxns(find(exp3(:,1) <= 1));
    % b1=tissueModel1.rxns(find(exp3(:,2) == 0));
    % g1=intersect(a1,b1);
    % tissueModel1= removeRxns(tissueModel,g1);
    tissueModel1= removeRxns(tissueModel1,removedRxns);
    fba=optimizeCbModel(tissueModel1);

    % prunnig the reactions without GPR

    exp4=zeros(length(tissueModel1.rxns),1);
    for q=1:length(tissueModel1.rxns)
        if	ismember(tissueModel1.rxns(q),model.rxns)
            exp4(q,1)=metagenomeset(find(strcmp(model.rxns,tissueModel1.rxns(q))),3);
        end
    end
    tissueModelF5 = KeepNecessaryRxns(tissueModel1, exp4, 0.0001, 0.1);

    AA2=setdiff(tissueModel1.rxns,tissueModelF5.rxns);

    % get all the reaction with score 1 excluding the transport and exchange
    % reactions
    FromMetaGenomics=model.rxns(find(exp(:,1) >= 1));
    indexEx=strfind(model.rxns,'Ex');
    IndexEx = find(not(cellfun('isempty',indexEx)));
    indexTr=strfind(model.rxns,'t_');
    IndexTr = find(not(cellfun('isempty',indexTr)));
    indexOfTrEx=union(IndexEx,IndexTr);
    TrEx=model.rxns(indexOfTrEx);
    FromMetaGenomics=setdiff(FromMetaGenomics,TrEx);
    % get the dead end reactions with score 1 
    tissueModelF5=tissueModel1;
    index=1:50:length(AA2);
    for h1=1:length(index)
        if h1 ~= length(index)
            [~, ~, removedRxns] = removeDeadEnds(tissueModelF5);
            d=AA2(index(h1):(index(h1)+49));
            T1=tissueModelF5;
            T1= removeRxns(T1,d);
            [~, ~, removedRxns1] = removeDeadEnds(T1);
            removed=setdiff(removedRxns1,removedRxns);
            hhh=intersect(removed,FromMetaGenomics);
            if  length(hhh)<10
                tissueModelF5=T1;
            else
                for h2=1:length(d)
                    T1=tissueModelF5;
                    T1= removeRxns(T1,d{h2});
                    [~, ~, removedRxns1] = removeDeadEnds(T1);
                    removed=setdiff(removedRxns1,removedRxns);
                    hhh=intersect(removed,FromMetaGenomics);
                    if  length(hhh)<10
                        tissueModelF5=T1;
                    end
                end
            end
        else
            [~, ~, removedRxns] = removeDeadEnds(tissueModelF5);
            d=AA2(index(h1):end);
            for h3=1:length(d)
                    T1=tissueModelF5;
                    T1= removeRxns(T1,d{h3});
                    [~, ~, removedRxns1] = removeDeadEnds(T1);
                    removed=setdiff(removedRxns1,removedRxns);
                    hhh=intersect(removed,FromMetaGenomics);
                    if  length(hhh)<10
                        tissueModelF5=T1;
                    end
            end
        end
    end

    exp4=zeros(length(tissueModelF5.rxns),1);
    for q=1:length(tissueModelF5.rxns)
        if	ismember(tissueModelF5.rxns(q),model.rxns)
            exp4(q,1)=metagenomeset(find(strcmp(model.rxns,tissueModelF5.rxns(q))),3);
        end
    end

    AA1=tissueModelF5.rxns(find(exp4<1 & exp4> -1));
    BB1=[];
    for q=1:length(AA1)
        if	ismember(AA1(q),model.rxns)
            BB1(q,1)=metagenomeset(find(strcmp(model.rxns,AA1{q})),3);
        else 
            BB1(q,1)=1;
        end
    end


    T=table(AA1,BB1);
    Sort_Table = sortrows(T,'BB1');
    AA2=table2cell(Sort_Table(:,1));
    %prune the reactions if they are not lethal for the network and functionality of the model
    fba=optimizeCbModel(tissueModelF5);
    index=1:50:length(AA2);
    for h1=1:length(index)
        if h1 ~= length(index)
            [~, ~, removedRxns] = removeDeadEnds(tissueModelF5);
            d=AA2(index(h1):(index(h1)+49));
            T1=tissueModelF5;
            T1= removeRxns(T1,d);
            [~, ~, removedRxns1] = removeDeadEnds(T1);
            removed=setdiff(removedRxns1,removedRxns);
            hhh=intersect(removed,FromMetaGenomics);
            fbaa=optimizeCbModel(T1);
            if  length(hhh)<10  && fbaa.f > fba.f*0.1
                tissueModelF5=T1;
            else
                for h2=1:length(d)
                    T1=tissueModelF5;
                    T1= removeRxns(T1,d{h2});
                    [~, ~, removedRxns1] = removeDeadEnds(T1);
                    removed=setdiff(removedRxns1,removedRxns);
                    hhh=intersect(removed,FromMetaGenomics);
                    fbaa=optimizeCbModel(T1);
                    if  length(hhh)<10  && fbaa.f > fba.f*0.1
                        tissueModelF5=T1;
                    end
                end
            end
        else
            [~, ~, removedRxns] = removeDeadEnds(tissueModelF5);
            d=AA2(index(h1):end);
            for h3=1:length(d)
                    T1=tissueModelF5;
                    T1= removeRxns(T1,d{h3});
                    [~, ~, removedRxns1] = removeDeadEnds(T1);
                    removed=setdiff(removedRxns1,removedRxns);
                    hhh=intersect(removed,FromMetaGenomics);
                    fbaa=optimizeCbModel(T1);
                    if  length(hhh)<10 && fbaa.f > fba.f*0.1
                        tissueModelF5=T1;
                    end
            end
        end
    end

    contextSpecificModel=tissueModelF5;
    %matrix_General(1,3)=length(contextSpecificModel.rxns);
    %matrix_General(1,4)=length(find(exp4>=1));
    %matrix_General(1,5)=length(find(exp4<=0));
    %matrix_General(1,6)=matrix_General(1,3)-(matrix_General(1,5)+matrix_General(1,4));
    contextSpecificModel.lb(find(strcmp(contextSpecificModel.rxns, 'rxn08173')))=0;
    contextSpecificModel.ub(find(strcmp(contextSpecificModel.rxns, 'rxn08173')))=1000;
    if exist('produced','var')
        contextSpecificModel.lb(find(ismember(contextSpecificModel.rxns,produced)))=0;
        contextSpecificModel.ub(find(ismember(contextSpecificModel.rxns,produced)))=1000;
    end
end

if strcmp(obic,'facultative')
	model.ub(find(strcmp(model.rxns, 'Ex_O2')))=0;
	model.lb(find(strcmp(model.rxns, 'Ex_O2')))=0;
	model.ub(find(strcmp(model.rxns, 'rxn08173')))=1000;
	model.lb(find(strcmp(model.rxns, 'rxn08173')))=0;
    % collect the reactions with score from other levels of taxonomy and check
    % the essentiality of the reactions.
    rxn2remove= model.rxns(find(metagenomeset(:,3)<threshold & metagenomeset(:,3)>=0),1);
    tissueModel= removeRxns(model,rxn2remove);
    cc=optimizeCbModel(tissueModel);
    matrix_General(1,1)=cc.f;
    exp=metagenomeset(:,3);
    TempModel = KeepNecessaryRxns(model, exp, 0.99, 0.1);
    ddd=setdiff(TempModel.rxns,tissueModel.rxns);
    ddd_general=ddd;
    matrix_General(1,2)=length(setdiff(TempModel.rxns,tissueModel.rxns));
    metagenomeset1=metagenomeset;
    for t=1:length(ddd)
        metagenomeset1(find(strcmp(model.rxns, ddd{t})),3)=1;
    end
    % remove the reactions with score from other levels of taxonomy
    rxn2remove=model.rxns(find(metagenomeset1(:,3)<threshold & metagenomeset1(:,3)>=0),1) ;
    tissueModel= removeRxns(model,rxn2remove);

    exp3=metagenomeset(:,3);
    tissueModel1=tissueModel;
    %remove deadEnd reactions that are not supported by metagenomics. 
    [~,~, removedRxns] = removeDeadEnds(tissueModel1);
    % exp3(:,2)=1;
    % exp3(find(ismember(tissueModel1.rxns,removedRxns)),2)=0;
    % a1=tissueModel1.rxns(find(exp3(:,1) <= 1));
    % b1=tissueModel1.rxns(find(exp3(:,2) == 0));
    % g1=intersect(a1,b1);
    % tissueModel1= removeRxns(tissueModel,g1);
    tissueModel1= removeRxns(tissueModel1,removedRxns);
    fba=optimizeCbModel(tissueModel1);

    % prunnig the reactions without GPR

    exp4=zeros(length(tissueModel1.rxns),1);
    for q=1:length(tissueModel1.rxns)
        if	ismember(tissueModel1.rxns(q),model.rxns)
            exp4(q,1)=metagenomeset(find(strcmp(model.rxns,tissueModel1.rxns(q))),3);
        end
    end
    tissueModelF5 = KeepNecessaryRxns(tissueModel1, exp4, 0.0001, 0.1);

    AA2=setdiff(tissueModel1.rxns,tissueModelF5.rxns);

    % get all the reaction with score 1 excluding the transport and exchange
    % reactions
    FromMetaGenomics=model.rxns(find(exp(:,1) >= 1));
    indexEx=strfind(model.rxns,'Ex');
    IndexEx = find(not(cellfun('isempty',indexEx)));
    indexTr=strfind(model.rxns,'t_');
    IndexTr = find(not(cellfun('isempty',indexTr)));
    indexOfTrEx=union(IndexEx,IndexTr);
    TrEx=model.rxns(indexOfTrEx)
    FromMetaGenomics=setdiff(FromMetaGenomics,TrEx);
    % get the dead end reactions with score 1 
    tissueModelF5=tissueModel1;
    index=1:50:length(AA2);
    for h1=1:length(index)
        if h1 ~= length(index)
            [~, ~, removedRxns] = removeDeadEnds(tissueModelF5);
            d=AA2(index(h1):(index(h1)+49));
            T1=tissueModelF5;
            T1= removeRxns(T1,d);
            [~, ~, removedRxns1] = removeDeadEnds(T1);
            removed=setdiff(removedRxns1,removedRxns);
            hhh=intersect(removed,FromMetaGenomics);
            if  length(hhh)<10
                tissueModelF5=T1;
            else
                for h2=1:length(d)
                    T1=tissueModelF5;
                    T1= removeRxns(T1,d{h2});
                    [~, ~, removedRxns1] = removeDeadEnds(T1);
                    removed=setdiff(removedRxns1,removedRxns);
                    hhh=intersect(removed,FromMetaGenomics);
                    if  length(hhh)<10
                        tissueModelF5=T1;
                    end
                end
            end
        else
            [~, ~, removedRxns] = removeDeadEnds(tissueModelF5);
            d=AA2(index(h1):end);
            for h3=1:length(d)
                    T1=tissueModelF5;
                    T1= removeRxns(T1,d{h3});
                    [~, ~, removedRxns1] = removeDeadEnds(T1);
                    removed=setdiff(removedRxns1,removedRxns);
                    hhh=intersect(removed,FromMetaGenomics);
                    if  length(hhh)<10
                        tissueModelF5=T1;
                    end
            end
        end
    end

    exp4=zeros(length(tissueModelF5.rxns),1);
    for q=1:length(tissueModelF5.rxns)
        if	ismember(tissueModelF5.rxns(q),model.rxns)
            exp4(q,1)=metagenomeset(find(strcmp(model.rxns,tissueModelF5.rxns(q))),3);
        end
    end

    AA1=tissueModelF5.rxns(find(exp4<1 & exp4> -1));
    BB1=[];
    for q=1:length(AA1)
        if	ismember(AA1(q),model.rxns)
            BB1(q,1)=metagenomeset(find(strcmp(model.rxns,AA1{q})),3);
        else 
            BB1(q,1)=1;
        end
    end


    T=table(AA1,BB1);
    Sort_Table = sortrows(T,'BB1');
    AA2=table2cell(Sort_Table(:,1));
    %prune the reactions if they dont collapse the network and functionality of the model
    fba=optimizeCbModel(tissueModelF5);
    index=1:50:length(AA2);
    for h1=1:length(index)
        if h1 ~= length(index)
            [~, ~, removedRxns] = removeDeadEnds(tissueModelF5);
            d=AA2(index(h1):(index(h1)+49));
            T1=tissueModelF5;
            T1= removeRxns(T1,d);
            [~, ~, removedRxns1] = removeDeadEnds(T1);
            removed=setdiff(removedRxns1,removedRxns);
            hhh=intersect(removed,FromMetaGenomics);
            fbaa=optimizeCbModel(T1);
            if  length(hhh)<10  && fbaa.f > fba.f*0.1
                tissueModelF5=T1;
            else
                for h2=1:length(d)
                    T1=tissueModelF5;
                    T1= removeRxns(T1,d{h2});
                    [~, ~, removedRxns1] = removeDeadEnds(T1);
                    removed=setdiff(removedRxns1,removedRxns);
                    hhh=intersect(removed,FromMetaGenomics);
                    fbaa=optimizeCbModel(T1);
                    if  length(hhh)<10  && fbaa.f > fba.f*0.1
                        tissueModelF5=T1;
                    end
                end
            end
        else
            [~, ~, removedRxns] = removeDeadEnds(tissueModelF5);
            d=AA2(index(h1):end);
            for h3=1:length(d)
                    T1=tissueModelF5;
                    T1= removeRxns(T1,d{h3});
                    [~, ~, removedRxns1] = removeDeadEnds(T1);
                    removed=setdiff(removedRxns1,removedRxns);
                    hhh=intersect(removed,FromMetaGenomics);
                    fbaa=optimizeCbModel(T1);
                    if  length(hhh)<10 && fbaa.f > fba.f*0.1
                        tissueModelF5=T1;
                    end
            end
        end
    end

    model_anaerobic=tissueModelF5;
    %matrix_General(1,3)=length(model_anaerobic.rxns);
    %matrix_General(1,4)=length(find(exp4>=1));
    %matrix_General(1,5)=length(find(exp4<=0));
    %matrix_General(1,6)=matrix_General(1,3)-(matrix_General(1,5)+matrix_General(1,4));
	model.ub(find(strcmp(model.rxns, 'Ex_O2')))=-1;
	model.lb(find(strcmp(model.rxns, 'Ex_O2')))=-19.2;
	model.ub(find(strcmp(model.rxns, 'rxn08173')))=1000;
	model.lb(find(strcmp(model.rxns, 'rxn08173')))=0.5;
	if metagenomeset(find(strcmp(reference_model.rxns,'rxn08173')),3) < 1
	metagenomeset(find(strcmp(reference_model.rxns,'rxn08173')),3)=1;
	end 
    % collect the reactions with score from other levels of taxonomy and check
    % the essentiality of the reactions.
    rxn2remove= model.rxns(find(metagenomeset(:,3)<threshold & metagenomeset(:,3)>=0),1);
    tissueModel= removeRxns(model,rxn2remove);
    cc=optimizeCbModel(tissueModel);
    matrix_General(1,1)=cc.f;
    exp=metagenomeset(:,3);
    TempModel = KeepNecessaryRxns(model, exp, 0.99, 0.1);
    ddd=setdiff(TempModel.rxns,tissueModel.rxns);
    ddd_general=ddd;
    matrix_General(1,2)=length(setdiff(TempModel.rxns,tissueModel.rxns));
    metagenomeset1=metagenomeset;
    for t=1:length(ddd)
        metagenomeset1(find(strcmp(model.rxns, ddd{t})),3)=1;
    end
    % remove the reactions with score from other levels of taxonomy
    rxn2remove=model.rxns(find(metagenomeset1(:,3)<threshold & metagenomeset1(:,3)>=0),1) ;
    tissueModel= removeRxns(model,rxn2remove);

    exp3=metagenomeset(:,3);
    tissueModel1=tissueModel;
    %remove deadEnd reactions that are not supported by metagenomics. 
    [~,~, removedRxns] = removeDeadEnds(tissueModel1);
    % exp3(:,2)=1;
    % exp3(find(ismember(tissueModel1.rxns,removedRxns)),2)=0;
    % a1=tissueModel1.rxns(find(exp3(:,1) <= 1));
    % b1=tissueModel1.rxns(find(exp3(:,2) == 0));
    % g1=intersect(a1,b1);
    % tissueModel1= removeRxns(tissueModel,g1);
    tissueModel1= removeRxns(tissueModel1,removedRxns);
    fba=optimizeCbModel(tissueModel1);

    % prunnig the reactions without GPR

    exp4=zeros(length(tissueModel1.rxns),1);
    for q=1:length(tissueModel1.rxns)
        if	ismember(tissueModel1.rxns(q),model.rxns)
            exp4(q,1)=metagenomeset(find(strcmp(model.rxns,tissueModel1.rxns(q))),3);
        end
    end
    tissueModelF5 = KeepNecessaryRxns(tissueModel1, exp4, 0.0001, 0.1);

    AA2=setdiff(tissueModel1.rxns,tissueModelF5.rxns);

    % get all the reaction with score 1 excluding the transport and exchange
    % reactions
    FromMetaGenomics=model.rxns(find(exp(:,1) >= 1));
    indexEx=strfind(model.rxns,'Ex');
    IndexEx = find(not(cellfun('isempty',indexEx)));
    indexTr=strfind(model.rxns,'t_');
    IndexTr = find(not(cellfun('isempty',indexTr)));
    indexOfTrEx=union(IndexEx,IndexTr);
    TrEx=model.rxns(indexOfTrEx)
    FromMetaGenomics=setdiff(FromMetaGenomics,TrEx);
    % get the dead end reactions with score 1 
    tissueModelF5=tissueModel1
    index=1:50:length(AA2);
    for h1=1:length(index)
        if h1 ~= length(index)
            [~, ~, removedRxns] = removeDeadEnds(tissueModelF5);
            d=AA2(index(h1):(index(h1)+49));
            T1=tissueModelF5;
            T1= removeRxns(T1,d);
            [~, ~, removedRxns1] = removeDeadEnds(T1);
            removed=setdiff(removedRxns1,removedRxns);
            hhh=intersect(removed,FromMetaGenomics);
            if  length(hhh)<10
                tissueModelF5=T1;
            else
                for h2=1:length(d)
                    T1=tissueModelF5;
                    T1= removeRxns(T1,d{h2});
                    [~, ~, removedRxns1] = removeDeadEnds(T1);
                    removed=setdiff(removedRxns1,removedRxns);
                    hhh=intersect(removed,FromMetaGenomics);
                    if  length(hhh)<10
                        tissueModelF5=T1;
                    end
                end
            end
        else
            [~, ~, removedRxns] = removeDeadEnds(tissueModelF5);
            d=AA2(index(h1):end);
            for h3=1:length(d)
                    T1=tissueModelF5;
                    T1= removeRxns(T1,d{h3});
                    [~, ~, removedRxns1] = removeDeadEnds(T1);
                    removed=setdiff(removedRxns1,removedRxns);
                    hhh=intersect(removed,FromMetaGenomics);
                    if  length(hhh)<10
                        tissueModelF5=T1;
                    end
            end
        end
    end

    exp4=zeros(length(tissueModelF5.rxns),1);
    for q=1:length(tissueModelF5.rxns)
        if	ismember(tissueModelF5.rxns(q),model.rxns)
            exp4(q,1)=metagenomeset(find(strcmp(model.rxns,tissueModelF5.rxns(q))),3);
        end
    end

    AA1=tissueModelF5.rxns(find(exp4<1 & exp4> -1));
    BB1=[];
    for q=1:length(AA1)
        if	ismember(AA1(q),model.rxns)
            BB1(q,1)=metagenomeset(find(strcmp(model.rxns,AA1{q})),3);
        else 
            BB1(q,1)=1;
        end
    end


    T=table(AA1,BB1);
    Sort_Table = sortrows(T,'BB1');
    AA2=table2cell(Sort_Table(:,1));
    %prune the reactions if they dont collapse the network and functionality of the model
    fba=optimizeCbModel(tissueModelF5);
    index=1:50:length(AA2);
    for h1=1:length(index)
        if h1 ~= length(index)
            [~, ~, removedRxns] = removeDeadEnds(tissueModelF5);
            d=AA2(index(h1):(index(h1)+49));
            T1=tissueModelF5;
            T1= removeRxns(T1,d);
            [~, ~, removedRxns1] = removeDeadEnds(T1);
            removed=setdiff(removedRxns1,removedRxns);
            hhh=intersect(removed,FromMetaGenomics);
            fbaa=optimizeCbModel(T1);
            if  length(hhh)<10  && fbaa.f > fba.f*0.1
                tissueModelF5=T1;
            else
                for h2=1:length(d)
                    T1=tissueModelF5;
                    T1= removeRxns(T1,d{h2});
                    [~, ~, removedRxns1] = removeDeadEnds(T1);
                    removed=setdiff(removedRxns1,removedRxns);
                    hhh=intersect(removed,FromMetaGenomics);
                    fbaa=optimizeCbModel(T1);
                    if  length(hhh)<10  && fbaa.f > fba.f*0.1
                        tissueModelF5=T1;
                    end
                end
            end
        else
            [~, ~, removedRxns] = removeDeadEnds(tissueModelF5);
            d=AA2(index(h1):end);
            for h3=1:length(d)
                    T1=tissueModelF5;
                    T1= removeRxns(T1,d{h3});
                    [~, ~, removedRxns1] = removeDeadEnds(T1);
                    removed=setdiff(removedRxns1,removedRxns);
                    hhh=intersect(removed,FromMetaGenomics);
                    fbaa=optimizeCbModel(T1);
                    if  length(hhh)<10 && fbaa.f > fba.f*0.1
                        tissueModelF5=T1;
                    end
            end
        end
    end

    model_aerobe=tissueModelF5;
	
	reactionInModel=union(model_anaerobic.rxns,model_aerobe.rxns);
	removeRxnsF=setdiff(model.rxns,reactionInModel);
	contextSpecificModel= removeRxns(model,removeRxnsF)
	contextSpecificModel.lb(find(strcmp(contextSpecificModel.rxns, 'rxn08173')))=0;
    contextSpecificModel.ub(find(strcmp(contextSpecificModel.rxns, 'rxn08173')))=1000;
    if exist('produced','var')
        contextSpecificModel.lb(find(ismember(contextSpecificModel.rxns,produced)))=0;
        contextSpecificModel.ub(find(ismember(contextSpecificModel.rxns,produced)))=1000;
    end
end
    
end

function reducedModelTemp = KeepNecessaryRxns(model, score, threshold, min)
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
function reducedModelTemp = KeepNecessaryRxnsLikeGIMME(model, score, threshold, min)
[modelIrrev,~,rev2irrev,irrev2rev] = convertToIrreversible(model);
expressionRxnsIrrev = zeros(length(modelIrrev.rxns),1);
for i1=1:length(modelIrrev.rxns)
    expressionRxnsIrrev(i1,1) = score(irrev2rev(i1,1),1);
end
cc=optimizeCbModel(model);
modelIrrev.lb(find(modelIrrev.c ==1),1)=cc.f*min; % minimum fraction of objective 
modelIrrev.c(:,1)=0;
for i1=1:length(modelIrrev.rxns)
    if (expressionRxnsIrrev(i1,1) > -1)   %if not absent reaction
        if (expressionRxnsIrrev(i1,1) < threshold)
            modelIrrev.c(i1,1) = threshold-expressionRxnsIrrev(i1,1); %FIX: use expression level as weight
        end
    end
end
gimmeSolution = optimizeCbModel(modelIrrev,'min');
reactionScoreTransition=zeros(length(expressionRxnsIrrev),1);
if (gimmeSolution.stat ~= 1)
    reactionScoreTransition(:,1) = 0;
else
    reactionScoreTransition(find(gimmeSolution.x>0),1)=1;
end
reactionScoreTransition(find(expressionRxnsIrrev>threshold))=1;

    %Translate reactionActivity to reversible model
    reactionActivity = zeros(nRxns,1);
    for i=1:nRxns
        for j=1:size(rev2irrev{i,1},2)
            if (reactionScoreTransition(rev2irrev{i,1}(1,j)) > reactionActivity(i,1))
                reactionActivity(i,1) = reactionScoreTransition(rev2irrev{i,1}(1,j));
            end
        end
    end 
    remove = model.rxns(reactionActivity == 0);
    reducedModelTemp = removeRxns(model,remove); 
end