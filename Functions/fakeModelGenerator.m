function [fakemodel, fakemodel1]=fakeModelGenerator(mets_art,Ex_art)
%inputs:
%   mets_art: 			list of exchange metabolites and suffix of the compartments
%   Ex_art:				list of exchange reactions and suffix of the compartments
%outputs:
%   fakemodel: 			a metabolic model
%   fakemodel1: 		a metabolic model
%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk
fakemodel.mets={mets_art{:,2} mets_art{:,3} mets_art{:,4}}';
fakemodel.rxns={Ex_art{:,2} Ex_art{:,3} Ex_art{:,4} Ex_art{:,5}}';
fakemodel.S=zeros(length(fakemodel.mets),length(fakemodel.rxns));
fakemodel.lb=zeros(length(fakemodel.rxns),1);
fakemodel.ub=zeros(length(fakemodel.rxns),1);
fakemodel.lb(find(ismember(fakemodel.rxns,Ex_art(:,2))))=-1000;
fakemodel.ub(find(ismember(fakemodel.rxns,Ex_art(:,2))))=0;

fakemodel.lb(find(ismember(fakemodel.rxns,Ex_art(:,3))),1)=0;
fakemodel.ub(find(ismember(fakemodel.rxns,Ex_art(:,3))),1)=1000;

fakemodel.lb(find(ismember(fakemodel.rxns,Ex_art(:,4))),1)=0;
fakemodel.ub(find(ismember(fakemodel.rxns,Ex_art(:,4))),1)=100000;

fakemodel.lb(find(ismember(fakemodel.rxns,Ex_art(:,5))),1)=0;
fakemodel.ub(find(ismember(fakemodel.rxns,Ex_art(:,5))),1)=100000;

for w1=1:size(mets_art,1)
    fakemodel.S(find(strcmp(fakemodel.mets,mets_art{w1,3})),find(strcmp(fakemodel.rxns,Ex_art{w1,2})))=-1;
    fakemodel.S(find(strcmp(fakemodel.mets,mets_art{w1,3})),find(strcmp(fakemodel.rxns,Ex_art{w1,3})))=-1;
    fakemodel.S(find(strcmp(fakemodel.mets,mets_art{w1,2})),find(strcmp(fakemodel.rxns,Ex_art{w1,3})))=1;
    
    fakemodel.S(find(strcmp(fakemodel.mets,mets_art{w1,2})),find(strcmp(fakemodel.rxns,Ex_art{w1,4})))=-1;
    fakemodel.S(find(strcmp(fakemodel.mets,mets_art{w1,4})),find(strcmp(fakemodel.rxns,Ex_art{w1,4})))=1;
    
    fakemodel.S(find(strcmp(fakemodel.mets,mets_art{w1,4})),find(strcmp(fakemodel.rxns,Ex_art{w1,5})))=-1;
end

fakemodel1.mets={mets_art{:,1} mets_art{:,2}}';
fakemodel1.rxns=Ex_art(:,1);
fakemodel1.S=zeros(length(fakemodel1.mets),length(fakemodel1.rxns));
fakemodel1.lb=zeros(length(fakemodel1.rxns),1);
fakemodel1.ub=zeros(length(fakemodel1.rxns),1);

fakemodel1.lb(find(ismember(fakemodel1.rxns,Ex_art(:,1))),1)=-1000;
fakemodel1.ub(find(ismember(fakemodel1.rxns,Ex_art(:,1))),1)=1000;

for w1=1:size(mets_art,1)  
    fakemodel1.S(find(strcmp(fakemodel1.mets,mets_art{w1,2})),find(strcmp(fakemodel1.rxns,Ex_art{w1,1})))=1;
    fakemodel1.S(find(strcmp(fakemodel1.mets,mets_art{w1,1})),find(strcmp(fakemodel1.rxns,Ex_art{w1,1})))=-1;
end
end

