function [constrainedModel]= DietConstrain(model,dietOption)
% this function constrain GEM using the dietOption.
%inputs:
%   model:				metabolic Model with COBRA or RAVEN format,
%   dietOption:			1 to 5 (% 1:high Fibre Plant Based, 2:high Fibre omnivore, 3:high Protein Plant based
%						4:high protein omnivore, 5:UK average.)
%outputs: 
%   constrainedModel:   constrained metabolic Model

%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk

%type of diet
index=dietOption;
% get path to where the MIGRENE Toolbox is located
MIGDIR = fileparts(which('MIGRENE_pipeline'));
% provide the path to the diets
load([MIGDIR filesep 'mat' filesep 'diets.mat']);

for i =1:length(diets.rxn)
    value= -(diets.value(i,index));
    if value ~= 0
        model.lb(find(strcmp(model.rxns, diets.rxn{i})))=value;
    end
end
constrainedModel=model;
end