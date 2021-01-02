function[micronutrients_diet_mmol, macronutrients_diet]= USDAcreatingDiet(food_id_item,food_grams_item)
% Input: 
%   food_grams_item:    the weight (in grams) for the food item 
%   food_id_item:       the ID of the food item (specified by USDA food ID)
%        NOTE: the food_grams_item and food_id_item need to be in the same order.
%        Also the food_grams_item needs to be specified in grams, so e.g 1.5 kg needs to be 1500 (g). 

%output: 
% macronutrients_diet:               total amount of macronutrients for the created diet 
% micronutrients_diet_mmol:          total amount of micronutrients for the created diet in mmol/gDW 

% #Authors Bouchra Ezzamouri.

% get path to where the MIGRENE Toolbox is located
MIGDIR = fileparts(which('MIGRENE_pipeline'));
%load USDA dataset from MIGRENE Toolbox.
USDA=[MIGDIR filesep 'mat' filesep 'USDAdataset.mat'];
load(USDA)

% find the common IDs in USDA dataset
food_id_members = ismember(food_item_USDA_id,food_id_item);%--> from the USDA foods (is a list of 8463 x 1) it will check if the specified input of the food_id_item is found the USDA food list. If so it is 1 otherwise 0  
food_macros = query_food_item_macros_values_1gDW(:,food_id_members==1); %--> from the list of food_id_members if it is equal to 1 (so the food is in the list) then we want the macro values from the matrix query_food_item_macros_values_1gDW
food_micros_mmolgDW = query_food_item_micros_mmol_gDW(:,food_id_members==1); % obtaining micronutrients from the matrix for the food of interest in mmol gDW

%output specifying matrix with a size of macronutrients/micronutrients x food items that were specified by the user. 
food_macros_diet = zeros(size(food_macros,1),size(food_id_item,1));
food_macros = transpose(food_macros);
food_macros_diet = transpose(food_macros_diet);

food_micros_mmolgDW_diet= zeros(size(food_micros_mmolgDW,1),size(food_id_item,1));
food_micros_mmolgDW = transpose(food_micros_mmolgDW);
food_micros_mmolgDW_diet= transpose(food_micros_mmolgDW_diet);

% for every food that the user wants in the diet multiply by the amount in
% grams that the user wants for that specific food. 
for i =1:length(food_id_item)
   food_macros_diet(i,:) = food_macros(i,:) * food_grams_item(i);
   food_micros_mmolgDW_diet(i,:) =  food_micros_mmolgDW(i,:) * food_grams_item(i)  ;
end

% the output will be a list of total macros in gDW and micronutrients (in gDW and mmol) for the specified diet 
 total_food_macros_diet = transpose(sum(food_macros_diet,1));
 total_food_micro_mmolgDW_diet= transpose(sum(food_micros_mmolgDW_diet,1));

 micronutrients_diet_mmol = table(mets_USDA_name,(total_food_micro_mmolgDW_diet));
 macronutrients_diet = table(macros_USDA_name,(total_food_macros_diet));
end
