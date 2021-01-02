function [catalogFileChecked]= checkCatalog(inputFile,numWorkers)
%inputs:
%   inputFile:             cell array contains two columns, first column is gene name and 
                            % second column provides KO annotation
%   numWorkers             integer indicating the number of cores to use for parallelization

%outputs: 
%   catalogFileChecked:    cell array contains two columns, first column is gene name and 
                             % second column provides KO annotation
%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk
 
% check the availability of KO annotation for each gene

expression = 'K\d\d\d\d\d';
matchStr = regexp(inputFile(:,2),expression,'match');
Index = find(not(cellfun('isempty',matchStr)));
if isempty(Index)
    error('the catalog is not annotated by KO or if KO-annotated, check the format of inputFile')
end

% genes can be annotated to more than one KO. Here, it splits the 
% annotation column and rearrange the catalog to long format
% including repeated genes with one KO annotation in each row 

% find the rows with one and more than one KO annotation 
index=strfind(inputFile(:,2),'K');
NumberOfKO = cellfun('length',index);
IndexAboveOne= find(NumberOfKO > 1);
IndexOne= find(NumberOfKO == 1);

Index = find(not(cellfun('isempty',index)));
if ~isempty(Index)
    % make a subset of gene catalog including rows with more than one KO annotation .
    if ~isempty(IndexAboveOne)
    output_1=inputFile(IndexAboveOne,:);
    end
    % make a subset of gene catalog including rows with one KO linked. 
     if ~isempty(IndexOne)
    output_2=inputFile(IndexOne,:);
     end
end

% check the number of workers for parallelization 
if numWorkers > 1
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers);
    end
else
    disp('You didnot specify the number of workers, so parallel mode is disabled. please dedicate number of workers')
end
% convert the subset of gene catalog with several annotated genes to long format. 
if ~isempty(IndexAboveOne)
    transition2={};
    parfor w=1:size(output_1,1)
        transition1=output_1{w,2};
        expression = 'K\d\d\d\d\d';
        matchStr = regexp(transition1,expression,'match');
        transition2{w,1}=matchStr;
    end
elseif ~isempty(IndexOne)
    catalogFileChecked=output_2;
end
output_1updated=cell([0 2]);
tic
if ~isempty(IndexAboveOne)
    parfor w=1:size(transition2,1)
        transition3={};
        transition1=transition2{w};
        transition3(:,2)=transition1';
        transition3(:,1)=output_1(w,1);
        output_1updated=vertcat(output_1updated,transition3);
    end
end
toc

if ~isempty(IndexAboveOne) & ~isempty(IndexOne)
% Concatenate the two arrays vertically to make a catalog file
catalogFileChecked=vertcat(output_2,output_1updated);
elseif ~isempty(IndexAboveOne) & isempty(IndexOne)
  catalogFileChecked=output_1updated;  
end