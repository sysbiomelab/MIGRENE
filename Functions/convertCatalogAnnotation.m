function [catalogFileConverted]= convertCatalogAnnotation(inputFile,mapping,numWorkers)
% this function automatically download and convert KO to KEGG reaction IDs in catalog.
%inputs:
%   inputFile:             cell array contains two columns, first column is gene name and 
%                          second column provides KO annotation
% OPTIONAL INPUTS:
%   mapping:               either cell array contains KOs (column 1) and their corresponding kegg 
%                          reaction ID (column 2) or empty cell array

%outputs: 
%   catalogFileConverted:   cell array contains two columns, first column is gene name and 
%                           second column provides kegg reaction annotation
%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk

% check the number of workers for parallelization 
if numWorkers > 1
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers);
    end
else
    disp('You didnot specify the number of workers, so parallel mode is disabled. please dedicate number of workers')
end
% if the mapping file for converting KO to reaction ID is not provided, it
% automatically download it from KEGG API. 
if isempty(mapping)
%   get path to where the MIGRENE Toolbox is located
    MIG = fileparts(which('MIGRENE_pipeline'));
    DATADIR=[MIG filesep 'data'];
%   save ko2rn file in the directory "data" of MIGRENE Toolbox
    urlwrite('http://rest.kegg.jp/link/reaction/ko',[DATADIR filesep 'ko2rn.txt']);
    Transition1 = readtable([DATADIR filesep 'ko2rn.txt'],'Format','%s%s','ReadVariableNames', false);
    ko2rn=table2cell(Transition1) ;
else
    ko2rn=mapping;
end

%remove prefix ko: and rn: from ko2rn if present 
ko2rn(:,1)=strrep(ko2rn(:,1),'ko:','');
ko2rn(:,2)=strrep(ko2rn(:,2),'rn:','');

%group ko2rn by KO so that the rn ID for the same KO were summerized. 
[~,~,ind]=unique(ko2rn(:,1));
uni=unique(ind);

File={};
for w=1:numel(uni)
    Index=find(ind == uni(w));
    transition2=unique(ko2rn(Index,2));
    str =strjoin(transition2,',');
    File(w,1)=unique(ko2rn(Index,1));
    File{w,2}=str;
end
inputFile1=inputFile;
for w=1:size(inputFile,1)
    transition1=inputFile{w,2};
    connection=File(find(strcmp(File(:,1),transition1)),2);
    if ~isempty(connection)
    inputFile1(w,2)=connection;
    end
end

expression = 'K\d\d\d\d\d';
matchStr = regexp(inputFile1(:,2),expression,'match');
Index = find(cellfun('isempty',matchStr));
catalogFileConverted=inputFile1(Index,:);
end

