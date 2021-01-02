function [microbiomeGEM]= microbiomeGEMgeneration(genericModel,cataloginfo,annotationFile,numWorkers)
% this function automatically download and convert KO to KEGG reaction IDs in your catalog.
%inputs:
%   genericModel:        metabolic Model with COBRA or RAVEN format, 
%   cataloginfo:          cell array contains two columns, first column is gene name and 
%                        second column provides kegg reaction ID
% OPTIONAL INPUTS:
%   annotationFile:       either cell array contains reaction IDs of the model(column 1),
%                         and column(s) for KO, EC and kegg RN annotations or empty cell array
%   numWorkers:           integer indicating the number of cores to use for parallelization

%outputs: 
%   microbiomeGEM:   cell array contains two columns, first column is gene name and 
%                           second column provides kegg reaction anntation
%#Author: Gholamreza Bidkori, KCL, UK, email: gbidkhori@gmail.com, gholamreza.bidkhori@kcl.ac.uk

if nargin<4
    numWorkers=1;
end
% get path to where the MIGRENE Toolbox is located
MIG = fileparts(which('MIGRENE_pipeline'));
DATADIR=[MIG filesep 'data'];
% check the number of workers for parallelization 
if numWorkers > 1
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(numWorkers);
    end
else
    disp('You didnot specify the number of workers, so parallel mode is disabled. please dedicate number of workers')
end
%% 
% if annotationFile is not provided, it automatically finds the annotations in the model. 
if isempty(annotationFile)
    disp('collecting the annotations from the model')
    %find fields including KO, rn or EC IDs
    KOexpression = 'K\d\d\d\d\d'; KOs={};
    RNexpression = 'R\d\d\d\d\d'; RNs={};
    ECexpression = 'EC:\d.\d.'; ECs={};
    fnames = fieldnames(genericModel);
    for w=1:numel(fnames)
        transition=strcat('genericModel.',fnames(w));
        transition=eval(char(transition));
        if iscellstr(transition)
            % find fields providing KO annotation for the reactions
            matchStr = regexp(transition,KOexpression,'match');
            Index = find(not(cellfun('isempty',matchStr)));
            if length(Index)/length(genericModel.rxns)*100>1
                KOs= vertcat(KOs,fnames(w));
            end
            % find fields providing kegg reaction IDs annotation for the reactions
            matchStr = regexp(transition,RNexpression,'match');
            Index = find(not(cellfun('isempty',matchStr)));
            if length(Index)/length(genericModel.rxns)*100>1
                RNs= vertcat(RNs,fnames(w));
            end
            % find fields providing EC annotation for the reactions
            matchStr = regexp(transition,ECexpression,'match');
            Index = find(not(cellfun('isempty',matchStr)));
            if length(Index)/length(genericModel.rxns)*100>1
                ECs= vertcat(ECs,fnames(w));
            end
        end
    end
else
   %find fields including KO, rn or EC IDs in provided annotationFile
   disp('collecting the annotations from the annotationFile')
   disp('as annotationFile is provided, the function ignores the available annotations in the model.')
   disp('please make sure annotationFile follows the following format:')
   disp('cell array contains reaction IDs of the model(column 1),and column(s) for KO, EC and kegg RN annotations')
    KOexpression = 'K\d\d\d\d\d'; KOs=[];
    RNexpression = 'R\d\d\d\d\d'; RNs=[];
    ECexpression = 'EC:\d.\d.'; ECs=[];
    numberOfannotation=size(annotationFile,2);
    for w=2:numberOfannotation
        transition=annotationFile(:,w);
        if iscellstr(transition)
            % find fields providing KO annotation for the reactions
            matchStr = regexp(transition,KOexpression,'match');
            Index = find(not(cellfun('isempty',matchStr)));
            if length(Index)/length(genericModel.rxns)*100>1
                KOs= [KOs w];
            end
            % find fields providing kegg reaction IDs annotation for the reactions
            matchStr = regexp(transition,RNexpression,'match');
            Index = find(not(cellfun('isempty',matchStr)));
            if length(Index)/length(genericModel.rxns)*100>1
                RNs=  [RNs w];
            end
            % find fields providing EC annotation for the reactions
            matchStr = regexp(transition,ECexpression,'match');
            Index = find(not(cellfun('isempty',matchStr)));
            if length(Index)/length(genericModel.rxns)*100>1
                ECs=  [ECs w];
            end
        end
    end
end

if length(RNs)>1 && isempty(annotationFile)
    error('there are more than one field in the metabolic model including kegg reaction IDs. It must be one cell array with kegg reaction ID') 
end

if length(ECs)>1 && isempty(annotationFile)
    error('there are more than one field in the metabolic model including EC annotation. It must be one cell array with EC annotation') 
end

if length(KOs)>1 && isempty(annotationFile)
    error('there are more than one field in the metabolic model including KO annotation. It must be one cell array with KO annotation') 
end

if length(RNs)>1 && ~isempty(annotationFile)
    error('there are more than one column in annotationFile including kegg reaction IDs. It must be one cell array with kegg reaction ID') 
end

if length(ECs)>1 && ~isempty(annotationFile)
    error('there are more than one column in annotationFile including EC annotation. It must be one cell array with EC annotation') 
end

if length(KOs)>1 && ~isempty(annotationFile)
    error('there are more than one column in annotationFile including KO annotation. It must be one cell array with KO annotation') 
end

if isempty(KOs)
    if isempty(RNs)
        if isempty(ECs)
            error('there are not any info in the metabolic model or annotationFile for KO,EC or kegg reaction ID annotation. It must be at least one cell array with one of the annotations') 
        else
            disp('EC annotation was found for integration of genes in catalog into model')
        end
    else
        if isempty(ECs)
            disp('kegg reaction annotation was found for integration of genes in catalog into model')
        else
            disp('kegg reaction and EC annotation were found for integration of genes in catalog into model')
        end
    end
else
     if isempty(RNs)
        if isempty(ECs)
            disp('KO annotation was found for integration of genes in catalog into model') 
        else
            disp('EC and KO annotations were found for integration of genes in catalog into model')
        end
    else
        if isempty(ECs)
            disp('KO and kegg reaction annotations were found for integration of genes in catalog into model')
        else
            disp('KO, kegg reaction and EC annotation were found for integration of genes in catalog into model')
        end
    end
end

% make the model ready for integration
%collect all the annotations in the model in one temporary field
genericModel.temporary=cell([numel(genericModel.rxns) 1]);
if ~isempty(RNs) && isempty(annotationFile)
    transition=strcat('genericModel.',RNs(1));
    genericModel.temporary=eval(char(transition));
elseif ~isempty(RNs) && ~isempty(annotationFile)
    for q=1:numel(genericModel.rxns)
      if ismember(genericModel.rxns{q},annotationFile(:,1))      
        genericModel.temporary(q,1)=annotationFile(find(strcmp(model.rxns,genericModel.rxns{q})),RNs);
      end
    end
end

if ~isempty(ECs)
    urlwrite('http://rest.kegg.jp/link/reaction/enzyme',[DATADIR filesep 'ec2rn.txt']);
    Transition1 = readtable([DATADIR filesep 'ec2rn.txt'],'Format','%s%s','ReadVariableNames', false);
    ec2rn=table2cell(Transition1);
end
if ~isempty(KOs)
    urlwrite('http://rest.kegg.jp/link/reaction/ko',[DATADIR filesep 'ko2rn.txt']);
    Transition1 = readtable([DATADIR filesep 'ko2rn.txt'],'Format','%s%s','ReadVariableNames', false);
    ko2rn=table2cell(Transition1) ;
end
check=who;

if ismember('ko2rn',check)
    ko2rn(:,1)=strrep(ko2rn(:,1),'ko:','');
    ko2rn(:,2)=strrep(ko2rn(:,2),'rn:','');
    %find not annotated reactions in the temporary file and fill in regarding available KO
    index=strfind(genericModel.temporary,'R');
    Index = find(cellfun('isempty',index));
    if isempty(annotationFile)
        transition=strcat('genericModel.',KOs(1));
        koInfo=eval(char(transition));
        koInfo=strrep(koInfo,' ',',');
        koInfo=strrep(koInfo,';',',');
        for qq=1:length(Index)     
            ec=char(koInfo(Index(qq)));
            if ~isempty(ec)
                tra1=strsplit(ec,',')';
                expression = 'K\d\d\d\d\d';
                matchStr = regexp(tra1,expression,'match');
                tra1=tra1(find(not(cellfun('isempty',matchStr))));
                anno={};
                if ~isempty(tra1)
                    for ee=1:length(tra1)
                         if ismember(tra1{ee},ec2rn(:,1))
                             % get all the kegg reaction IDs related for the KO
                             dd=ec2rn(ismember(ec2rn(:,1),tra1{ee}),2);
                             anno=vertcat(anno,dd);
                         end
                    end
                    genericModel.temporary{Index(qq),1}=strjoin(unique(anno),',');
                end
            end
        end
    else
        koInfo=annotationFile(:,KOs);
        koInfo=strrep(koInfo,' ',',');
        koInfo=strrep(koInfo,';',',');
        for qq=1:length(Index)
			tran1=koInfo(find(strcmp(annotationFile(:,1),genericModel.rxns(Index(qq)))));
            ec=char(tran1);
            if ~isempty(ec)
                tra1=strsplit(ec,',')';
                expression = 'K\d\d\d\d\d';
                matchStr = regexp(tra1,expression,'match');
                tra1=tra1(find(not(cellfun('isempty',matchStr))));
                if ~isempty(tra1)
                    anno={};
                    for ee=1:length(tra1)
                         if ismember(tra1{ee},ec2rn(:,1))
                             dd=ec2rn(ismember(ec2rn(:,1),tra1{ee}),2);
                             anno=vertcat(anno,dd);
                         end
                    end
                    genericModel.temporary{Index(qq),1}=strjoin(unique(anno),',');
                end
            end
        end 
    end
end

   
if ismember('ec2rn',check)
    %remove prefix rn: from ko2rn if present 
    ec2rn(:,2)=strrep(ec2rn(:,2),'rn:','');
    %find not annotated reactions in the temporary file and fill in regarding available EC
    index=strfind(genericModel.temporary,'R');
    Index = find(cellfun('isempty',index));
    if isempty(annotationFile)
        transition=strcat('genericModel.',ECs(1));
        ecInfo=eval(char(transition));
        ecInfo=strrep(ecInfo,' ',',');
        ecInfo=strrep(ecInfo,';',',');
        ecInfo=strrep(ecInfo,'EC','ec');
        for qq=1:length(Index)     
            ec=char(ecInfo(Index(qq)));
            if ~isempty(ec)
                tra1=strsplit(ec,',')';
                expression = 'ec:\d.\d.';
                matchStr = regexp(tra1,expression,'match');
                tra1=tra1(find(not(cellfun('isempty',matchStr))));
                anno={};
                if ~isempty(tra1)
                    for ee=1:length(tra1)
                         if ismember(tra1{ee},ec2rn(:,1))
                             % get all the kegg reaction IDs related for the ec
                             dd=ec2rn(ismember(ec2rn(:,1),tra1{ee}),2);
                             anno=vertcat(anno,dd);
                         end
                    end
                    genericModel.temporary{Index(qq),1}=strjoin(unique(anno),',');
                end
            end
        end
    else
        ecInfo=annotationFile(:,ECs);
        ecInfo=strrep(ecInfo,' ',',');
        ecInfo=strrep(ecInfo,';',',');
        ecInfo=strrep(ecInfo,'EC','ec');
        for qq=1:length(Index)
			tran1=ecInfo(find(strcmp(annotationFile(:,1),genericModel.rxns(Index(qq)))));
            ec=char(tran1);
            if ~isempty(ec)
                tra1=strsplit(ec,',')';
                expression = 'ec:\d.\d.';
                matchStr = regexp(tra1,expression,'match');
                tra1=tra1(find(not(cellfun('isempty',matchStr))));
                if ~isempty(tra1)
                    anno={};
                    for ee=1:length(tra1)
                         if ismember(tra1{ee},ec2rn(:,1))
                             dd=ec2rn(ismember(ec2rn(:,1),tra1{ee}),2);
                             anno=vertcat(anno,dd);
                         end
                    end
                    genericModel.temporary{Index(qq),1}=strjoin(unique(anno),',');
                end
            end
        end 
    end
end

%% summerize the catalog based on kegg reaction IDs
disp('start getting the catalog info for integration')
cataloginfoLong=ConvertTOLongFormat(cataloginfo,numWorkers);
[~,~,ind]=unique(cataloginfoLong(:,2));
uni=unique(ind);
catalogForInteg={};
for w=1:numel(uni)
    Index=find(ind == uni(w));
    transition2=unique(cataloginfoLong(Index,1));
    str =strjoin(transition2,',');
    catalogForInteg(w,1)=unique(cataloginfoLong(Index,2));
    catalogForInteg{w,2}=str;
end

%% add GPRs in the model based on catalog data and generate the GEM
disp('start adding GPR into the model')
% this section fill all the following empty cell arrays
genericModel.grRules=cell([numel(genericModel.rxns) 1]);
genericModel.genes={};
genericModel.rules=cell([numel(genericModel.rxns) 1]);
genericModel.geneNames={};
genericModel.rxnGeneMat={};

for w=1:numel(genericModel.temporary)
    rnInfo=genericModel.temporary{w,1};
    if ~isempty(rnInfo)
        rnInfo=strrep(rnInfo,' ',',');
        rnInfo=strrep(rnInfo,';',',');
        tra1=strsplit(rnInfo,',')';
        expression = 'R\d\d\d\d\d';
        matchStr = regexp(tra1,expression,'match');
        tra1=tra1(find(not(cellfun('isempty',matchStr))));
        if ~isempty(tra1)
            annotation={};
            for ee=1:length(tra1)
                 if ismember(tra1{ee},catalogForInteg(:,1))
                     dd=catalogForInteg(ismember(catalogForInteg(:,1),tra1{ee}),2);
                     annotation=vertcat(annotation,dd);
                 end
            end
            genericModel.grRules(w,1)={strjoin(unique(annotation),',')};
        end
    else
        genericModel.grRules(w,1)={''};    
    end
end

temp1={};
for i=1:length(genericModel.grRules)
    if ~isempty(genericModel.grRules{i})
     tra2=strsplit(genericModel.grRules{i},',')';
     temp1=vertcat(temp1,tra2);
    end
end
genericModel.genes=unique(temp1); % genericModel.genes generated
genericModel.geneNames=genericModel.genes; % genericModel.geneNames generated

c=[];
for h=1:length(genericModel.grRules)
    if  ~isempty(genericModel.grRules{h,1})
        tra2=strsplit(genericModel.grRules{h,1},',')';
        for k=1:length(tra2)
        index=find(strcmp(genericModel.genes,tra2{k}));        
        c=[c;[h index 1]];
        end
    end
end
checkTheSize=intersect(find(c(:,1)== length(genericModel.rxns)),find(c(:,2)== length(genericModel.genes)));
if isempty(checkTheSize)
    c=[c;[length(genericModel.rxns) length(genericModel.genes) 0]];
end
temp1= sparse(c(:,1),c(:,2),c(:,3));
genericModel.rxnGeneMat=temp1; % genericModel.rxnGeneMat generated


for h=1:length(genericModel.grRules)
    if  ~isempty(genericModel.grRules{h,1})
        tra2=strsplit(genericModel.grRules{h,1},',')';
        collect={};
        for k=1:length(tra2)
            index=find(strcmp(genericModel.genes,tra2{k})) ;
            converted=['x(' num2str(index) ')'];
            collect=vertcat(collect,converted);
        end
        genericModel.rules(h,1)={strjoin(unique(collect),',')}; % generating genericModel.rules
    else
        genericModel.rules(h,1)={''}; 
    end
end

for i=1:length(genericModel.rules)
    if ~isempty(regexp(genericModel.rules{i},',','match'))
     genericModel.rules{i}=['( ' genericModel.rules{i} ' )'];
    end
end
genericModel.rules=strrep(genericModel.rules,',',' | '); % genericModel.rules generated


for i=1:length(genericModel.grRules)
    if ~isempty(regexp(genericModel.grRules{i},',','match'))
     genericModel.grRules{i}=['(' genericModel.grRules{i} ')'];
    end
end
genericModel.grRules=strrep(genericModel.grRules,',',' or '); % genericModel.grRules generated

genericModel.rxnAnnotation=genericModel.temporary; %consider temporary file as extended annotation file for the model
genericModel = rmfield( genericModel,'temporary'); 

microbiomeGEM=genericModel;
end
%%
function [cataloginfoLong]= ConvertTOLongFormat(cataloginfo,numWorkers)
% find the rows with one and more than one kegg rn annotation 
index=strfind(cataloginfo(:,2),'R');
NumberOfKO = cellfun('length',index);
IndexAboveOne= find(NumberOfKO > 1);
IndexOne= find(NumberOfKO == 1);

Index = find(not(cellfun('isempty',index)));
if ~isempty(Index)
    % make a subset of gene catalog including rows with more than one KO annotation .
    if ~isempty(IndexAboveOne)
    output_1=cataloginfo(IndexAboveOne,:);
    end
    % make a subset of gene catalog including rows with one KO linked. 
     if ~isempty(IndexOne)
    output_2=cataloginfo(IndexOne,:);
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
        expression = 'R\d\d\d\d\d';
        matchStr = regexp(transition1,expression,'match');
        transition2{w,1}=matchStr;
    end
elseif ~isempty(IndexOne)
    cataloginfoLong=output_2;
end
output_1updated=cell([0 2]);
%tic
if ~isempty(IndexAboveOne)
    parfor w=1:size(transition2,1)
        transition3={};
        transition1=transition2{w};
        transition3(:,2)=transition1';
        transition3(:,1)=output_1(w,1);
        output_1updated=vertcat(output_1updated,transition3);
    end
end
%toc

if ~isempty(IndexAboveOne) & ~isempty(IndexOne)
% Concatenate the two arrays vertically to make a catalog file
cataloginfoLong=vertcat(output_2,output_1updated);
elseif ~isempty(IndexAboveOne) & isempty(IndexOne)
  cataloginfoLong=output_1updated;  
end
end
