%% tblParse.m
% Created 20-Jun-2022 WN: parse files and write to xlsx
% Edited 21-Jun-2022 WN: simplified targets and agonist/antagonist
% Edited 25-Jun-2022 WN: cellIndex to handle 2d cell, prevent SimilarDrugs
% column from overwrites
% Edited 27-Jun-2022 WN: remove same ATC Codes, prevent duplicates in SimilarDrugs
% Edited 01-Jul-2022 WN: incorperate competitive/non-competitive data
% Edited 13-Jul-2022 WN: change delimiters from '; ' to ';'
% Edited 13-Jul-2022 WN: Include notes column

%% housekeeping
close all;
clear
clc

%% file read-ins
fileDir = 'drugbank_results.csv';
tbl1 = readtable(fileDir);

fileDir = 'drugbank_results (1).csv';
tbl2 = readtable(fileDir);

fileDir = 'drugbank_results (2).csv';
tbl3 = readtable(fileDir);

fileDir = 'Drug Actions.xlsx';
actions = readtable(fileDir);

%% data parsing
masterTbl = [tbl1;tbl2;tbl3]; % combine csv files
masterTbl = masterTbl(:,[2:7,1,8]); % reorder columns
masterTbl = unique(masterTbl,'rows'); % Sort by drug name

masterTbl.TargetGeneName = cellfun(@(x) replace(x,' ',''),masterTbl.TargetGeneName,'UniformOutput',0);

% give Saucerman Drugs Target Action
masterTbl.TargetActions(cellIndex(masterTbl.Name,'Fasudil|Becaplermin|Trabedersen',1)) = {'antagonist'};
masterTbl.Description(cellIndex(masterTbl.Name,'Trabedersen',0))={'N/A'};

% give Saucerman Drugs ATC Codes
masterTbl.ATCCodes(cellIndex(masterTbl.Name,'Gallium nitrate',0)) = {'M05BA'};
masterTbl.ATCCodes(cellIndex(masterTbl.Name,'Marimastat',0)) = {'L01XX'};
masterTbl.ATCCodes(cellIndex(masterTbl.Name,'Trabedersen',0)) = {'L01XX.'};

%% get rows with competitive/non-competitive data
blank = cellIndex(table2cell(actions(:,1:3)),'',0);
actions = actions(~blank,:); % remove rows without data 

blank = cellIndex(table2cell(actions(:,1:3)),'NA',0);
actions = actions(~blank,:); % remove rows with NA

idx = ismember(masterTbl.Name,actions.Name);
masterTbl = masterTbl(idx,:); % remove drugs without action data


targets = readlines('targets.txt'); % read in model target nodes

sz = size(masterTbl,1);
gene = cell(sz);
for r = 1:sz
    %% remove targets not found in model
    names = split(string(masterTbl.TargetGeneName(r)),';');
    for i_name = 1:length(names)
        if ismember(names{i_name},targets)
            gene{r,i_name} = [names{i_name},';'];
        else
            gene{r,i_name} = '';
        end
    end

    geneName = [gene{r,:}];
    geneName = stripsemi(geneName);
    masterTbl.TargetGeneName(r) = {geneName};
    
    %% convert action column to competitive/non-competitive and agonist/antagonist
    idx = cellIndex(actions.Name,masterTbl.Name(r),0);
    tmp = append(actions.Competitive_Non_Competitive{idx},' ',actions.Agonist_Antagonist{idx});
    masterTbl.TargetActions(r) = {tmp};
    masterTbl.Notes(r) = actions.Notes(idx);
end

% get all similar ATC codes
sz = size(masterTbl,1);
M = uniqueMap(sz,masterTbl.ATCCodes);

% remove same ATC codes
masterTbl = removeSimilar(M,masterTbl);


% get all similar drugs - **Can't do until get competitive or non competitve
sz = size(masterTbl,1);
M = uniqueMap(sz,masterTbl.TargetGeneName,masterTbl.TargetActions);

% remove same drugs
masterTbl = removeSimilar(M,masterTbl);

tmp = cellfun(@(x) string(split(x)),masterTbl.TargetActions,'UniformOutput',0);
C_NC = cellfun(@(x) x(1),tmp);
masterTbl.Binding = C_NC;
Ag_Ant = cellfun(@(x) x(2),tmp);
masterTbl.TargetActions = Ag_Ant;
masterTbl = masterTbl(:,[1:8,10,9]);

writetable(masterTbl,'drugbank_resultsComb.xlsx','WriteMode','overwritesheet','AutoFitWidth',0);

%% strip all semicolons and spaces
function x = stripsemi(x)
    if isempty(x)||strcmp(x,'')
        return
    end
    x = char(x);
    test = strcmp(x(end),' ') || strcmp(x(end),';') || strcmp(x(1),' ') || strcmp(x(1),';');
    while test
        x = strip(strip(x),';');
        if isempty(x)
            return
        else
            test = strcmp(x(end),' ') || strcmp(x(end),';') || strcmp(x(1),' ') || strcmp(x(1),';');
        end
    end
end

%% create map of Key:{names} Value:{[indexes],...} where name is comb of [x,' ',y]
function M = uniqueMap(sz,x,y)
    if nargin==2
        y=repelem({''},sz);
    end
    for r = 1:sz
        names = split(string(x(r)),';');
        for i_name = 1:length(names)
                allcomb{r,i_name} = append(names{i_name},' ',y{r});
        end
    end

    M = containers.Map('KeyType','char','ValueType','any');
    
    [row,col] = size(allcomb);
    for r = 1:row
        for c = 1:col
            val = allcomb{r,c};
            if isempty(val)
                continue
            end
            if isKey(M,val)
                M(val) = [M(val),r];
            else
                M(val) = r;
            end
        end
    end
end

%% returns table t with similar rows as given by map M removed
function C = removeSimilar(M,t)
    newTbl={};
    names = repmat("",size(M));
    key = keys(M);
    for i_k = 1:length(key)
        idx = cell2mat(values(M,key(i_k)));
        str = string(t.Name(idx(2:end)));
        if isempty(str)
            str = "";
        end
        % prevent potential overwrite of SimilarDrugs column
        if ismember('SimilarDrugs', t.Properties.VariableNames)
            sim = string(t.SimilarDrugs(idx));
            sim = stripsemi(regexprep(join(sim,';'),';{2,}',';'));
            if ~isempty(sim)&&~strcmp(sim,'')
                str = [str;sim];
            end
        end
            
        % handle multiple strings when combining
        str1 = arrayfun(@(x) split(x,';'),str,UniformOutput=0);
        ind = cellfun(@(x) numel(x)>1,str1);
        if any(ind)
            dat = str1{ind};
            str1(ind)=[];
            str = [str1;dat];
        end
        names(i_k) = stripsemi(join(str,';'));
        newTbl = [newTbl;t(idx(1),:)];
        
    end
    newTbl.SimilarDrugs = names;
    
    % handle removal of duplicate rows
    [C,~,ic] = unique(newTbl(:,1:end-1),'rows');
    temp = strings(max(ic),1);
    for k = 1:length(ic)
        temp(ic(k)) = append(temp(ic(k)),';',string(newTbl{k,"SimilarDrugs"}));
        temp(ic(k)) = stripsemi(temp(ic(k)));
    end
    C.SimilarDrugs = temp;

    % remove duplicate SimilarDrugs
    for r = 1:size(C,1)
        C.SimilarDrugs(r) = join(unique(split(C.SimilarDrugs(r),';')),';');
    end
end

%% get index of cell rows
function idx = cellIndex(C,str,regex)
    idx = zeros(size(C));
    for i_col = 1:size(C,2)
        if regex
            idx(:,i_col) = ~cellfun(@(x) isempty(regexp(string(x),str, 'once')),C(:,i_col));
        else
            idx(:,i_col) = strcmp(string(C(:,i_col)),str);
        end
    end
    idx = any(idx,2);
end