%% makeSimulated Sheet.m
% Created 12-Jul-2022 WN: agonist/antagonist and gene columns, gene mapping
% Edited 13-Jul-2022 WN: add index mapping, handle dual agonist/antagonist
% case

%% housekeeping
close all;
clear
clc

%% file read-ins
fileDir = 'drugbank_resultsComb.xlsx';
tbl = readtable(fileDir);

% parse data
tb1 = tbl(:,[1,7:10]);

%% create Y/N columns
isAgonist = strcmp('Agonist',tb1.TargetActions);
isAntagonist = strcmp('Antagonist',tb1.TargetActions);
isBoth = strcmp('Both',tb1.TargetActions);

% agregate both into agonist/antagonist
isAgonist = any([isAgonist,isBoth],2);
isAntagonist = any([isAntagonist,isBoth],2);

% set Y/N
Agonist(isAgonist,1) = "Yes";
Agonist(~isAgonist,1) = "No";
Antagonist(isAntagonist,1) = "Yes";
Antagonist(~isAntagonist,1) = "No";

%% Create target columns
AgonistTargetGene = tb1.TargetGeneName;
AntagonistTargetGene = tb1.TargetGeneName;
AgonistTargetGene(~isAgonist) = {''};
AntagonistTargetGene(~isAntagonist) = {''};

% handle both case
tmp = cellfun(@(x) split(x,'/'),tb1.Notes(isBoth),'UniformOutput',0);
for r = 1:size(tmp,1)
    tmp1 = tmp{r};
    tmp(r,1) = tmp1(1);
    tmp(r,2) = tmp1(2);
end
AgonistTargetGene(isBoth) = tmp(:,1);
AntagonistTargetGene(isBoth) = tmp(:,2);

% split to each type and store in column

%% grab Gene IDs
fileDir = 'F:\Creative Inquiry (BIOE)\Models\FibroblastSNM_20220623.xlsx';
tbl = readtable(fileDir,'Sheet','species');
tbl = [tbl.ID,tbl.geneName];

%% get gene names
addpath('F:\Creative Inquiry (BIOE)\Netflux Model');
[params,~] = NetfluxODE_loadParams();
speciesNames = params{1,4};

% remove empty rows
blank = cellfun(@isempty,tbl(:,2));
tbl = tbl(~blank,:);

%% Map Creation/Processing
M = containers.Map;
I = containers.Map;
sz = size(tbl,1);

for r = 1:sz
    % get species name => index map
    idx = find(strcmp(tbl{r,1},speciesNames));
    I(tbl{r,1}) = idx;

    % get gene name => species name map
    gene_key = split(tbl{r,2},';');
    for i_k = 1:length(gene_key)
        if isKey(M, gene_key{i_k})
            continue
        end
        M(gene_key{i_k}) = tbl{r,1};
    end
end

%% Target Gene and Gene Index creation
sz = size(tb1,1);
AgonistIndex = repmat({''},sz,1);
AntagonistIndex = AgonistIndex;

for r = 1:sz
    if strcmp(Agonist(r),'Yes')
        genes = split(AgonistTargetGene(r),';');
        specs = values(M,genes);
        AgonistTargetGene(r) = {join(unique(specs),';')};
        idx = string(values(I,specs));
        AgonistIndex(r) = {join(unique(idx),';')};
    end
    if strcmp(Antagonist(r),'Yes')
        genes = split(AntagonistTargetGene(r),';');
        specs = values(M,genes);
        AntagonistTargetGene(r) = {join(unique(specs),';')};
        idx = string(values(I,specs));
        AntagonistIndex(r) = {join(unique(idx),';')};
    end
end

%% Create Agonist Target Columns given Gene ID
% split each row of Agonist target genes and get their index in the conversion table
geneNames = cellfun(@(x) split(x,';'),AgonistTargetGene,'UniformOutput',0);
[r,~] = cellfun(@(x) cellfun(@(y)find(string(y)==tbl(:,1)),x,'UniformOutput',0),geneNames,'UniformOutput',0);
%process the row index to get the Gene Names
r = cellfun(@(x) cell2mat(x),r,'UniformOutput',0);
idNames = cellfun(@(x) tbl(x,2),r,'UniformOutput',0);
% join the gene names by row
AgonistGeneName = cellfun(@(x)join(x,';'),idNames);

%% Create Antagonist Target Columns given Gene ID
% split each row of Antagonist target genes and get their index in the conversion table
geneNames = cellfun(@(x) split(x,';'),AntagonistTargetGene,'UniformOutput',0);
[r,~] = cellfun(@(x) cellfun(@(y)find(string(y)==tbl(:,1)),x,'UniformOutput',0),geneNames,'UniformOutput',0);
% process the row index to get the Gene Names
r = cellfun(@(x) cell2mat(x),r,'UniformOutput',0);
idNames = cellfun(@(x) tbl(x,2),r,'UniformOutput',0);
% join the gene names by row
AntagonistGeneName = cellfun(@(x)join(x,';'),idNames);

%% Target Action search to convert Name to Competitive/NonCompetitive
targetAct = repmat({''},sz,1);
actSearch = readtable("Drug Actions.xlsx");
[r,~] = cellfun(@(x) find(x==string(actSearch.Name)),tb1.Name);
targetAct = table2cell(actSearch(r,3));

%% Create Finalized Table
finalTbl = [tb1.Name,Agonist,AgonistTargetGene,AgonistGeneName,AgonistIndex,...
    Antagonist,AntagonistTargetGene,AntagonistGeneName,AntagonistIndex,...
    targetAct,tb1.SimilarDrugs];

colNames = {'Drug','IsAgonist','AgonistTargetGeneID','AgonistTarget',...
    'AgonistTargetIndex','IsAntagonist','AntagonistTargetGeneID',...
    'AntagonistTarget','AntagonistTargetIndex','DrugAction','SimilarDrugs'};

simulatedDrugs = array2table(finalTbl,'VariableNames',colNames);

% Clean up missing Drug actions
simulatedDrugs(simulatedDrugs.DrugAction=="",:)=[];
simulatedDrugs(simulatedDrugs.DrugAction=="NA",:)=[];

writetable(simulatedDrugs,'SimulatedDrugsUpdate.xlsx','WriteMode','overwritesheet','AutoFitWidth',0);
