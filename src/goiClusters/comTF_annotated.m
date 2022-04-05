%% 
% To find common TF that are expressed in the Transcription Cluster
% Input
% edgeTable with annotated gene expression from annoted_expressed.m
% Gene Expression Data, set threshold >0 or >2 
%
% Output 
% edgeTable with gene expression annotated & comTF expression annotated.

%% Load in Data
%addpath C:\Pore-C\01062022
%addpath C:\Pore-C\01142022\
addpath C:\Pore-C\GeneExpression\GenesExpressed
%addpath C:\Pore-C\01172022
%{
load('GM_ActiveGenes_Annotated.mat');
load("GM12878_Expressed_Genes.mat");
%}

%{
load('GM_Active_GenesAnnotated')
load("BJ_Genes_Expressed.mat");
load("IR_Genes_Expressed.mat");
load('GM_Active_GenesAnnotated')
load('IR_table_w_GE_annotated_01142022.mat');
load('BJ_table_w_GE_annotated_01142022.mat');
%load('BJClusters_ActiveGenesAnnotated.mat')
%load("IRClusters_ActiveGenesAnnotated.mat")
%}


%% Load BJ Single Cell Data
load('scBJ_GE_Annotated_03252022.mat')
addpath C:\Pore-C\GeneExpression\GenesExpressed
load("BJ_GenesTPMgreater0.mat");

%{
%% process IR
%IR=addvars(IR_edgeTable_w_genes_TF,IR_edgeTable_w_genes_TF.genes,'Before','genes','NewVariableNames',{'geneNetwork'});
idx_IR = ~cellfun(@isempty,IR.comTF); % to index by whether cell is empty or not
IR=IR(idx_IR,:);
IR.comTF=cellfun(@unique2cell,IR.comTF,'UniformOutput',false);
[active_genes_ir,active_expression_ir]=cellfun(@(x) find_active(x,IR_Genes_Expressed), IR.comTF,'UniformOutput',false);
IR.active_comTF=active_genes_ir;
IR.active_expression_comTF=active_expression_ir;
%}

%% process BJ
%BJ=addvars(BJ_edgeTable_w_genes_TF,BJ_edgeTable_w_genes_TF.genes,'Before','genes','NewVariableNames',{'geneNetwork'});
idx_BJ = ~cellfun(@isempty,scBJ.comTF); % to index by whether cell is empty or not
BJ=scBJ(idx_BJ,:);
BJ.comTF=cellfun(@unique2cell,scBJ.comTF,'UniformOutput',false);
[active_genes_bj,active_expression_bj]=cellfun(@(x) find_active(x,BJ_Genes_Expressed), BJ.comTF,'UniformOutput',false);
BJ.active_comTF=active_genes_bj;
BJ.active_expression_comTF=active_expression_bj;

%{
%% process GM
%BJ=addvars(BJ_edgeTable_w_genes_TF,BJ_edgeTable_w_genes_TF.genes,'Before','genes','NewVariableNames',{'geneNetwork'});
idx_GM = ~cellfun(@isempty,GM.comTF); % to index by whether cell is empty or not
GM=GM(idx_GM,:);
GM.comTF=cellfun(@unique2cell,GM.comTF,'UniformOutput',false);
[active_genes_gm,active_expression_gm]=cellfun(@(x) find_active(x,GM12878_GenesExpressed), GM.comTF,'UniformOutput',false);
GM.active_comTF=active_genes_gm;
GM.active_expression_comTF=active_expression_gm;
%}

function [active_genes, active_expression]=find_active(x,gene_table)
    active_genes=intersect(x,gene_table.Properties.RowNames);
    active_expression=gene_table.NormTPM(active_genes);  
end