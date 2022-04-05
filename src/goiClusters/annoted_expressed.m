%% 
% To find Genes that are Expressed in the Transcription Cluster
% Input
% edgeTable from porecAnalysisTransFac3.m script
% Gene Expression Data, set threshold >0 or >2 
%
% Output 
% edgeTable with gene expression annotated.
%% typical startup stuff
disp('annotated expression')
clear
clc
close all

%% load data

% BJ Single Cell Data
addpath C:\Pore-C\BJ_SingleCell
load('scBJ_edgeTable_w_genes_TF_noenhancer_03252022.mat')

addpath C:\Pore-C\GeneExpression\GenesExpressed
load("BJ_GenesTPMgreater0.mat");

%{
%% process IR
IR=addvars(edgeTable_w_genes_TF,edgeTable_w_genes_TF.genes,'Before','genes','NewVariableNames',{'geneNetwork'});
idx_IR = ~cellfun(@isempty,IR.comTF); % to index by whether cell is empty or not
IR=IR(idx_IR,:);
IR.genes=cellfun(@unique2cell,IR.genes,'UniformOutput',false);
[active_genes_ir,active_expression_ir]=cellfun(@(x) find_active(x,IR_Genes_Expressed), IR.genes,'UniformOutput',false);
IR.active_genes=active_genes_ir;
IR.active_expression=active_expression_ir;
%}

%% process BJ
BJ=addvars(edgeTable_w_genes_TF,edgeTable_w_genes_TF.genes,'Before','genes','NewVariableNames',{'geneNetwork'});
idx_BJ = ~cellfun(@isempty,BJ.comTF); % to index by whether cell is empty or not
BJ=BJ(idx_BJ,:);
BJ.genes=cellfun(@unique2cell,BJ.genes,'UniformOutput',false);
[active_genes_bj,active_expression_bj]=cellfun(@(x) find_active(x,BJ_Genes_Expressed), BJ.genes,'UniformOutput',false);
BJ.active_genes=active_genes_bj;
BJ.active_expression=active_expression_bj;

%{
%% process GM12878
GM_1=addvars(GM_edgeTable_w_genes,GM_edgeTable_w_genes.genes,'Before','genes','NewVariableNames',{'geneNetwork'});
%GM = edgeTable_w_genes_TF;
idx_GM = ~cellfun(@isempty,GM_1.comTF); % to index by whether cell is empty or not
GM_1=GM_1(idx_GM,:);
GM_1.genes=cellfun(@unique2cell,GM_1.genes,'UniformOutput',false);
[active_genes_gm,active_expression_gm]=cellfun(@(x) find_active(x,GM12878_Genes_Expresssed), GM_1.genes,'UniformOutput',false);
GM_1.active_genes=active_genes_gm;
GM_1.active_expression=active_expression_gm;
%}
%% core function
function [active_genes, active_expression]=find_active(x,gene_table)
    active_genes=intersect(x,gene_table.Properties.RowNames);
    active_expression=gene_table.NormTPM(active_genes);  
end
