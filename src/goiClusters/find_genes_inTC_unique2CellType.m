%% 
% Find TC containing unique genes to the cell Type
% Input
% edgeTable with gene expression annotated & comTF expression annotated. 
% Gene Expression Data, genes that are unique to cell type
%
% Output 
% Set of Transcriptional Clusters that are unique to cell type.
clear
clc
close all
%% Load Data
addpath C:\Pore-C\GeneExpression\GenesExpressed
load('Unique_IR_Genes_corrected.mat');
load('UniqueBJ_Genes_corrected.mat');
load('UniqueGM_Genes_corrected.mat');
addpath C:\Pore-C\DATE\01142022
load('BJ_table_w_GE_ComTF_annotated_01142022.mat');
%load('IR_table_w_GE_ComTF_annotated_01142022.mat');
%load('GM_annoted_genes_comTF_01172022.mat');
genes=UniqueBJ;   %% Sub Unique Gene List here

%% RunDATA

%IR=addvars(IR_edgeTable_w_genes_TF,IR_edgeTable_w_genes_TF.genes,'Before','genes','NewVariableNames',{'geneNetwork'});
%IR = IR_edgeTable_w_genes_TF;
%IR.genes=cellfun(@unique2cell,IR.genes,'UniformOutput',false);
%idx1=cellfun(@(x) find_common(x,genes),IR.genes);
%nnz(idx1)



%BJ=addvars(BJ_edgeTable_w_genes_TF,BJ_edgeTable_w_genes_TF.genes,'Before','genes','NewVariableNames',{'geneNetwork'});
BJ.genes=cellfun(@unique2cell,BJ.genes,'UniformOutput',false);
idx2=cellfun(@(x) find_common(x,genes),BJ.genes);
nnz(idx2)

%{
%GM=addvars(GM12878_edgeTable_w_genes_TF,GM12878_edgeTable_w_genes_TF.genes,'Before','genes','NewVariableNames',{'geneNetwork'});
GM.genes=cellfun(@unique2cell,GM.genes,'UniformOutput',false);
idx3=cellfun(@(x) find_common(x,genes),GM.genes);
nnz(idx3)
%}
%setIR=IR(idx1,:);
setBJ=BJ(idx2,:);
%setGM=GM(idx3,:);

%BJ_IR_Clusters = intersect(setIR,setBJ);
%Unique_Potential_Clusters = setBJ; % Change this depedent on cell type
%idx = cellfun(@isempty,Unique_Potential_Clusters.comTF); % to index by whether cell is empty or not
%TClusters = Unique_Potential_Clusters(idx == 0,:); % get rid of the empty cells

function idx=find_common(x,genes)
    idx=any(ismember(x,genes));
end



