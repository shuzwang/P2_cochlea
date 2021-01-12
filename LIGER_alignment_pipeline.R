## Liger Alignment
library(liger)
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)

# Load scRNA-seq data
rna_seurat = readRDS("~/Dropbox/P2_project/Manuscript/P2_analysis/P2_scRNAseq_seurat_object.RDS")
rna_data = rna_seurat@assays$RNA@counts

# Load scATAC-seq data
atac_snapatac = readRDS("~/Dropbox/P2_project/Manuscript/P2_analysis/P2_scATACseq_snap_object_pmat.RDS")
atac_data = atac_snapatac@gmat

# Manually define variable genes as scRNA-seq cluster specific DEGs
de_genes_df = read.delim("~/Dropbox/P2_project/Manuscript/P2_analysis/P2_scRNAseq_DE_genes_clusters.txt",sep = "\t",header = T,stringsAsFactors = F)
de_genes_list = unique(de_genes_df$gene)
de_genes_list = de_genes_list[which(de_genes_list %in% liger_object@raw.data$rna@Dimnames[[1]])]
de_genes_list = de_genes_list[which(de_genes_list %in% liger_object@raw.data$atac@Dimnames[[1]])]

# Create liger object
liger_object = createLiger(list(rna = rna_data, atac = t(atac_data)))
liger_object = liger::normalize(liger_object)
#liger_object = selectGenes(liger_object, datasets.use = 1)
liger_object@var.genes = de_genes_list
liger_object = scaleNotCenter(liger_object)

# Add given cell cluster info
liger_object@cell.data$known_ID = NA
liger_object@cell.data$known_ID[1:695] = rna_seurat@meta.data$RNA_snn_res.0.5
liger_object@cell.data$known_ID[696:nrow(liger_object@cell.data)] = atac_snapatac@cluster

liger_object@cell.data$known_cluster = NA
liger_object@cell.data$known_cluster[1:695] = as.character(rna_seurat@active.ident)
atac_cluster = factor(atac_snapatac@cluster,levels = 1:6,labels = c("HC","PC/DC/IBC","RF","MES","END","IMM"))
liger_object@cell.data$known_cluster[696:nrow(liger_object@cell.data)] = as.character(atac_cluster)

# Dimension reduction in Liger
liger_object = optimizeALS(liger_object, k = 20,rand.seed = 2)
liger_object = quantile_norm(liger_object)
liger_object = louvainCluster(liger_object, resolution = 0.2)

liger_object = runUMAP(liger_object, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
plotByDatasetAndCluster(liger_object,axis.labels = c('UMAP 1', 'UMAP 2'),pt.size=1)

anno1 = liger_object@cell.data$known_cluster
names(anno1) = rownames(liger_object@cell.data)
plotByDatasetAndCluster(liger_object,cluster = anno1,axis.labels = c('UMAP 1', 'UMAP 2'),pt.size=1)