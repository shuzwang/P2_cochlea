## Seurat scRNA-seq cluster analysis pipeline
library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)

# Load scRNA-seq data
# Combine 10x runs
rna_seurat = Read10X(data.dir = "~/Dropbox/P2_project/P2_RNAseq_ApexBase_2856_aggr/outs/filtered_feature_bc_matrix")
rna_seurat = CreateSeuratObject(counts = rna_seurat, project = "P2_scRNAseq")

# Basic QC for scRNA-seq data
rna_seurat[["percent.mt"]] = PercentageFeatureSet(object = rna_seurat, pattern = "mt-")

VlnPlot(object = rna_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(object = rna_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = rna_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

rna_seurat = subset(rna_seurat, subset = nFeature_RNA > 600 & nFeature_RNA < 8000 & percent.mt < 10)

# Normalization
rna_seurat = NormalizeData(object = rna_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
rna_seurat = FindVariableFeatures(object = rna_seurat, selection.method = "vst", nfeatures = 2000)

# Scaling the data
all_genes = rownames(x = rna_seurat)
rna_seurat = ScaleData(object = rna_seurat, features = all_genes)

# Perform linear dimensional reduction
rna_seurat = RunPCA(object = rna_seurat, features = VariableFeatures(object = rna_seurat),verbose =F)

# Examine and visualize PCA results a few different ways
# VizDimLoadings(object = rna_seurat, dims = 1:2, reduction = "pca")
DimPlot(object = rna_seurat, reduction = "pca")

# Add spatial label (apex/base) in meta data
rna_seurat@meta.data$Experiment = rownames(rna_seurat@meta.data)
rna_seurat@meta.data$Experiment = gsub(".*-","",rna_seurat@meta.data$Experiment)
rna_seurat@meta.data$Experiment[which(rna_seurat@meta.data$Experiment == 1)] = "Apex"
rna_seurat@meta.data$Experiment[which(rna_seurat@meta.data$Experiment == 2)] = "Base"

DimPlot(object = rna_seurat, reduction = "pca",group.by = "Experiment")

# Clustering
ElbowPlot(object = rna_seurat,ndims = 50)
rna_seurat = FindNeighbors(object = rna_seurat, dims = 1:10,k.param = 20)
rna_seurat = FindClusters(object = rna_seurat, resolution = 0.5,random.seed = 10)

rna_seurat = RunUMAP(object = rna_seurat, dims = 1:10,seed.use = 20)
DimPlot(object = rna_seurat, reduction = "umap",label = T)
FeaturePlot(rna_seurat,features = "Oc90")
DimPlot(rna_seurat,reduction = "umap",group.by = "Experiment")

FeaturePlot(rna_seurat,features = c("Twist1", "Plp1", "Sox2", "Fgfr3", "Cdh5", "Pou4f3", "Trpm1", "Fgfr2", "Mki67", "Cdh4", "Oc90", "C1qa"))
FeaturePlot(rna_seurat,features = c("Insm1", "Ikzf2","Slc26a5","Otof"))
FeaturePlot(rna_seurat,features = c("Slc26a5","Otof"))

# Cell type annotations
current_cluster_ids = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
new_cluster_ids = c("MES","SAT","PC/DC/IBC","END","SCH","HC","latSC","PRO","medSC","RF","IMM")
names(new_cluster_ids) = levels(new_cluster_ids)
rna_seurat@active.ident = plyr::mapvalues(x = rna_seurat@active.ident, from = current_cluster_ids, to = new_cluster_ids)
DimPlot(object = rna_seurat, reduction = "umap",label = T)

# saveRDS(rna_seurat,"~/Dropbox/P2_project/Manuscript/P2_analysis/P2_scRNAseq_seurat_object.RDS")
# Heatmap for clustering visualization
gene_markers =  c("Twist1","Plp1","Fgfr3","Cdh5","Pou4f3","Trpm1","Fgfr2","Mki67","Cdh4","Oc90","C1qa")
DoHeatmap(rna_seurat,features = gene_markers,raster = F,size = 4.5,angle = 45,draw.lines = F,group.colors = c("#E6242E","#F28429","#F1E911","#97C83E","#65C0A5","#43BBEC","#3952A1","#863F93","#E61B60","#A4A3A3","black"))+
  scale_fill_gradientn(colors = c("white", "yellow","red"))+guides(colour=FALSE)

# Look at the read depth distribution across experiments and clusters to idnetify the doublets
VlnPlot(object = rna_seurat, features = "nCount_RNA", group.by = "Experiment")
VlnPlot(object = rna_seurat, features = "nCount_RNA", group.by = "seurat_clusters")
VlnPlot(object = rna_seurat, features = "nFeature_RNA", group.by = "seurat_clusters")

## DE analysis across cell types
# subset the genes for DE analysis
rna_markers = FindAllMarkers(rna_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "wilcox")
de_rna_markers = subset(rna_markers,rna_markers$p_val_adj < 0.05)
write.table(de_rna_markers,"~/Dropbox/P2_project/Manuscript/P2_analysis/P2_scRNAseq_DE_genes_clusters.txt",sep = "\t",col.names = T,row.names = F,quote = F)
