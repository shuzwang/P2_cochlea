### chromVAR pipeline for scATAC-seq analysis
library(TFBSTools)
library(motifmatchr)
library(chromVAR)
library(SnapATAC)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Matrix)
library(BSgenome.Mmusculus.UCSC.mm10)
library(SummarizedExperiment)
library(chromVARmotifs)
library(BiocParallel)

# Prepare the dataset to run chromVAR
snap_atac = readRDS("~/Desktop/P2_analysis/P2_scATACseq_snap_object_pmat.RDS")
atac_matrix = as.matrix(snap_atac@pmat)
atac_matrix = t(atac_matrix)
peaks = snap_atac@peak

# Load HOCOMOCOv10 database to prepare for matchMotifs() function
motif_pwm_list = readRDS("~/Desktop/P2_analysis/hocomoco_v10_mouse_motifs.RDS")

# Generate chromVAR object
fragment_counts = SummarizedExperiment(assays = list(counts = atac_matrix),rowRanges = peaks)

# Add GC content
fragment_counts = addGCBias(fragment_counts,genome = BSgenome.Mmusculus.UCSC.mm10)

# Filter peaks
counts_filtered = filterPeaks(fragment_counts, non_overlapping = TRUE,min_fragments_per_peak = 3)

# Annotation
motif_ix = matchMotifs(motif_pwm_list, counts_filtered,genome = BSgenome.Mmusculus.UCSC.mm10)

# Compute deviations
set.seed(10)
register(BiocParallel::SerialParam())

dev = computeDeviations(object = counts_filtered,annotations = motif_ix)
# saveRDS(dev,"~/Desktop/P2_analysis/P2_scATACseq_chromVAR_dev.RDS")

# TF-by-cell matrix
zscore_matrix = deviationScores(dev)

## Cluster from tsne in chromVAR
set.seed(10)

tsne_results = deviationsTsne(dev, threshold = 1.5, perplexity = 30, shiny = FALSE)
dev@colData$CellType = feature_plot_df$cluster

tsne_plots = plotDeviationsTsne(dev, tsne_results,sample_column = "CellType", shiny = FALSE)
tsne_plots[[1]]
