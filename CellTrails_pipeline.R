### CellTrails pipeline to analyze scATAC-seq trajectory
library(CellTrails)
library(chromVAR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Seurat)

# Load the chromVAR results
chromvar_results = readRDS("~/Desktop/P2_analysis/P2_scATACseq_chromVAR_dev_hocomoco_v10.RDS")
zscore = deviationScores(chromvar_results)

# Subset the HC single cells
snap_atac = readRDS("~/Desktop/P2_analysis/P2_scATACseq_snap_object_pmat.RDS")
cellbarcode = (snap_atac@metaData$barcode[which(snap_atac@cluster == 1)])

hc_zscore = zscore[,(colnames(zscore) %in% cellbarcode)]

# Start the Celltrails algrithm
hc_celltrails = SingleCellExperiment(assays = list(counts = hc_zscore,logcounts = hc_zscore))

se = embedSamples(hc_celltrails)
d = findSpectrum(se$eigenvalues, frac=100)
latentSpace(hc_celltrails) = se$components[, d]

# Find states
cl = findStates(hc_celltrails, min_size=0.1, min_feat=5, max_pval=1e-4, min_fc=2)
states(hc_celltrails) = cl
plotManifold(hc_celltrails, color_by="phenoName", name="state")

# Add spatial label in phenoNames
colData(hc_celltrails)$Experiment = rownames(colData(hc_celltrails))
colData(hc_celltrails)$Experiment = gsub(".*-","",colData(hc_celltrails)$Experiment)
colData(hc_celltrails)$Experiment[which(colData(hc_celltrails)$Experiment == 1)] = "Apex"
colData(hc_celltrails)$Experiment[which(colData(hc_celltrails)$Experiment == 2)] = "Base"

plotManifold(hc_celltrails, color_by="phenoName", name="Experiment")

# Sample ordering
hc_celltrails = connectStates(hc_celltrails, l=4)
plotStateTrajectory(hc_celltrails, color_by="phenoName", name="Experiment", component=1, point_size=1.5, label_offset=4)
plotStateTrajectory(hc_celltrails, color_by="featureName", name="ATOH1_MOUSE.H10MO.C", component=1, point_size=5)
plotStateTrajectory(hc_celltrails, color_by="featureName", name="INSM1_MOUSE.H10MO.C", component=1, point_size=5)

# Make trajectory
hc_celltrails = selectTrajectory(hc_celltrails, component=1)
hc_celltrails = fitTrajectory(hc_celltrails)

write.ygraphml(sce=hc_celltrails,file='~/Desktop/P2_analysis/P2_scATACseq_chromVAR_celltrails_hc.graphml',color_by='phenoName',name='state',node_label='state')

tl = read.ygraphml("~/Desktop/P2_analysis/P2_scATACseq_chromVAR_celltrails_hc.graphml")
plot(tl[,1:2], axes=FALSE, xlab="", ylab="", pch=20, cex=.25)
trajLayout(hc_celltrails, adjust=TRUE) = tl

# saveRDS(hc_celltrails,"~/Desktop/Project_scATACseq/P2_analysis/P2_scATACseq_celltrails_hc.RDS")
