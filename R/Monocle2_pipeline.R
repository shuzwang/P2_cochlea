## Monocle 2 trajectory pipeline
library(monocle)
library(chromVAR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Matrix)

# Load the TF-by-cell matrix from chromVAR
chromvar_results = readRDS("~/Dropbox/P2_project/Manuscript/P2_analysis/P2_scATACseq_chromVAR_dev_hocomoco_v10.RDS")
zscore = deviationScores(chromvar_results)

# Subset the HC single cells
snap_atac = readRDS("~/Dropbox/P2_project/Manuscript/P2_analysis/P2_scATACseq_snap_object_pmat.RDS")
cellbarcode = (snap_atac@metaData$barcode[which(snap_atac@cluster == 1)])
hc_zscore = zscore[,(colnames(zscore) %in% cellbarcode)]

# Start the Monocle algrithm
# format cell info
cellinfo = data.frame(cells = colnames(hc_zscore),row.names = colnames(hc_zscore))
# format TF info
geneinfo = data.frame(genes = rownames(hc_zscore),row.names = rownames(hc_zscore))

# Preparation
input_cds = suppressWarnings(newCellDataSet(hc_zscore,phenoData = new("AnnotatedDataFrame", data=cellinfo),featureData = new("AnnotatedDataFrame", data=geneinfo),expressionFamily = gaussianff()))

# Trajectory analysis
input_cds = reduceDimension(input_cds, max_components = 4,reduction_method = 'DDRTree',norm_method = "none",scaling = T)
input_cds = orderCells(input_cds)
plot_cell_trajectory(input_cds, color_by = "State")
plot_cell_trajectory(input_cds, color_by = "Pseudotime")

# Project CellTrails states onto the Monocle trajectory
hc_celltrails = readRDS("~/Dropbox/P2_project/Manuscript/P2_analysis/P2_scATACseq_celltrails_hc.RDS")
input_cds$celltrails = hc_celltrails$CellTrails.state

# Visualization
tiff("~/Dropbox/P2_project/Manuscript/P2_figures/Trajectory_P2_scATACseq_Monocle2_without_legend.tiff",res = 300,width = 8, height = 8, units = 'in')
plot_cell_trajectory(input_cds, 1, 2, color = "celltrails",cell_link_size=1,show_branch_points=F) + 
  scale_color_manual(values = c("grey80","#377eb8","#4daf4a","#984ea3")) +
  theme_classic() +
  theme(axis.title = element_text(size = 30),axis.text = element_text(size = 30))+
  theme(aspect.ratio=1) + 
  theme(axis.line.x  = element_line(size = 0.8),axis.line.y  = element_line(size = 0.8),axis.ticks = element_line(size = 0.8))+
  theme(legend.position = "none",legend.title = element_text(size = 25),legend.text = element_text(size = 25))
dev.off()
