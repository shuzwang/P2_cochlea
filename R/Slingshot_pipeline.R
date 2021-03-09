### Slingshot Trajectory Analysis
library(slingshot)
library(uwot)
library(mclust)
library(RColorBrewer)
library(chromVAR)
library(SnapATAC)
library(SingleCellExperiment)
library(tidyr)
library(dplyr)
library(ggplot2)
library(CellTrails)

# Load the chromVAR results
chromvar_results = readRDS("~/Dropbox/P2_project/Manuscript/P2_analysis/P2_scATACseq_chromVAR_dev_hocomoco_v10.RDS")
zscore = deviationScores(chromvar_results)

# Subset the HC single cells
snap_atac = readRDS("~/Dropbox/P2_project/Manuscript/P2_analysis/P2_scATACseq_snap_object_pmat.RDS")
cellbarcode = (snap_atac@metaData$barcode[which(snap_atac@cluster == 1)])

hc_zscore = zscore[,(colnames(zscore) %in% cellbarcode)]

# Load celltrails states
hc_celltrails = readRDS("~/Dropbox/P2_project/Manuscript/P2_analysis/P2_scATACseq_celltrails_hc.RDS")

# Start the Celltrails algrithm
hc_slingshot = SingleCellExperiment(assays = List(counts = hc_zscore,norm = hc_zscore))

# Dimensionality Reduction using PCA
pca = prcomp(t((assays(hc_slingshot)$norm)), scale. = T)
rd1 = pca$x[,1:2]
plot(rd1,col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(hc_slingshot) = SimpleList(PCA = rd1)

# Cluster Identification
cl1 = Mclust(rd1)$classification
colData(hc_slingshot)$GMM = cl1
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

# Run Slingshot
hc_slingshot = slingshot(hc_slingshot, clusterLabels = 'GMM', reducedDim = 'PCA')

colors = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol = colors[cut(hc_slingshot$slingPseudotime_1, breaks=100)]

# Get trajectory
lin1 = getLineages(rd1, cl1, start.clus = c("1"),end.clus = c('2',"3"))
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(lin1, lwd = 3, col = 'black')

crv1 = getCurves(lin1)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], asp = 1, pch = 16)
lines(crv1, lwd = 3, col = 'black')

feature_df = as.data.frame(reducedDims(hc_slingshot)$PCA)
colnames(feature_df)[1:2] = c("PC1","PC2")
feature_df$Experiment = rownames(feature_df)
feature_df$Experiment = gsub(".*-","",feature_df$Experiment)
feature_df$Experiment[which(feature_df$Experiment == 1)] = "Apex"
feature_df$Experiment[which(feature_df$Experiment == 2)] = "Base"
feature_df$cluster1 = cl1
feature_df = cbind(feature_df,t(hc_zscore))

all.equal(rownames(feature_df),colnames(hc_celltrails))
feature_df$celltrails = hc_celltrails$CellTrails.state

## Make formal plot for slingshot results
s1 = as.data.frame(crv1@curves$curve1$s)
s1 = s1[crv1@curves$curve1$ord,]
rownames(s1) = 1:nrow(s1)
s2 = as.data.frame(crv1@curves$curve2$s)
s2 = s2[crv1@curves$curve2$ord,]
rownames(s2) = 1:nrow(s2) 

tiff("~/Dropbox/P2_project/Manuscript/P2_figures/Trajectory_P2_scATACseq_Slingshot_without_legend.tiff",res = 300,width = 8, height = 8, units = 'in')
ggplot(feature_df) + geom_point(aes(PC1,PC2,color = celltrails)) +
  scale_color_manual(values = c("grey80","#377eb8","#4daf4a","#984ea3")) +
  geom_path(data = s1,aes(PC1,PC2),size=1) +
  geom_path(data = s2,aes(PC1,PC2),size=1) +
  theme_classic() +
  theme(axis.title = element_text(size = 30),axis.text = element_text(size = 30))+
  theme(aspect.ratio=1) + 
  theme(axis.line.x  = element_line(size = 0.8),axis.line.y  = element_line(size = 0.8),axis.ticks = element_line(size = 0.8))+
  theme(legend.position = "none",legend.title = element_text(size = 25),legend.text = element_text(size = 25))
dev.off()
