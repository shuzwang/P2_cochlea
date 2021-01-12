## SnapATAC pipepline for scATAC-seq analysis
library(SnapATAC)
library(tidyr)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(reticulate)
library(leiden)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(reshape2)

# Load files
snap_files = "~/Desktop/P2_analysis/atac_peak_cell_20000.snap"
sample_names = "scATACseq_5K"
barcode_files = "~/Desktop/P2_ATAC_ApexBase_aggr_2872/outs/singlecell.csv"

# Create snap object
snap_atac = createSnap(file=snap_files,sample=sample_names)

########## QC processes
barcode_ls = lapply(seq(snap_files), function(i){
  barcodes = read.csv(
    barcode_files[i], 
    head=TRUE
  );
  # remove NO BAROCDE line
  barcodes = barcodes[2:nrow(barcodes),];
  barcodes$logUMI = log10(barcodes$passed_filters + 1);
  barcodes$promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
  barcodes
})
barcode_ls = as.data.frame(barcode_ls[[1]])
ggplot(barcode_ls, aes(x=logUMI, y=promoter_ratio)) + geom_point(size=0.3, col="grey") +theme_classic()	+ggtitle(sample_names) +ylim(0, 1) + xlim(0, 6) + labs(x = "log10(UMI)", y="promoter ratio")

barcode_ls = barcode_ls[which(barcode_ls$logUMI >= 3 & barcode_ls$logUMI <= 5 & barcode_ls$promoter_ratio >= 0.2 & barcode_ls$promoter_ratio <= 0.8),]

barcode_shared = intersect(snap_atac@barcode, barcode_ls$barcode)
snap_atac = snap_atac[match(barcode_shared, snap_atac@barcode),]
barcode_ls = barcode_ls[match(barcode_shared, barcode_ls$barcode),]
snap_atac@metaData = barcode_ls

# Binarization & add binary cell by bin matrix to the snap object
snap_atac = addBmatToSnap(snap_atac, bin.size=5000)
bin_coverage = snap_atac@bmat
snap_atac = makeBinary(snap_atac, mat="bmat")

# Counts for each bin
bin_counts = apply(snap_atac@bmat, 2, sum)
hist(bin_counts)

# Remove blacklist
black_list = read.table("~/Desktop/P2_analysis/mm10.blacklist.bed.gz");
black_list_gr = GRanges(black_list[,1], IRanges(black_list[,2], black_list[,3]))
idy = queryHits(findOverlaps(snap_atac@feature, black_list_gr))
if(length(idy) > 0){snap_atac = snap_atac[,-idy, mat="bmat"]}

# Remove unwant chromosomes
chr_exclude = seqlevels(snap_atac@feature)[grep("random|chrM", seqlevels(snap_atac@feature))]
idy = grep(paste(chr_exclude, collapse="|"), snap_atac@feature)
if(length(idy) > 0){snap_atac = snap_atac[,-idy, mat="bmat"]}

########## Dimension reduction
# Dimension reduction & clustering
set.seed(20)
snap_atac = runDiffusionMaps(obj=snap_atac,input.mat="bmat", num.eigs=30)
plotDimReductPW(obj=snap_atac, eigs.dims=1:30,point.size=0.3,point.color="grey", point.shape=19, point.alpha=0.6,down.sample=5000,pdf.file.name=NULL, pdf.height=7, pdf.width=7)

set.seed(10)
snap_atac = runKNN(obj=snap_atac,eigs.dims=1:15,k=15)
snap_atac=runCluster(obj=snap_atac,tmp.folder=tempdir(),louvain.lib="leiden",seed.use=10,resolution=0.2)
snap_atac = runViz(obj=snap_atac, tmp.folder=tempdir(),dims=2,eigs.dims=1:15, method="umap",seed.use=10)

########## Gene-based annotation
txdb_mm10 = TxDb.Mmusculus.UCSC.mm10.knownGene
txdb_gene = as.data.frame(genes(txdb_mm10))

anno_gene = select(org.Mm.eg.db,keys = txdb_gene$gene_id,columns = c("SYMBOL","ENTREZID"),keytype = "ENTREZID")
colnames(anno_gene) = c("gene_id","name")
txdb_gene = left_join(txdb_gene,anno_gene,by = "gene_id")
txdb_gene_gr = makeGRangesFromDataFrame(txdb_gene,keep.extra.columns=T)

# Create gmat (cell by gene matrix)
snap_atac = createGmatFromMat(obj = snap_atac,input.mat = "bmat",genes = txdb_gene_gr,do.par = T,num.cores = 2)
# saveRDS(snap_atac,"~/Desktop/P2_analysis/P2_scATACseq_snap_object_gmat.RDS")

########## Peak calling for each cluster
peaks = runMACS(obj=snap_atac[which(snap_atac@cluster=="1"),], output.prefix="P2_scATACseq_cluster1", path.to.snaptools="/Users/shuzewang/miniconda3/bin/snaptools",path.to.macs="/Users/shuzewang/miniconda3/bin/macs2",gsize="mm",buffer.size=500, num.cores=2,macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",tmp.folder=tempdir())
peaks = runMACS(obj=snap_atac[which(snap_atac@cluster=="2"),], output.prefix="P2_scATACseq_cluster2", path.to.snaptools="/Users/shuzewang/miniconda3/bin/snaptools",path.to.macs="/Users/shuzewang/miniconda3/bin/macs2",gsize="mm",buffer.size=500, num.cores=2,macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",tmp.folder=tempdir())
peaks = runMACS(obj=snap_atac[which(snap_atac@cluster=="3"),], output.prefix="P2_scATACseq_cluster3", path.to.snaptools="/Users/shuzewang/miniconda3/bin/snaptools",path.to.macs="/Users/shuzewang/miniconda3/bin/macs2",gsize="mm",buffer.size=500, num.cores=2,macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",tmp.folder=tempdir())
peaks = runMACS(obj=snap_atac[which(snap_atac@cluster=="4"),], output.prefix="P2_scATACseq_cluster4", path.to.snaptools="/Users/shuzewang/miniconda3/bin/snaptools",path.to.macs="/Users/shuzewang/miniconda3/bin/macs2",gsize="mm",buffer.size=500, num.cores=2,macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",tmp.folder=tempdir())
peaks = runMACS(obj=snap_atac[which(snap_atac@cluster=="5"),], output.prefix="P2_scATACseq_cluster5", path.to.snaptools="/Users/shuzewang/miniconda3/bin/snaptools",path.to.macs="/Users/shuzewang/miniconda3/bin/macs2",gsize="mm",buffer.size=500, num.cores=2,macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",tmp.folder=tempdir())
peaks = runMACS(obj=snap_atac[which(snap_atac@cluster=="6"),], output.prefix="P2_scATACseq_cluster6", path.to.snaptools="/Users/shuzewang/miniconda3/bin/snaptools",path.to.macs="/Users/shuzewang/miniconda3/bin/macs2",gsize="mm",buffer.size=500, num.cores=2,macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",tmp.folder=tempdir())

# Find DARs
peaks_names = c("P2_scATACseq_cluster1_peaks.narrowPeak","P2_scATACseq_cluster2_peaks.narrowPeak","P2_scATACseq_cluster3_peaks.narrowPeak","P2_scATACseq_cluster4_peaks.narrowPeak","P2_scATACseq_cluster5_peaks.narrowPeak","P2_scATACseq_cluster6_peaks.narrowPeak")
folder_path = "~/Desktop/P2_analysis/P2_scATACseq_cluster_peaks/"
peaks_names = paste(folder_path,peaks_names,sep = "")

peak_gr_ls = lapply(peaks_names, function(x){
  peak_df = read.table(x)
  GRanges(peak_df[,1], IRanges(peak_df[,2], peak_df[,3]))
})
peak_df = reduce(Reduce(c, peak_gr_ls))
peak_df = as.data.frame(peak_df)[,1:3]
# write.table(peak_df,file = "~/Desktop/P2_analysis/P2_scATACseq_overall_peaks.bed",append=FALSE,quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),fileEncoding = "")

########## Add cell-by-peak matrix
snap_atac = createPmat(snap_atac, peaks=peak_df,do.par=TRUE,num.cores=2)
# saveRDS(snap_atac,"~/Desktop/P2_analysis/P2_scATACseq_snap_object_pmat.RDS")

