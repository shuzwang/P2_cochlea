## Cicero pipeline to predict cis-regulatory interactions
library(cicero)
library(monocle)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Gviz)
library(SnapATAC)
library(GenomicRanges)

## Generate a function to run Cicero 
cicero_connections = function(input_cds,chrom,mouse_mm10_genome,num_dimension = 6,tsne_perplexity = "NULL"){
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(cicero)
  library(monocle)
  
  #Ensure there are no peaks included with zero reads
  input_cds = input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
  
  # Dimention reduction
  set.seed(10)
  atac_cds = detectGenes(input_cds)
  atac_cds = estimateSizeFactors(atac_cds)
  
  if(is.null(tsne_perplexity) == TRUE){
    atac_cds = reduceDimension(atac_cds, max_components = 2, num_dim=num_dimension,reduction_method = 'tSNE', norm_method = "none")
    
  }else{
    atac_cds = reduceDimension(atac_cds, max_components = 2, num_dim=num_dimension,reduction_method = 'tSNE', norm_method = "none",perplexity = 20)
  }
  
  tsne_coords = t(reducedDimA(atac_cds))
  row.names(tsne_coords) = row.names(pData(atac_cds))
  cicero_cds = make_cicero_cds(atac_cds, reduced_coordinates = tsne_coords)
  
  # Co-accessbility
  sample_genome = subset(mouse_mm10_genome, V1 == chrom)
  conns = run_cicero(cicero_cds, sample_genome, sample_num = 2)
  
  # Visualization
  # Get the annotation 
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  mm10_anno_gene_track_regions = GeneRegionTrack(txdb,chromosome = chrom)
  mm10_anno = as.data.frame(mm10_anno_gene_track_regions@range)
  mm10_anno = mm10_anno[,-c(7,8,12)]
  mm10_anno$gene_ensembl = mapIds(org.Mm.eg.db,mm10_anno$gene,column = "ENSEMBL",keytype = "ENTREZID")
  mm10_anno$gene_symbol = mapIds(org.Mm.eg.db,mm10_anno$gene,column = "SYMBOL",keytype = "ENTREZID")
  
  mm10_anno_visual = mm10_anno[,c("seqnames","start","end","strand","feature","gene_ensembl","transcript","gene_symbol")]
  colnames(mm10_anno_visual) = c("chromosome","start","end","strand","feature","gene","transcript","symbol")
  
  return(list(conns,mm10_anno_visual))
}

## Generate function to identify TFBSs around TSS of a particular gene
motif_detection_cicero = function(connection_df,upstream = 50000,downstream = 50000,tss_position,tss_direction = c("positive","negative")){
  library(tidyr)
  library(dplyr)
  
  # Subset the connection data frame based on TSS position
  connect_df = separate(connection_df,col = Peak1,into = c("chr","bp1","bp2"),sep = "_")
  connect_df$bp1 = as.numeric(connect_df$bp1)
  connect_df$bp2 = as.numeric(connect_df$bp2)
  
  if(tss_direction == "negative"){
    connect_df = subset(connect_df,(connect_df$bp2 <= tss_position + upstream) & (connect_df$bp1 >= (tss_position - downstream)))
  }
  
  if(tss_direction == "positive"){
    connect_df = subset(connect_df,(connect_df$bp2 <= tss_position + downstream) & (connect_df$bp1 >= (tss_position - upstream)))
  }
  
  # Find the unique peaks
  final_df = connect_df[,1:3]
  final_df = final_df[!duplicated(final_df),]
  
  return(final_df)
  
}

# Load mm10 chrom size
mouse_mm10_genome = read.delim("~/Library/Mobile Documents/com~apple~CloudDocs/Manuscript/P2_analysis/mm10.chrom.sizes.txt",sep="\t",header = F,stringsAsFactors = F)

######## For HC cluster
# Create HC input data object for Cicero
hc_snap_atac = readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Manuscript/P2_analysis/P2_scATACseq_HC_snap_object_pmat.RDS")
hc_snap_atac_pmat = hc_snap_atac@pmat
hc_peak_df = as.data.frame(hc_snap_atac@peak)
hc_peak_df = hc_peak_df[,1:3]
hc_peak_df[,4] = unite(hc_peak_df,col = "site_name",seqnames:end, sep = "_")
colnames(hc_peak_df) = c("chr","bp1","bp2","site_name")
rownames(hc_peak_df) = hc_peak_df$site_name

colnames(hc_snap_atac_pmat) = hc_peak_df$site_name
rownames(hc_snap_atac_pmat) = hc_snap_atac@barcode
hc_snap_atac_pmat = t(hc_snap_atac_pmat)
hc_snap_atac_pmat@x[hc_snap_atac_pmat@x > 0] = 1

hc_cellinfo = data.frame(cells = colnames(hc_snap_atac_pmat),row.names = colnames(hc_snap_atac_pmat))
hc_input_cds =  suppressWarnings(newCellDataSet(hc_snap_atac_pmat,phenoData = new("AnnotatedDataFrame", data=hc_cellinfo),
                                                featureData = new("AnnotatedDataFrame", data=hc_peak_df)))

# Run pre-defined function
hc_conns = cicero_connections(hc_input_cds,chrom = "chr10",mouse_mm10_genome)

# Define a TSS position for S100b gene
tss_position = 76253853
hc_peaks = motif_detection_cicero(hc_conns[[1]],upstream = 50000,downstream = 50000,tss_position = tss_position,tss_direction = "positive")
# write.table(hc_peaks,"~/Library/Mobile Documents/com~apple~CloudDocs/Manuscript/P2_analysis/P2_cicero/S100b_hc_peaks_up50kb_down50kb.bed",sep = "\t",col.names = F,row.names = F,quote = F) 

# Read S100b putative TFBSs from FIMO results
hc_fimo = read.delim("~/Library/Mobile Documents/com~apple~CloudDocs/Manuscript/P2_analysis/P2_cicero/S100b_hc_peaks_up50kb_down50kb_fimo.txt.gz",sep = "\t",header = T,stringsAsFactors = F)

gfi1_hc_fimo = subset(hc_fimo,hc_fimo$motif_id == "GFI1_MOUSE.H10MO.C")
sox9_hc_fimo = subset(hc_fimo,hc_fimo$motif_id == "SOX9_MOUSE.H10MO.B")

hc_plot = plot_connections(hc_conns[[1]], "chr10", 76200000, 76320000, gene_model = hc_conns[[2]], coaccess_cutoff = 0, connection_width = 2, collapseTranscripts = "longest",
                           include_axis_track = T,return_as_list = T)

plotTracks(hc_plot,from = 76200000,to = 76320000,title.width = 0.3, showTitle = TRUE, chromosome = "chr10", sizes = c(1,0.2,0.2,0.2),
           transcriptAnnotation = "symbol",col.border.title="transparent",background.title = "transparent", lwd.border.title = "transparent",
           col.axis = "black", fontsize.group = 6,fontsize=10)

########## For SC cluster
# Create SC input data object for Cicero
sc_snap_atac = readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Manuscript/P2_analysis/P2_scATACseq_sc_snap_object_pmat.RDS")
sc_snap_atac_pmat = sc_snap_atac@pmat
sc_peak_df = as.data.frame(sc_snap_atac@peak)
sc_peak_df = sc_peak_df[,1:3]
sc_peak_df[,4] = unite(sc_peak_df,col = "site_name",seqnames:end, sep = "_")
colnames(sc_peak_df) = c("chr","bp1","bp2","site_name")
rownames(sc_peak_df) = sc_peak_df$site_name

colnames(sc_snap_atac_pmat) = sc_peak_df$site_name
rownames(sc_snap_atac_pmat) = sc_snap_atac@barcode
sc_snap_atac_pmat = t(sc_snap_atac_pmat)
sc_snap_atac_pmat@x[sc_snap_atac_pmat@x > 0] = 1

sc_cellinfo = data.frame(cells = colnames(sc_snap_atac_pmat),row.names = colnames(sc_snap_atac_pmat))

# make CDS
sc_input_cds =  suppressWarnings(newCellDataSet(sc_snap_atac_pmat,phenoData = new("AnnotatedDataFrame", data=sc_cellinfo),
                                                featureData = new("AnnotatedDataFrame", data=sc_peak_df)))
# saveRDS(sc_input_cds,"~/Library/Mobile Documents/com~apple~CloudDocs/Manuscript/P2_analysis/P2_cicero/sc_input_cds.RDS")
sc_conns = cicero_connections(sc_input_cds,chrom = "chr10",mouse_mm10_genome)

tss_position = 76253853
sc_peaks = motif_detection_cicero(sc_conns[[1]],upstream = 50000,downstream = 50000,tss_position = tss_position,tss_direction = "positive")
write.table(sc_peaks,"~/Library/Mobile Documents/com~apple~CloudDocs/Manuscript/P2_analysis/P2_cicero/S100b_sc_peaks_up50kb_down50kb.bed",sep = "\t",col.names = F,row.names = F,quote = F) 

sc_plot = plot_connections(sc_conns[[1]], "chr10", 76200000, 76320000, gene_model = sc_conns[[2]], coaccess_cutoff = 0, connection_width = 2, collapseTranscripts = "longest",
                           include_axis_track = T,return_as_list = T)

sc_fimo = read.delim("~/Library/Mobile Documents/com~apple~CloudDocs/Manuscript/P2_analysis/P2_cicero/S100b_sc_peaks_up50kb_down50kb_fimo.txt.gz",sep = "\t",header = T,stringsAsFactors = F)
sox9_sc_fimo = subset(sc_fimo,sc_fimo$motif_id == "SOX9_MOUSE.H10MO.B")
sox9_sc_granges = sox9_sc_fimo[,3:6]
colnames(sox9_sc_granges) = c("motif_chr","motif_start","motif_end","motif_strand")
sox9_sc_granges = makeGRangesFromDataFrame(sox9_sc_granges,keep.extra.columns = T,ignore.strand = T,seqnames.field = "motif_chr",start.field = "motif_start",end.field = "motif_end",strand.field = "motif_strand")

sc_plot[[5]] = Gviz::AnnotationTrack(sox9_sc_granges, name="Sox9_TFBSs", fill = "#fb9a99",lwd = .0000001, fontsize.group=6,
                                     fontsize=6, cex.feature = 0.5)

gfi1_sc_fimo = subset(sc_fimo,sc_fimo$motif_id == "GFI1_MOUSE.H10MO.C")
gfi1_sc_granges = gfi1_sc_fimo[,3:6]
colnames(gfi1_sc_granges) = c("motif_chr","motif_start","motif_end","motif_strand")
gfi1_sc_granges = makeGRangesFromDataFrame(gfi1_sc_granges,keep.extra.columns = T,ignore.strand = T,seqnames.field = "motif_chr",start.field = "motif_start",end.field = "motif_end",strand.field = "motif_strand")

sc_plot[[6]] = Gviz::AnnotationTrack(gfi1_sc_granges, name="Gfi1_TFBSs", fill = "#fdbf6f",lwd = .0000001, fontsize.group=6,
                                     fontsize=6, cex.feature = 0.5)
