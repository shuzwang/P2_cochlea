#### Jaccard Index Similarity Matrix Functional Code
### This code is to determine the cell identities for each scATAC-seq cluster by calculating the Jaccard index for each pair
### of clusters for scATAC-seq and scRNA-seq. 

### Identify differentially expressed TF genes between 2 populations using scRNA-seq data
#' @param rna_genes A data frame of DE genes of from scRNA-seq. For instance, DE results from FindAllMarkers function in Seurat.
#' @param atac_genes A data frame of annotated DARs (annotated genes) from scATAC-seq. For instance, the output of annotate_atac_peaks
#' function could be used here.
#' @return A Jaccard index similarity matrix with scRNA-seq clusters in the rows and scATAC-seq clusters in the columns. 

# Calculate Jaccard Index Matrix between DE genes from scRNA-seq and annotated genes from scATAC-seq
calculate_jaacard_index_matrix = function(rna_genes,atac_genes){
  
  # Separate the data frame based on clusters
  rna_cluster_list = list()
  rna_cluster_length = list()
  for (i in 1:length(unique(rna_genes$cluster))) {
    rna_cluster_list[[i]] = subset(rna_genes,rna_genes$cluster == unique(rna_genes$cluster)[i])
    rna_cluster_length[i] = nrow(rna_cluster_list[[i]])
  }
  rna_cluster_length = unlist(rna_cluster_length)
  
  atac_cluster_list = list()
  atac_cluster_length = list()
  for(i in 1:length(unique(atac_genes$cluster))){
    atac_cluster_list[[i]] = subset(atac_genes,atac_genes$cluster == unique(atac_genes$cluster)[i])
    atac_cluster_length[i] = length(unique(atac_cluster_list[[i]]$SYMBOL))
  }
  atac_cluster_length = unlist(atac_cluster_length)
  
  # Calculate Jaccard Index Matrix 
  jaccard_index_matrix = matrix(data = NA,nrow = length(rna_cluster_length),ncol = length(atac_cluster_length))
  for(i in 1:length(rna_cluster_length)){
    for(j in 1:length(atac_cluster_length)){
      intersected_genes = intersect(rna_cluster_list[[i]]$gene,unique(atac_cluster_list[[j]]$SYMBOL))
      jaccard_index_matrix[i,j] = length(intersected_genes)/(rna_cluster_length[i]+atac_cluster_length[j]-length(intersected_genes))
    }
  }
  jaccard_index_matrix = as.data.frame(jaccard_index_matrix)
  colnames(jaccard_index_matrix) = unique(atac_genes$cluster)
  rownames(jaccard_index_matrix) = unique(rna_genes$cluster)
  
  return(jaccard_index_matrix)
}


# Annotate peaks to genomic regions using proximity-based method
#' @param peak_file Standard peak file in .narrowPeak format from MACS2. 
#' @param cluster A meta information user wants to add in the file. For example, "cluster" = "HC" represents the cell type identities 
#' for the peak file.
#' @param txdb Genome reference annotation database generated from UCSC wrapped up in TxDb object. By default, txdb is 
#' TxDb.Mmusculus.UCSC.mm10.knownGene. 
#' @param annodb Genomewide annotation for mouse based on mapping using gene information. By default, annodb is org.Mm.eg.db. 
#' @return An annotated data frame using nearest gene method with meta data. 

annotate_atac_peaks = function(peak_file,cluster,txdb = TxDb.Mmusculus.UCSC.mm10.knownGene,annodb = 'org.Mm.eg.db'){
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  
  ## Prepare the peak file
  peak_df = peak_file[,1:5]
  
  ## Keep the meta data for future use
  if(ncol(peak_file) > 5){
    meta_df = peak_file[,6:ncol(peak_file)]
  }
  
  colnames(peak_df) = c("chrom","start","end","width","strand")
  regions = makeGRangesFromDataFrame(peak_df,keep.extra.columns = T,seqnames.field = "chrom",start.field = "start",end.field = "end",strand.field = "strand")
  
  # Add cluster information to the peaks
  elementMetadata(regions)[["cluster"]] = cluster
  
  # Annotate the peaks
  regionAnno = as.data.frame(annotatePeak(regions, TxDb = txdb, annoDb = annodb))
  regionAnno[grep("Exon", regionAnno$annotation), "annotation"] = "Exon"
  regionAnno[grep("Intron", regionAnno$annotation), "annotation"] = "Intron"
  
  # Add the meta data for future use
  if(ncol(peak_file) > 5){
    regionAnno = cbind(regionAnno,meta_df)
  }
  
  return(regionAnno)
}
