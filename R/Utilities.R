#### Utilities
### This code is utility code. 

### Calculate the number of overlaps and show the overlapped genes
#' @param rna_de_genes A data frame of DE genes of from scRNA-seq. For instance, DE results from FindAllMarkers function in Seurat.
#' @param atac_da_genes A data frame of annotated DARs (annotated genes) from scATAC-seq. For instance, the output of annotate_atac_peaks
#' function could be used here.
#' @return A list of overlapped genes results. The first element is a data frame containing the number of overlapped genes between 
#' scRNA-seq cluster and scATAC-seq cluster. The scRNA-seq clusters are in columns and scATAC-seq clusters are in rows. The second 
#' element is a list of overlapped genes. 

overlapped_genes = function(rna_de_genes,atac_da_genes){
  library(stringr)
  
  rna_list = list()
  atac_list = list()
  
  # Separate the data based on clusters
  rna_cluster = unique(rna_de_genes$cluster)
  for (i in 1:length(unique(rna_de_genes$cluster))) {
    rna_list[[i]] = subset(rna_de_genes,rna_de_genes$cluster == rna_cluster[i])
  }
  
  for(j in 1:length(unique(atac_da_genes$cluster))){
    atac_list[[j]] = subset(atac_da_genes,atac_da_genes$cluster == j)
  }
  
  # Calculate the number of overlapped genes
  overlapped_matrix = data.frame()
  overlapped_gene_list = data.frame()
  
  for(i in 1:length(rna_list)){
    rna_gene_list = unique(as.data.frame(rna_list[[i]])$gene)
    for(j in 1:length(atac_list)){
      atac_gene_list = unique(as.data.frame(atac_list[[j]])$SYMBOL)
      overlapped_matrix[j,i] = length(intersect(rna_gene_list,atac_gene_list))
      overlapped_gene_list[j,i] = str_c(intersect(rna_gene_list,atac_gene_list),collapse = ",")
    }
  }
  colnames(overlapped_matrix) = rna_cluster
  colnames(overlapped_gene_list) = rna_cluster
  
  return(list(overlapped_matrix,overlapped_gene_list))
}


### Visualize the QC results from scATAC-seq with smoothing
#' @param df A data frame with QC scores (e.g. TSS enrichment scores, fragment length proportion, etc.) with each row represents 
#' for a particular position. 
#' @param window A size of window to aggregate for taking the average to smooth the data points. By default, window is 10bp.
#' @param spline_method Splines are a smooth and flexible way of fitting nonlinear models and learning the nonlinear 
#' interactions from the data. spline_method is the method users want to choose smoothing the data. By default, spline_method
#' is "natural". There are other spline methods inherited from spline function in stats package. For example, "fmm", "periodic",
#' "monoH.FC", and "hyman" are provided. 
#' @param df_length An interval of data points along the x-axis. Here the df_length is the length of genomic regions (e.g. 1kb 
#' around the TSS). By default, df_length is 1001bp. 
#' @return A data frame of smoothed data after aggregation. 

smooth_fragment_length = function(df,window = 10,spline_method = "natural",df_length = 1001){
  
  # Aggregate windows and take the average
  return_df = data.frame(matrix(ncol = 7, nrow = df_length))
  return_df[,1] = df[,1]
  for(j in 2:ncol(df)){
    for(i in 1:nrow(df)){
      if(i != 1){
        sum_read = sum_read + df[i,j]
      }else{
        sum_read = df[i,j]
      }
      
      if((i-1) %% window == 0){
        if((i-1) == 0){
          return_df[i,j] = df[i,j]
        }else{
          return_df[i,j] = sum_read/window
          sum_read = 0
        }
      }
    }
  }
  return_df = return_df[complete.cases(return_df),]
  
  # Spline for smoothing
  smooth_df = list()
  for (i in 2:ncol(return_df)) {
    smooth_df[[i-1]] = as.data.frame(spline(return_df[,1], return_df[,i],n = df_length,method = spline_method))
    smooth_df[[i-1]]$y[which(smooth_df[[i-1]]$y < 0)] = 0
  }
  
  smooth_return = cbind(smooth_df[[1]],smooth_df[[2]]$y,smooth_df[[3]]$y,smooth_df[[4]]$y,smooth_df[[5]]$y,smooth_df[[6]]$y)
  colnames(smooth_return) = c("length","cluster1","cluster2","cluster3","cluster4","cluster5","cluster6")
  
  return(smooth_return)
}


### Make a single curtain plot using scATAC-seq data
#' The function is 
#' @param atac_snap An object of SnapATAC with pmat (peak information). 
#' @param rotated_data_frame A data frame with rank for each individual cell. For instance, the output of rotate_data function. 
#' @param atac_hc_annotated_peak A data frame of annotated peaks. It is to select the corresponding genes of DARs or peaks. 
#' @param feature A feature (e.g DAR, peak) name to show in the single curtain plot. 
#' @param seed Set seed for geom_jitter. By default, seed is 10.
#' @param is_regression A boolean parameter to show the fitting line next to the single curtain plot. By default, is_regression 
#' is FALSE. 
#' @param is_experiment A boolean parameter to show the meta-data, Library ID, in single curtain plot. If FALSE, the value of the
#' feature will be shown in single curtain plot. By default, experiment is TRUE. 
#' @return A single scATAC-seq curtain plot. 

atac_single_curtain_plot = function(atac_snap,rotated_data_frame,atac_hc_annotated_peak,is_experiment = TRUE,feature,seed = 10,is_regression = FALSE){
  library(gridExtra)
  library(ggrepel)
  library(ggplot2)
  
  set.seed(seed)
  
  if(is_experiment == TRUE){
    g1 = ggplot(rotated_atac_hc,aes(0,Rank,color = Experiment)) + 
      geom_jitter(size = 3) + 
      coord_fixed(ratio = 0.01) +
      scale_color_manual(values = c("red","blue"))
    theme_bw() + 
      theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())
  }
  else{
    atac_peak_library = annotate_atac_peaks(as.data.frame(atac_snap@peak),cluster = "All")
    atac_peak_library$peak_name = paste(atac_peak_library$seqnames,atac_peak_library$start,atac_peak_library$end,sep = "_")
    
    # DARs
    dar = subset(atac_hc_annotated_peak,atac_hc_annotated_peak$SYMBOL == feature)
    dar_index = paste(dar$seqnames,dar$start,dar$end,sep = "_")
    index = rownames(atac_peak_library)[atac_peak_library$peak_name == dar_index]
    index = as.numeric(index)
    
    peak_index = as.data.frame(atac_snap@pmat[,index])
    colnames(peak_index) = feature
    peak_gradient_df = cbind(rotated_atac_hc,peak_index)
    
    peak_gradient_df$label = feature
    
    # Visualization
    set.seed(seed)
    
    peak_gradient_df[which(peak_gradient_df[,feature] == 0),feature] = NA
    quarter_quantile = quantile(peak_gradient_df[,feature],probs = 0.75,na.rm = T)
    
    g1 = ggplot(peak_gradient_df,aes(0,Rank,color = peak_gradient_df[,feature])) + 
      geom_jitter(size = 3) + coord_fixed(ratio = 0.01) + 
      scale_color_gradient2(low = "#ffff00",high = "#371865",mid = "#1dd41d",midpoint = quarter_quantile,na.value = "grey80")+
      theme_bw() + 
      theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),legend.title = element_text(size = 20),legend.text = element_text(size = 20),legend.position = "none") +
      facet_grid(~label)+
      theme(strip.text.x = element_text(size = 20))+
      labs(color = feature)+
      theme(plot.margin=unit(c(1,-4,1,1), "cm")) 
    
    if(is_regression == FALSE){
      return(g1)
    }
    else{
      g2 = ggplot(data = peak_gradient_df,aes(x = Rank,y = peak_gradient_df[,feature])) + 
        stat_smooth(method = "auto",se = F,color="black") + 
        coord_flip() +
        theme_bw() + 
        theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.x = element_text(size = 20),axis.text.x = element_text(size = 12))+facet_grid(~label) +
        theme(strip.text.x = element_text(size = 20))+
        theme(axis.text.x = element_blank(),axis.ticks = element_blank(),axis.title = element_blank()) +
        theme(aspect.ratio = 5)+
        theme(plot.margin=unit(c(1,1,1,-4), "cm"))
      
      tiff(paste("~/Desktop/spatial_plot_regression_line_",feature,".tiff",sep = ""),res = 300,width = 8, height = 8, units = 'in')
      grid.arrange(g1,g2,ncol = 2, nrow = 1)
      dev.off()
    }
    
  }
  return(g1)
}


## Make dotplot for DE genes in scRNA-seq
#' The function is to show the DE genes in a dotplot. The color of the dots represent the gene expression level. The size of the dots
#' represents percentage of cells expressing a given transcript for the clusters. 
#' @param rna_seurat A Seurat object of scRNA-seq data after normalization.
#' @param gene_list A list of genes users want to show in the dotplot. 
#' @return A data frame to prepare for dotplot visualization.

make_dot_plot = function(rna_seurat,gene_list){
  feature_rna_plot = rna_seurat@meta.data
  gene_df = as.data.frame(t(rna_seurat@assays$RNA@data[(rownames(rna_seurat@assays$RNA@data) %in% gene_list),]))
  
  feature_rna_plot = cbind(feature_rna_plot,gene_df)
  
  # Calculate ave_expressing values for each cluster separately and pct for each cluster separately
  cluster_gene_df = matrix(NA,nrow = length(unique(feature_rna_plot$seurat_clusters)),ncol = 2*ncol(gene_df))
  clusters = c(0:10)
  
  for(i in 1:length(clusters)){
    subcluster_df = subset(feature_rna_plot,feature_rna_plot$seurat_clusters == clusters[i])
    for(j in 8:ncol(feature_rna_plot)){
      cluster_gene_df[i,j-7] = sum(subcluster_df[,j])/nrow(subcluster_df)
      cluster_gene_df[i,(j-7+ncol(gene_df))] = length(which(subcluster_df[,j] != 0))/nrow(subcluster_df) * 100
    }
  }
  
  cluster_gene_df = as.data.frame(cluster_gene_df)
  rownames(cluster_gene_df) = as.character(clusters)
  colnames(cluster_gene_df)[1:ncol(gene_df)] = colnames(feature_rna_plot)[8:(8+ncol(gene_df)-1)]
  colnames(cluster_gene_df)[((ncol(gene_df))+1): (2*ncol(gene_df))] = paste("pct",colnames(feature_rna_plot)[8:(8+ncol(gene_df)-1)],sep = "_")
  
  # Convert the wide format to long format
  ave_expression_df = cluster_gene_df[,1:ncol(gene_df)]
  ave_expression_df$cluster = clusters
  
  pct_df = cluster_gene_df[,(ncol(gene_df)+1): (2*ncol(gene_df))]
  pct_df$cluster = clusters
  colnames(pct_df) = colnames(ave_expression_df)
  
  return_df1 = gather(ave_expression_df,gene,ave_expression,(colnames(ave_expression_df)[1]):(colnames(ave_expression_df)[ncol(gene_df)]))
  return_df2 = gather(pct_df,gene,pct_gene,(colnames(ave_expression_df)[1]):(colnames(ave_expression_df)[ncol(gene_df)]))
  
  return_df = left_join(return_df1,return_df2,by = c("cluster","gene"))
  
  # return a data frame
  return(return_df)
}


## ChIP-seq peak file annotation
#' The function is to annotate ChIP-seq peak file to the nearest genes based on the defined window size. 
#' @param chipseq_data A ChIP-seq data in .bed or .bedGraph, or in .txt format. 
#' @param upstream_window The size of upstream window of the TSS. For example, a peak within the 1kb upstream of the TSS will be 
#' annotated to the gene. Otherwise, it will not be annotated to the gene. 
#' @param downstream_window The same as upstream_window. But downstream_window is the size of downstream window of the TSS.
#' @return A list of ChIP-seq peak annotation results. The first element is the gene list to present all annotated genes within 
#' the defined windows. The second element is a data frame containing details of the annotated genes, like distance to the TSS, gene
#' symbol, entrezIDs, etc. 

chipseq_annotation_sox2 = function(chipseq_data,upstream_window = 1000,downstream_window = 1000){
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(GenomicRanges)
  
  regions = makeGRangesFromDataFrame(chipseq_data,keep.extra.columns = T,seqnames.field = "chrom",start.field = "start",end.field = "end",strand.field = "strand")
  regionAnno = as.data.frame(annotatePeak(regions, TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, annoDb = "org.Mm.eg.db"))
  
  # Consider the window threshold
  promoter_regions = subset(regionAnno,(regionAnno$distanceToTSS <= downstream_window & regionAnno$distanceToTSS >= 0) | (regionAnno$distanceToTSS >= upstream_window & regionAnno$distanceToTSS < 0))
  output_gene_list = promoter_regions$SYMBOL
  
  return(list(output_gene_list,promoter_regions))
}
