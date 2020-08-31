##### Gene Regulatory Network (GRN) Inference Function Code
### A 3-step pipeline to reconstruct GRN. 

### Step1: identify co-expression modules using GENIE3
#' GENIE3 decomposes the network into a lot of regressions. For each regression, thegene expression level of a TF 
#' is predicted from the remaining genes using tree-basedmethods such as Random Forest or Extra-Trees.

#' @param gene_by_cell_matrix A scRNA-seq matrix (gene x cells). The cell barcodes include only two cell types for comparison. 
#' @param gene_list A list of interested genes
#' @param motif_list A list of TFs/motifs. Each TF/motif will be treated as dependent variables in regression model. 
#' @param gene_anno (Optional) Annotations of genes to represent different cell types (clusters) shown in columns of heatmap
#' @param motif_anno (Optional) Annotations of motifs to represent cell type specific TFs shown in rows of heatmap
#' @param ncore Number of core to use running GENIE3. By default, ncore is 2.
#' @param link_threshold The threshold of importance score from Random forest in GENIE3. By default, we only keep the TF-target gene 
#' links with importance score greater than 0.001. 
#' @param corr_threshold The threshold of correlation between any pair of genes is to define positive correlation and negative correlation. 
#' By default, corr_threshold is 0.03. The correlation less than 0.03 is defined as 0. 
#' @param filter_gene_thres Number of expressed gene in each cell. By default, filter_gene_thres is 1. 
#' @return A list of results. The first element is the finalized matrix (gene x cell) as input for GENIE3 algorithm. 
#' The second element is gene correlation squared matrix. The third element is a long matrix containing TF-target gene pairs
#' passed the link threshold. The fourth element is a wide matrix (gene x TF) for visualization. The fifth element is a 
#' heatmap with TFs in the rows and genes in the columns.

single_cell_GENIE3 = function(gene_by_cell_matrix = rna_matrix,gene_list,motif_list,gene_anno=NULL,motif_anno=NULL,ncore = 2,link_threshold = 0.001,corr_threshold = 0.03,filter_gene_thres = 1){
  library(GENIE3)
  library(tidyr)
  library(dplyr)
  library(pheatmap)
  library(Hmisc)
  
  # Prepare inputs for GENIE3 
  gene_tf_list = unique(c(gene_list,motif_list))
  input_matrix = gene_by_cell_matrix[which(rownames(gene_by_cell_matrix) %in% gene_tf_list),]
  input_matrix = as.matrix(input_matrix)
  
  # Remove genes that only express in a few cells 
  gene_expressed_list = rowSums(input_matrix != 0)
  gene_kept = rownames(input_matrix)[which(gene_expressed_list >= filter_gene_thres)]
  input_matrix = input_matrix[which(rownames(input_matrix) %in% gene_kept),]
  
  inputTFs = intersect(motif_list,rownames(input_matrix))
  
  # Correlation analysis
  corrMat = cor(t(input_matrix), method="spearman")
  
  # GENIE3
  set.seed(10)
  weightMatrix = GENIE3::GENIE3(input_matrix, regulators=inputTFs, nCores=ncore, targets=rownames(input_matrix), treeMethod = "RF",verbose = T)
  linkList_df = GENIE3::getLinkList(weightMatrix,threshold = link_threshold)
  
  # Add correlation into the link data frame
  tf_gene_corrMat = corrMat[,which(colnames(corrMat) %in% inputTFs)]
  tf_gene_corrMat = as.data.frame(tf_gene_corrMat)
  tf_gene_corrMat$targetGene = rownames(tf_gene_corrMat)
  tf_gene_corrDf = gather(tf_gene_corrMat,regulatoryGene,correlation,1:(ncol(tf_gene_corrMat)-1),factor_key=FALSE)
  tf_gene_corrDf$correlation = as.numeric(tf_gene_corrDf$correlation >= corr_threshold) - as.numeric(tf_gene_corrDf$correlation <= -corr_threshold) 
  
  linkList_df$regulatoryGene = as.character(linkList_df$regulatoryGene)
  linkList_df$targetGene = as.character(linkList_df$targetGene)
  corr_link_df = left_join(linkList_df,tf_gene_corrDf,by = c("regulatoryGene","targetGene"))
  
  # Generate updated weights to include the direction of regulation
  corr_link_df$mode_weight = corr_link_df$weight * corr_link_df$correlation
  
  # Generate heatmap to show the co-expression results
  gene_tf_visual_df = spread(corr_link_df[,c(1,2,5)],regulatoryGene,mode_weight,fill=0)
  rownames(gene_tf_visual_df) = gene_tf_visual_df$targetGene
  gene_tf_visual_df = gene_tf_visual_df[,-1]
  
  # Add annotation on x-axis and y-axis
  if(is.null(gene_anno) == FALSE & is.null(motif_anno) == FALSE){
    p1 = pheatmap(t(gene_tf_visual_df),scale = "none",annotation_row = motif_anno,annotation_col = gene_anno,fontsize = 4)
  }else if(is.null(gene_anno) == FALSE & is.null(motif_anno) == TRUE){
    p1 = pheatmap(t(gene_tf_visual_df),scale = "none",annotation_col = gene_anno,fontsize = 4)
  }else if(is.null(gene_anno) == TRUE & is.null(motif_anno) == FALSE){
    p1 = pheatmap(t(gene_tf_visual_df),scale = "none",annotation_row = motif_anno,fontsize = 4)
  }else{
    p1 = pheatmap(t(gene_tf_visual_df),scale = "none",fontsize = 4)
  }
  
  return(list(input_matrix,corrMat,corr_link_df,gene_tf_visual_df,p1))
}

## Step2: integrate scATAC-seq peak information 
#' To identify direct target genes, scATAC-seq data has been used. Specifically, FIMO has been applied for each cluster, separately.
#' The TF-target gene links which lacking putative TF binding sites within a fixed window of the target gene TSS were filtered out.

#' @param fimo_results FIMO output files for two clusters in .txt format.
#' @param translation_df A translation table which includes motif IDs and corresponding gene names.
#' @param gene_links The TF-target gene link data frame, which is the thrid element of step1 output. 
#' @param window_size The fixed window size. By default, window_size is 50kb, which is 50kb upstream of the TSS and 50kb downstream the TSS.
#' In total, the fixed window size is 100kb. 
#' @param gene_anno (Optional) Annotations of genes to represent different cell types (clusters) shown in columns of heatmap
#' @param motif_anno (Optional) Annotations of motifs to represent cell type specific TFs shown in rows of heatmap
#' @param consider_negative_correlation A boolean parameter. If TRUE, negative correlation will be considered. If FALSE, only positive correlation
#' will be considered. 
#' @param negative_correlation_filter A boolean parameter. If TRUE, only negative correlated TF-target gene links will be kept.
#' @return A list of results. The first element is the finalized matrix (gene x TF) after removing TF-target genes links which lacking
#'  putative TF binding sites within a fixed window. The second element is is a heatmap with TFs in the rows and genes in the columns.
#'  The third element is a motif list only including self-regulating genes.

run_atacseq_peak = function(fimo_results = fimo_results,translation_df = translation_table,gene_links = link_data,window_size = 50000,motif_anno,gene_anno,consider_negative_correlation = TRUE,negative_correlation_filter = TRUE){
  library(tidyr)
  library(dplyr)
  library(pheatmap)
  
  fimo_df = fimo_results[,c("motif_id","cluster","annotation","geneId","distanceToTSS","SYMBOL")]
  
  # Filter TFBS based on window size around the TSS
  fimo_tss1 = subset(fimo_df,fimo_df$distanceToTSS >= 0 & fimo_df$distanceToTSS <= window_size)
  fimo_tss2 = subset(fimo_df,fimo_df$distanceToTSS < 0 & fimo_df$distanceToTSS >= -window_size)
  
  fimo_tss_combined = rbind(fimo_tss1,fimo_tss2)
  
  # Extract motif name from motif_id
  fimo_tss_combined$motif_id = gsub("_MOUSE.H10MO.*","",fimo_tss_combined$motif_id)
  colnames(fimo_tss_combined)[1] = "HOCO" 
  fimo_tss = left_join(fimo_tss_combined,translation_table,by = "HOCO")
  
  # Generate TF-target gene data frame
  fimo_tss = fimo_tss[,c("gene_name","SYMBOL")]
  colnames(fimo_tss) = c("regulatoryGene","targetGene")
  
  # Subset the positive correlation links from GENIE3
  link_data_pos = subset(gene_links,gene_links$correlation == 1)
  link_data_neg = subset(gene_links,gene_links$correlation == -1)
  
  # Integrate the weight matrix with fimo results
  peak_link = left_join(fimo_tss,link_data_pos[,1:3],by=c("regulatoryGene","targetGene"))
  peak_link = peak_link[complete.cases(peak_link),]
  peak_link = peak_link[!duplicated(peak_link),]
  
  # Consider negative correlation (repression mode)
  if(consider_negative_correlation == TRUE){
    if(negative_correlation_filter != TRUE){
      link_data_neg = link_data_neg[,c(1,2,5)]
      colnames(link_data_neg) = c("regulatoryGene","targetGene","weight")
      peak_link = rbind(peak_link,link_data_neg)
    }
    
    if(negative_correlation_filter == TRUE){
      peak_link_neg = left_join(fimo_tss,link_data_neg[,c(1,2,5)],by=c("regulatoryGene","targetGene"))
      peak_link_neg = peak_link_neg[complete.cases(peak_link_neg),]
      peak_link_neg = peak_link_neg[!duplicated(peak_link_neg),]
      colnames(peak_link_neg) = c("regulatoryGene","targetGene","weight")
      peak_link = rbind(peak_link,peak_link_neg)
    }
  }
  
  # Keep self-regulation
  fimo_tss_hoco = left_join(fimo_tss_combined,translation_table,by = "HOCO")
  fimo_tss_hoco = fimo_tss_hoco[which(fimo_tss_hoco$SYMBOL == fimo_tss_hoco$gene_name),]
  motif_self_regulation_list = unique(fimo_tss_hoco$SYMBOL)
  
  # Generate gene-by-TF matrix
  output_matrix = spread(peak_link,regulatoryGene,weight,fill=0)
  rownames(output_matrix) = output_matrix$targetGene
  output_matrix = output_matrix[,-1]
  
  # Visualization
  if(is.null(gene_anno) == FALSE & is.null(motif_anno) == FALSE){
    p1 = pheatmap(t(output_matrix),scale = "none",annotation_row = motif_anno,annotation_col = gene_anno,fontsize = 4)
  }else if(is.null(gene_anno) == FALSE & is.null(motif_anno) == TRUE){
    p1 = pheatmap(t(output_matrix),scale = "none",annotation_col = gene_anno,fontsize = 4)
  }else if(is.null(gene_anno) == TRUE & is.null(motif_anno) == FALSE){
    p1 = pheatmap(t(output_matrix),scale = "none",annotation_row = motif_anno,fontsize = 4)
  }else{
    p1 = pheatmap(t(output_matrix),scale = "none",fontsize = 4)
  }
  
  return(list(output_matrix,p1,motif_self_regulation_list))
}

## Step2: Identify regulons
#' Identify activator regulons (downstream target genes) and repressor regulons. 

#' @param step2_output_matrix A gene-by-TF matrix identified from run_atacseq_peak function. 
#' @param min_thres The minimum number of downstream target genes to define a regulon. By default, a TF regulon with less than 10 target genes 
#' will be removed for further analysis. 
#' @param max_thres The maximun number of downstream target genes to define a regulon. By default, a TF regulon with greater than 900 target genes 
#' will be removed for further analysis. 
#' @param only_repression A boolean parameters to determine ativator or repressor regulons . By default, only_repression is FALSE which identify 
#' activator regulons. If TRUE, only identified repressor regulons.
#' @param self_regulation_list (Optional) By default, auto-regulation will not be considered. However, with an input of a TF list, auto-regulation
#' could be considered. The input list could be determined from the third element of run_atacseq_peak function. 

#' @return A list of TF regulons. For each TF, it shows the downstream target genes. 

run_regulon = function(step2_output_matrix,min_thres = 10, max_thres = 900,only_repression = FALSE,self_regulation_list = NULL){
  
  # Generate regulon list
  if(only_repression == FALSE){
    regulon_list = list()
    gene_length = apply(step2_output_matrix, 2, function(c) sum(c > 0))
    gene_length_index = which(gene_length  >= min_thres & gene_length <= max_thres)
    
    expr_matrix = step2_output_matrix[,gene_length_index]
    
    # Check whether there is only 1 column in the matrix
    if(ncol(expr_matrix) == 1){
      
      expr_matrix = as.matrix(expr_matrix)
      rownames(expr_matrix) = rownames(step2_output_matrix)
      colnames(expr_matrix) = colnames(step2_output_matrix)
      
      if(is.null(self_regulation_list) == TRUE){
        gene_index = which(expr_matrix[,1] > 0)
        
        regulon_list[[1]] = rownames(expr_matrix)[gene_index]
        names(regulon_list)[[1]] = colnames(expr_matrix)[1]
      }
      if(is.null(self_regulation_list) != TRUE){
        gene_index = which(expr_matrix[,1] > 0)
        
        regulon_list[[1]] = rownames(expr_matrix)[gene_index]
        if((colnames(expr_matrix))[1] %in% self_regulation_list == TRUE){
          regulon_list[[1]][length(regulon_list[[1]]) +1] = (colnames(expr_matrix))[1]
        }
        names(regulon_list)[[1]] = colnames(expr_matrix)[1]
      }
    }else{
      for(i in 1:ncol(expr_matrix)){
        if(is.null(self_regulation_list) == TRUE){
          gene_index = which(expr_matrix[,i] > 0)
          
          regulon_list[[i]] = rownames(expr_matrix)[gene_index]
          names(regulon_list)[[i]] = colnames(expr_matrix)[i]
        }
        if(is.null(self_regulation_list) != TRUE){
          gene_index = which(expr_matrix[,i] > 0)
          
          regulon_list[[i]] = rownames(expr_matrix)[gene_index]
          if((colnames(expr_matrix))[i] %in% self_regulation_list == TRUE){
            regulon_list[[i]][length(regulon_list[[i]]) +1] = (colnames(expr_matrix))[i]
          }
          names(regulon_list)[[i]] = colnames(expr_matrix)[i]
        }
      }
    }
    
  }
  
  if(only_repression == TRUE){
    
    regulon_list = list()
    gene_length = apply(step2_output_matrix, 2, function(c) sum(c < 0))
    gene_length_index = which(gene_length  >= min_thres & gene_length <= max_thres)
    
    expr_matrix = step2_output_matrix[,gene_length_index]
    
    if(ncol(expr_matrix) == 1){
      
      expr_matrix = as.matrix(expr_matrix)
      rownames(expr_matrix) = rownames(step2_output_matrix)
      colnames(expr_matrix) = colnames(step2_output_matrix)
      
      if(is.null(self_regulation_list) == TRUE){
        gene_index = which(expr_matrix[,1] < 0)
        
        regulon_list[[1]] = rownames(expr_matrix)[gene_index]
        names(regulon_list)[[1]] = colnames(expr_matrix)[1]
      }
      
      if(is.null(self_regulation_list) != TRUE){
        gene_index = which(expr_matrix[,1] < 0)
        
        regulon_list[[1]] = rownames(expr_matrix)[gene_index]
        if((colnames(expr_matrix))[1] %in% self_regulation_list == TRUE){
          regulon_list[[1]][length(regulon_list[[1]]) +1] = (colnames(expr_matrix))[1]
        }
        names(regulon_list)[[1]] = colnames(expr_matrix)[1]
      }
    }else{
      for(i in 1:ncol(expr_matrix)){
        if(is.null(self_regulation_list) == TRUE){
          gene_index = which(expr_matrix[,i] < 0)
          
          regulon_list[[i]] = rownames(expr_matrix)[gene_index]
          names(regulon_list)[[i]] = colnames(expr_matrix)[i]
        }
        
        if(is.null(self_regulation_list) != TRUE){
          gene_index = which(expr_matrix[,i] < 0)
          
          regulon_list[[i]] = rownames(expr_matrix)[gene_index]
          if((colnames(expr_matrix))[i] %in% self_regulation_list == TRUE){
            regulon_list[[i]][length(regulon_list[[i]]) +1] = (colnames(expr_matrix))[i]
          }
          names(regulon_list)[[i]] = colnames(expr_matrix)[i]
        }
        
      }
    }
    
  }
  
  return(regulon_list)
}

## Step3: Regulon enrichment for each indivual cell using AUCell
#' Calculate regulon enrichment score for each individual cell using AUCell package. It uses the "Area Under the Curve" (AUC) of the recovery curve 
#' to determine the enrichment of regulons for individual cells. 

#' @param regulon_list A list of TF regulons. For each TF, it shows the downstream target genes. For, example, the output of run_regulon function.
#' @param expression_matrix A matrix (gene x cell) from scRNA-seq data. Usually, the output of first element of single_cell_GENIE3 function was used.
#' @param set_seed Set seed. By default, seed is 123. 

#' @return A matrix (TF regulon x cell) with AUC enrichment scores. 

aucell_function = function(regulon_list,expression_matrix,set_seed = 123){
  library(AUCell)
  
  set.seed(set_seed)
  # Step1: Build gene-expression ranks for each individual cells
  aucellRankings = AUCell_buildRankings(expression_matrix, nCores=1, plotStats=TRUE,verbose = FALSE)
  
  # Step2: Calculate AUC
  regulonAUC = AUCell_calcAUC(regulon_list, aucellRankings, aucMaxRank=aucellRankings@nGenesDetected["1%"], nCores=1)
  
  variableRegulons = names(which(apply(getAUC(regulonAUC), 1, sd) > 0))
  reguDist = as.dist(1-cor(t(getAUC(regulonAUC)[variableRegulons,]), method="spear"))
  reguClust = hclust(reguDist, method="ward.D2")
  regulonClusters = setNames(dynamicTreeCut::cutreeDynamic(reguClust, distM=as.matrix(reguDist), verbose = FALSE), reguClust$labels)
  regulonOrder = reguClust$labels[reguClust$order]
  regulonOrder = regulonOrder[order(regulonClusters[regulonOrder], decreasing = TRUE)]
  
  regulonAUC = regulonAUC[regulonOrder,]
  
  # Step3: Default threshold
  cells_AUCellThresholds = NULL
  cells_AUCellThresholds = AUCell_exploreThresholds(regulonAUC, smallestPopPercent=0.25,assignCells=TRUE, plotHist=FALSE, verbose=FALSE, nCores=1)
  regulonsCells = getAssignments(cells_AUCellThresholds)

  return(getAUC(regulonAUC))
}
