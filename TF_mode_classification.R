#### TF mode classification Functional Code
### This code is to classify the mode of TFs into activators and repressors based on the relation between scRNA-seq and scATAC-seq.


### Identify differentially expressed TF genes between 2 populations using scRNA-seq data
#' @param rna_seurat A Seurat object for scRNA-seq data.
#' @param comparison_group Two groups (populations) you want to compare. For example, we want to compare HCs and PC/DC/IBCs. 
#' @param run_test The test of differential analysis between 2 groups. student t-test, wilcoxon sum rank test, and median value 
#' comparison are provided. 
#' @param input_genes (Optional) A list of genes you want to concentrate on to speed the calculations. For example, here only TF
#' genes are of interest. By default, input_genes is NULL, which includes all input genes for differential analysis.        
#' @return A data frame of differential analysis results. It includes mean values of each group, P_values, log2FC, BH, target1_results, 
#' and target2_results. "target1_results" is a binary variable with 1 represents group1 is significantly different (highly expressed) 
#' compared with group2.

run_regulation_mode_scRNA = function(rna_seurat,comparison_group = c("HC","PC/DC/IBC"),run_test = c("t-test","median","wilcox"),input_genes = NULL){
  library(Seurat)
  library(dplyr)
  library(tidyr)
  
  # Subset the cells based on the comparison group
  target1_seurat = SubsetData(rna_seurat,ident.use = comparison_group[1])
  target2_seurat = SubsetData(rna_seurat,ident.use = comparison_group[2])
  
  # Subset the genes
  if(is.null(input_genes) == TRUE){
    target1_expression_matrix = target1_seurat@assays$RNA@data
    target2_expression_matrix = target2_seurat@assays$RNA@data
  }
  
  if(is.null(input_genes) == FALSE){
    target1_expression_matrix = target1_seurat@assays$RNA@data[which(rownames(target1_seurat@assays$RNA@data) %in% input_genes),]
    target2_expression_matrix = target2_seurat@assays$RNA@data[which(rownames(target2_seurat@assays$RNA@data) %in% input_genes),]
  }
  
  # Remove non-expression TF genes before running tests
  expression_matrix = cbind(target1_expression_matrix,target2_expression_matrix)
  gene_length = apply(expression_matrix, 1, function(c) sum(c != 0))
  gene_length_index = which(gene_length  > 1 )
  
  target1_expression_matrix = target1_expression_matrix[gene_length_index,]
  target2_expression_matrix = target2_expression_matrix[gene_length_index,]
  
  # Calculate the log2FC first
  target1_expression_matrix = as.matrix(t(target1_expression_matrix))
  target2_expression_matrix = as.matrix(t(target2_expression_matrix))
  
  log2FC_1 = log2(x = colMeans(x = expm1(x = target1_expression_matrix)) + 1)
  log2FC_2 = log2(x = colMeans(x = expm1(x = target2_expression_matrix)) + 1)
  log2FC = log2FC_2 - log2FC_1
  
  # Compare two groups using statistical tests
  output_df = data.frame()
  
  if(run_test == "median"){
    for(i in 1:ncol(target1_expression_matrix)){
      output_df[i,1] = median(target1_expression_matrix[,i])
      output_df[i,2] = median(target2_expression_matrix[,i])
      output_df[i,3] = sum((output_df[i,1]>output_df[i,2]) - (output_df[i,1]<output_df[i,2]))
      output_df[i,4] = sum((output_df[i,1]<output_df[i,2]) - (output_df[i,1]>output_df[i,2]))
    }
    output_df[,5] = log2FC
    
    colnames(output_df) = c("target1_median","target2_median","target1_results","target2_results","log2FC")
    rownames(output_df) = colnames(target1_expression_matrix)
    
  }
  
  if(run_test == "t-test"){
    for(i in 1:ncol(target1_expression_matrix)){
      t_results = t.test(target1_expression_matrix[,i],target2_expression_matrix[,i],alternative = "two.sided",paired = FALSE,var.equal = F)
      output_df[i,1] = t_results$estimate[1]
      output_df[i,2] = t_results$estimate[2]
      output_df[i,3] = t_results$p.value
    }
    output_df[,4] = log2FC
    output_df[,5] = p.adjust(output_df[,3],method = "BH")
    
    for(i in 1:nrow(output_df)){
      if(output_df[i,1] > output_df[i,2] & output_df[i,5] <= 0.05){
        output_df[i,6] = 1
        output_df[i,7] = -1
      }
      else if(output_df[i,1] < output_df[i,2] & output_df[i,5] <= 0.05){
        output_df[i,6] = -1
        output_df[i,7] = 1
      }
      else{
        output_df[i,6] = 0
        output_df[i,7] = 0
      }
    }
    colnames(output_df) = c("target1_mean","target2_mean","t_test_pvalue","log2FC","BH","target1_results","target2_results")
    rownames(output_df) = colnames(target1_expression_matrix)
    
  }
  
  if(run_test == "wilcox"){
    for(i in 1:ncol(target1_expression_matrix)){
      wilcox_results = wilcox.test(target1_expression_matrix[,i],target2_expression_matrix[,i],alternative = "two.sided",paired = FALSE)
      output_df[i,1] = mean(target1_expression_matrix[,i])
      output_df[i,2] = mean(target2_expression_matrix[,i])
      output_df[i,3] = wilcox_results$p.value
    }
    output_df[,4] = log2FC
    output_df[,5] = p.adjust(output_df[,3],method = "BH")
    
    for(i in 1:nrow(output_df)){
      if(output_df[i,1] > output_df[i,2] & output_df[i,5] <= 0.05){
        output_df[i,6] = 1
        output_df[i,7] = -1
      }
      else if(output_df[i,1] < output_df[i,2] & output_df[i,5] <= 0.05){
        output_df[i,6] = -1
        output_df[i,7] = 1
      }
      else{
        output_df[i,6] = 0
        output_df[i,7] = 0
      }
    }
    colnames(output_df) = c("target1_mean","target2_mean","wilcox_pvalue","log2FC","BH","target1_results","target2_results")
    rownames(output_df) = colnames(target1_expression_matrix)
  }
  
  return(output_df)
}


### Identify differentially accessible TF motifs between 2 populations using scATAC-seq z-scores
#' @param zscore_matrix A zscore matrix (cell x TF/motif) identified from chromVAR. 
#' @param target1_cell_barcode A list of cell barcodes of group1. 
#' @param target2_cell_barcode A list of cell barcodes of group2. 
#' @param run_test The test of differential analysis between 2 groups. student t-test, wilcoxon sum rank test, and median value 
#' comparison are provided. 
#' @param translation_table A translation table which matches HOCOMOCO v10 motif names with their corresponding gene symbols. The 
#' translation table is provided in the github. Please check HOCOMOCO_v10_gene_translation_table.txt in the github.      
#' @return A data frame of differential analysis results. It includes mean values of each group, P_values, log2FC, BH, target1_results, 
#' and target2_results. "target1_results" is a binary variable with 1 represents group1 is significantly different (highly expressed) 
#' compared with group2.

run_regulation_mode_scATAC = function(zscore_matrix,target1_cell_barcode, target2_cell_barcode, run_test = c("median","t-test","wilcox"),translation_table){
  library(chromVAR)
  library(SnapATAC)
  
  # Figure out the duplicated motifs
  zscore_motif = data.frame(motif = rownames(zscore_matrix))
  zscore_motif$HOCO = gsub("_MOUSE.H10MO.*","",zscore_motif$motif)
  zscore_motif = left_join(zscore_motif,translation_table,by = "HOCO")
  
  motif_index = which(duplicated(zscore_motif$gene_name) == FALSE)
  zscore_matrix = zscore_matrix[motif_index,]
  rownames(zscore_matrix) = zscore_motif$gene_name[motif_index]
  
  # Subset the cells based on target1 and target2
  target1_zscore_matrix = zscore_matrix[,which(colnames(zscore_matrix) %in% target1_cell_barcode)]
  target2_zscore_matrix = zscore_matrix[,which(colnames(zscore_matrix) %in% target2_cell_barcode)]
  
  target1_zscore_matrix = as.data.frame(t(target1_zscore_matrix))
  target2_zscore_matrix = as.data.frame(t(target2_zscore_matrix))
  
  # Calculate the log2FC first
  log2FC_1 = log2(x = colMeans(x = expm1(x = target1_zscore_matrix)) + 1)
  log2FC_2 = log2(x = colMeans(x = expm1(x = target2_zscore_matrix)) + 1)
  log2FC = log2FC_2 - log2FC_1
  
  # Compare two groups using statistical tests
  output_df = data.frame()
  
  if(run_test == "median"){
    for(i in 1:ncol(target1_zscore_matrix)){
      output_df[i,1] = median(target1_zscore_matrix[,i])
      output_df[i,2] = median(target2_zscore_matrix[,i])
      output_df[i,3] = sum((output_df[i,1]>output_df[i,2]) - (output_df[i,1]<output_df[i,2]))
      output_df[i,4] = sum((output_df[i,1]<output_df[i,2]) - (output_df[i,1]>output_df[i,2]))
    }
    output_df[,5] = log2FC
    colnames(output_df) = c("target1_median","target2_median","target1_results","target2_results","log2FC")
    rownames(output_df) = colnames(target1_zscore_matrix)
  }
  
  if(run_test == "t-test"){
    for(i in 1:ncol(target1_zscore_matrix)){
      t_results = t.test(target1_zscore_matrix[,i],target2_zscore_matrix[,i],alternative = "two.sided",paired = FALSE,var.equal = F)
      output_df[i,1] = t_results$estimate[1]
      output_df[i,2] = t_results$estimate[2]
      output_df[i,3] = t_results$p.value
    }
    output_df[,4] = log2FC
    output_df[,5] = p.adjust(output_df[,3],method = "BH")
    
    for(i in 1:nrow(output_df)){
      if(output_df[i,1] > output_df[i,2] & output_df[i,5] <= 0.05){
        output_df[i,6] = 1
        output_df[i,7] = -1
      }
      else if(output_df[i,1] < output_df[i,2] & output_df[i,5] <= 0.05){
        output_df[i,6] = -1
        output_df[i,7] = 1
      }
      else{
        output_df[i,6] = 0
        output_df[i,7] = 0
      }
    }
    colnames(output_df) = c("target1_mean","target2_mean","t_test_pvalue","log2FC","BH","target1_results","target2_results")
    rownames(output_df) = colnames(target1_zscore_matrix)
    
  }
  
  if(run_test == "wilcox"){
    for(i in 1:ncol(target1_zscore_matrix)){
      wilcox_results = wilcox.test(target1_zscore_matrix[,i],target2_zscore_matrix[,i],alternative = "two.sided",paired = FALSE)
      output_df[i,1] = mean(target1_zscore_matrix[,i])
      output_df[i,2] = mean(target2_zscore_matrix[,i])
      output_df[i,3] = wilcox_results$p.value
    }
    output_df[,4] = log2FC
    output_df[,5] = p.adjust(output_df[,3],method = "BH")
    
    for(i in 1:nrow(output_df)){
      if(output_df[i,1] > output_df[i,2] & output_df[i,5] <= 0.05){
        output_df[i,6] = 1
        output_df[i,7] = -1
      }
      else if(output_df[i,1] < output_df[i,2] & output_df[i,5] <= 0.05){
        output_df[i,6] = -1
        output_df[i,7] = 1
      }
      else{
        output_df[i,6] = 0
        output_df[i,7] = 0
      }
    }
    colnames(output_df) = c("target1_mean","target2_mean","wilcox_pvalue","log2FC","BH","target1_results","target2_results")
    rownames(output_df) = colnames(target1_zscore_matrix)
  }
  
  return(output_df)
}


### Classify the TF mode into Activators and repressors
#' @param RNA_output An output from run_regulation_mode_scRNA function. It contains the differential analysis results of scRNA-seq between
#' 2 populations. 
#' @param ATAC_output An output from run_regulation_mode_scATAC function. It contains the differential analysis results of scATAC-seq zscores 
#' between 2 populations. 
#' @return A list of classification results. The first element is a list of all activators and repressors. The second element is a data frame
#' to integrate and determine the relationship between gene expression level and chromatin accessibility. Detailed information could be 
#' traced back in the data frame. 

run_activator_repressor = function(RNA_output,ATAC_output){
  # Find the intersect TF genes
  intersected_TF_genes = intersect(rownames(RNA_output),rownames(ATAC_output))
  
  # Concentenate two data frames
  rna_df = RNA_output[which(rownames(RNA_output) %in% intersected_TF_genes),]
  atac_df = ATAC_output[which(rownames(ATAC_output) %in% intersected_TF_genes),]
  
  rna_df$gene = rownames(rna_df)
  atac_df$gene = rownames(atac_df)
  
  combined_df = left_join(rna_df,atac_df,by = "gene")
  
  # Determine the activators and repressors
  ativators_list = combined_df$gene[which(combined_df$target1_results.x * combined_df$target1_results.y == 1)]
  repressor_list = combined_df$gene[which(combined_df$target1_results.x * combined_df$target1_results.y == -1)]
  
  output_list = list(ativators = ativators_list,repressors = repressor_list)
  
  # Add regulation mode in the data frame
  combined_df$mode = "Undertermined"
  combined_df$mode[which(combined_df$gene %in% ativators_list)] = "Activator"
  combined_df$mode[which(combined_df$gene %in% repressor_list)] = "Repressor"
  
  
  return(list(output_list,combined_df))
}


