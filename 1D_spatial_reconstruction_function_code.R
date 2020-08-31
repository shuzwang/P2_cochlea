#### ID Spatial Reconstruction Functional Code
### This code is to reconstruct 1D spatial map using scRNA-seq and scTAC-seq data. 

### Quantile normalization
#' @param df A feature-by-sample matrix (e.g. gene x cell).
#' @return A normalized feature-by-sample matrix after quantile normalization.

quantile_normalisation = function(df){
  df_rank = apply(df,2,rank,ties.method="min")
  df_sorted = data.frame(apply(df, 2, sort))
  df_mean = apply(df_sorted, 1, mean)
  
  index_to_mean = function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final = apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) = rownames(df)
  return(df_final)
}

## Rotate the data based on centriod vector
#' Centroid for two clusters (Apex and Base) are calculated in Euclidean space. Then a vector from base centriod to apex centriod
#' is generated. The 2D coordinate system is rotated based on the base-to-apex vector, with apex facing up. 
#' 
#' @param data A data frame with cells in rows and PC1 and PC2 coordinates in columns. For instance, the data frame is the second
#' element of projection_cell function. 
#' @param norm Normlization method to order the cells. rank normalization arranges the cells based on PC1. quantile normalization
#' assumes PC_1 and PC_2 across cells follow same distribution.
#' @return A rank normlized data frame for further analysis and visualization.

rotate_data=function(data,norm = c("rank","quantile")){
  # calculate the centriod
  data$Experiment = factor(data$Experiment,levels = c("Apex","Base"))
  centriod = data %>% group_by(Experiment) %>% summarise(ave_v1 = mean(PC_1),ave_v2 = mean(PC_2)) %>% as.matrix()
  
  c = dist(centriod[,2:3], method = "euclidean")
  a = as.numeric(centriod[1,2:3])
  b = as.numeric(centriod[2,2:3])
  
  y_vec = c(0,1)
  
  # Vector projection to y-vec
  y_proj = abs(y_vec %*% (a-b))
  
  # Generate rotation matrix (pi+theta)
  cosine_theta = y_proj/c
  sine_theta = sqrt(c^2 - y_proj^2)/c
  
  if(a[1] >= b[1] & a[2] >= b[2]){
    Q_matrix = matrix(c(cosine_theta,-sine_theta,sine_theta,cosine_theta),nrow = 2, ncol = 2)
    
    # Rotate the data
    R_matrix = as.matrix(data[,1:2]) %*% Q_matrix
    R_matrix = as.data.frame(R_matrix)
    
  }else if(a[1] >= b[1] & a[2] < b[2]){
    #Q_matrix = matrix(c(-cosine_theta,sine_theta,-sine_theta,-cosine_theta),nrow = 2, ncol = 2)
    Q_matrix = matrix(c(-cosine_theta,sine_theta,-sine_theta,-cosine_theta),nrow = 2, ncol = 2)
    # Rotate the data
    R_matrix = as.matrix(data[,1:2]) %*% t(Q_matrix)
    R_matrix = as.data.frame(R_matrix)
  }else if(a[1] < b[1] & a[2] >= b[2]){
    Q_matrix = matrix(c(cosine_theta,sine_theta,-sine_theta,cosine_theta),nrow = 2, ncol = 2)
    
    # Rotate the data
    R_matrix = as.matrix(data[,1:2]) %*% (Q_matrix)
    R_matrix = as.data.frame(R_matrix)
  }else if(a[1] < b[1] & a[2] < b[2]){
    Q_matrix = matrix(c(-cosine_theta,sine_theta,-sine_theta,-cosine_theta),nrow = 2, ncol = 2)
    
    # Rotate the data
    R_matrix = as.matrix(data[,1:2]) %*% Q_matrix
    R_matrix = as.data.frame(R_matrix)
  }
  
  # Rank/Quantile normalization
  if(norm == "rank"){
    # Rank the cells based on PC_1, which contains spatial information 
    R_matrix$Rank = rank(R_matrix$V2)
    Rank_normalized_matrix = cbind(R_matrix,data[,3:ncol(data)])
  }else if(norm == "quantile"){
    # Assume PC_1 and PC_2 across cells follow same distribution
    Rank_normalized_matrix = quantile_normalisation(R_matrix)
    Rank_normalized_matrix = cbind(Rank_normalized_matrix,data[,3:ncol(data)])
  }
  
  return(Rank_normalized_matrix)
}

## Project the cells along apex to base using seurat object. 
#' To project the cells along tonotopic axis, differential features (e.g. DEGs and DARs) between apex and base are identified. 
#' Then PCA is applied for the DEGs/DARs. 2D PCA coordinate system is setup for further analysis (e.g. rotate_data function).
#' @param object A Seurat object for one specific cluster
#' @param meta A gene list which will be kept in output for visualization.Each gene will be shown in a individual colunmn.
#' @param cell_type A character to add annotations for each cell (e.g. HC, PC/DC).
#' @param npcs The number of selected PCs in PCA analysis. By default, it is 50 PCs.
#' @param is_positive A boolean parameter for differential analysis. If TRUE, only return positive markers in FindAllMarkers function.
#' @param logFC_threshold A threshold to define differentially expressed genes. By default, logFC_threshold is 0.25 that genes at least
#' 0.25-fold difference between apex and base are kept.
#' @param p_threshold P_value threshold of differential analysis. By default, p_threshold is 1e-2. 
#' @return A list with three elements. The first element is a Seurat object. The second element is a data frame with PCA 
#' coordinates, which is the input for rotate_data function. The third element is a data frame with DEGs between apex and
#' base.

projection_cell = function(object,meta,cell_type,npcs = 50,is_positive=TRUE,logFC_threshold = 0.25,p_threshold = 1e-2){
  # Convert the active ident to Experiment (apex and base)
  object@meta.data$Experiment = rownames(object@meta.data)
  object@meta.data$Experiment = gsub(".*-","",object@meta.data$Experiment)
  
  new_active_ident = object@meta.data$Experiment
  new_active_ident = as.factor(new_active_ident)
  names(new_active_ident) = rownames(object@meta.data)
  levels(new_active_ident)[levels(new_active_ident)==1] = "Apex"
  levels(new_active_ident)[levels(new_active_ident)==2] = "Base"
  object@active.ident = new_active_ident
  #object@active.ident = factor(object@active.ident,levels = c("Apex","Base"))
  
  # Find DE genes between apex and base
  de_genes = FindAllMarkers(object,only.pos = is_positive, min.pct = 0.25, logfc.threshold = logFC_threshold,test.use = "wilcox",return.thresh = p_threshold)
  de_genes = subset(de_genes,de_genes$p_val <= 0.005)
  object = ScaleData(object,features = de_genes$gene)
  mito_gene = de_genes$gene[grepl("^mt",de_genes$gene)]
  de_gene = setdiff(de_genes$gene,mito_gene)
  object = RunPCA(object, features = de_gene,npcs = npcs,verbose = F)
  pca_df = as.data.frame(object@reductions$pca@cell.embeddings)
  pca_df = pca_df[,1:2]
  
  # Add meta data
  pca_df$Experiment = new_active_ident
  pca_df$CellType = cell_type
  if(length(meta) > 0){
    for (i in 1:length(meta)) {
      gene = meta[i]
      pca_df[,gene] = object@assays$RNA@data[gene,]
    }
  }
  
  
  return(list(object,pca_df,de_genes))
}


## Make a single curtain plot
#' @param rotated_df A data frame with rank for each individual cell. For instance, the output of rotate_data function. 
#' @param feature A feature (e.g. gene, DAR) name to show in the single curtain plot. 
#' @param seed Set seed for geom_jitter. By default, seed is 10.
#' @param ratio The aspect ratio of x-axis and y-axis. By default, ratio is 0.05.
#' @param experiment A boolean parameter to show the meta-data, Library ID, in single curtain plot. If FALSE, the value of the
#' feature will be shown in single curtain plot. By default, experiment is TRUE. 
#' @param color_fill The color of strip background. By default, color_fill is grey. 
#' @param legend_position The position of legends. By default, legend_position is right. There are other options: left, none.
#' @return A single curtain plot.

make_single_curtain_plot = function(rotated_df,feature,seed = 10,ratio = 0.05,experiment = T,color_fill = "grey",legend_position = "right"){
  
  if(experiment == T){
    # Set seed to make jitter plot reproducible
    set.seed(seed)
    
    p1 = ggplot(rotated_df,aes(0,Rank,color = Experiment)) + 
      geom_jitter(size = 3,alpha = 1) + 
      coord_fixed(ratio = ratio) +
      scale_colour_manual(values = c("#ec1c24","#2e3191"))+
      theme_bw() + 
      theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),legend.position = legend_position)+
      facet_grid(~CellType)+
      theme(strip.text.x = element_text(size = 25),strip.background = element_rect(fill = color_fill))
  }else{
    
    feature_vec = rotated_df[,feature]
    feature_vec[which(feature_vec == 0)] = NA
    rotated_df[,feature] = feature_vec
    
    rotated_df$label = feature
    
    mid = mean(rotated_df[,feature],na.rm = T)
    # Set seed to make jitter plot reproducible
    set.seed(seed)
    p1 = ggplot(rotated_df,aes(0,Rank,color = rotated_df[,feature])) + 
      geom_jitter(size = 3) + 
      coord_fixed(ratio = ratio) +
      scale_color_gradient2(low = "white",high = "#E65948",mid = "#FEF6B1",midpoint = mid,na.value = "grey80")+
      theme_bw() + 
      theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),legend.text = element_text(size = 15),legend.position = legend_position) +
      labs(shape = "Experiment",color = feature)+
      facet_grid(~label)+
      theme(strip.text.x = element_text(size = 25),strip.background = element_rect(fill = color_fill))
  }
  
  return(p1)
}