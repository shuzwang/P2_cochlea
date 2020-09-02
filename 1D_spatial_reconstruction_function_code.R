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

## Rotate the data based on centroid vector
#' Centroid for two clusters (Apex and Base) are calculated in Euclidean space. Then a vector from base centroid to apex centroid
#' is generated. The 2D coordinate system is rotated based on the base-to-apex vector, with apex facing up. 
#' 
#' @param data A data frame with cells in rows and PC1 and PC2 coordinates in columns. For instance, the data frame is the second
#' element of projection_cell function. 
#' @param norm Normlization method to order the cells. rank normalization arranges the cells based on PC1. quantile normalization
#' assumes PC_1 and PC_2 across cells follow same distribution.
#' @return A rank normlized data frame for further analysis and visualization.

rotate_data=function(data,norm = c("rank","quantile")){
  # calculate the centroid
  data$Experiment = factor(data$Experiment,levels = c("Apex","Base"))
  centroid = data %>% group_by(Experiment) %>% summarise(ave_v1 = mean(PC_1),ave_v2 = mean(PC_2)) %>% as.matrix()
  
  c = dist(centroid[,2:3], method = "euclidean")
  a = as.numeric(centroid[1,2:3])
  b = as.numeric(centroid[2,2:3])
  
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


## Make single curtain plot for IHC and OHC separately
#' @param rotated_hc A data frame with rank for each individual cell including IHCs and OHCs. For instance, the output of rotate_data function. 
#' @param rotated_ihc A data frame with rank for each individual cell only including IHCs.
#' @param rotated_ohc A data frame with rank for each individual cell only including OHCs. 
#' @param point_color A feature will be presented in the curtain plot. If "cluster", it shows the IHC/OHC cluster IDs. If "Experiment", it shows
#' meta data information. If "gene_name", you need to provide gene symbols to show a specific feature expression in the curtain plot. 
#' @param point_shape A feature you want to distinguish in shape of data points. For instance, different shape of data points represent subclusters
#' of HCs (e.g. IHCs, OHCs)
#' @param shift_width In a single curtain plot, two subclusters are separate. The shift_width is to adjust the width of interval between two 
#' subclusters. By default, the shift_width is 0.5. 
#' @param seed Set seed. By default, seed is 10.
#' @param plot_ratio The aspect ratio of x-axis and y-axis. By default, ratio is 0.07.
#' @param plot_legend A boolean parameter to show legends. By default, plot_legend is TRUE to show the legends. 
#' @param mid_value The mid value of gene expression used in scale_color_gradient2 function to define the color scale. "mean" is the average 
#' expression level. "median" is the median of gene expression level across cells. Any fixed numeric number could be used. For example, 
#' mid_value is 1.5. By default, mid_value is "mean".  
#' @param dot_size The size of data points. By default, dot_size is 3. 
#' @return A single curtain plot with IHCs locates on the left and OHCs on the right.

generate_curtain_plot_IHC_OHC = function(rotated_hc,rotated_ihc,rotated_ohc,point_color=c("cluster","Experiment","gene_name"),point_shape, shift_width = 0.5,plot_ratio = 0.07,plot_legend = TRUE,mid_value = "mean",dot_size = 3,seed = 10){
  library(ggplot2)
  
  # Setup jitter position
  set.seed(20)
  jitter_ohc = position_jitter(width = 0.2)
  set.seed(20)
  jitter_ihc = position_jitter(width = 0.1)
  
  # Visualization
  if(point_color == "cluster"){
    rotated_hc$label = "Cell Type"
    if(plot_legend == FALSE){
      set.seed(seed)
      p1 = ggplot(rotated_hc) + geom_point(data = rotated_ohc,aes(shift_width,Rank,color = rotated_ohc[,point_color],shape = rotated_ohc[,point_shape]),position = jitter_ohc,size = dot_size,color = "#CC79A7") +
        geom_point(data = rotated_ihc,aes(0,Rank,color = rotated_ihc[,point_color],shape = rotated_ihc[,point_shape]),position = jitter_ihc,size = dot_size,color = "#E69F00") + 
        coord_fixed(ratio = plot_ratio) + theme_bw() + 
        theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+
        facet_grid(~label)+
        theme(strip.text.x = element_text(size = 25)) +
        theme(legend.title = element_text(size = 30),legend.text = element_text(size = 25),legend.position = "none") +
        labs(shape = "Cluster")
    }
    if(plot_legend == TRUE){
      set.seed(seed)
      p1 = ggplot(rotated_hc) + geom_point(data = rotated_ohc,aes(shift_width,Rank,color = rotated_ohc[,point_color],shape = rotated_ohc[,point_shape]),position = jitter_ohc,size = dot_size,color = "#CC79A7") +
        geom_point(data = rotated_ihc,aes(0,Rank,color = rotated_ihc[,point_color],shape = rotated_ihc[,point_shape]),position = jitter_ihc,size = dot_size,color = "#E69F00") + 
        coord_fixed(ratio = plot_ratio) + theme_bw() + 
        theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+
        facet_grid(~label)+
        theme(strip.text.x = element_text(size = 25)) +
        theme(legend.title = element_text(size = 30),legend.text = element_text(size = 25)) +
        labs(shape = "Cluster")
    }
    
  }
  
  if(point_color == "Experiment"){
    rotated_hc$label = "Library ID"
    if(plot_legend == FALSE){
      set.seed(seed)
      p1 = ggplot(rotated_hc) + geom_point(data = rotated_ohc,aes(shift_width,Rank,color = rotated_ohc[,point_color],shape = rotated_ohc[,point_shape]),position = jitter_ohc,size = dot_size) +
        geom_point(data = rotated_ihc,aes(0,Rank,color = rotated_ihc[,point_color],shape = rotated_ihc[,point_shape]),position = jitter_ihc,size = dot_size) + 
        coord_fixed(ratio = plot_ratio) + theme_bw() + 
        scale_color_manual(values = c("#ec1c24","#2e3191")) +
        theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+
        facet_grid(~label)+
        theme(strip.text.x = element_text(size = 25)) +
        theme(legend.title = element_text(size = 30),legend.text = element_text(size = 25),legend.position = "none") +
        labs(shape = "Cluster",color = "Library ID")
    }
    
    if(plot_legend == TRUE){
      set.seed(seed)
      p1 = ggplot(rotated_hc) + geom_point(data = rotated_ohc,aes(shift_width,Rank,color = rotated_ohc[,point_color],shape = rotated_ohc[,point_shape]),position = jitter_ohc,size = dot_size) +
        geom_point(data = rotated_ihc,aes(0,Rank,color = rotated_ihc[,point_color],shape = rotated_ihc[,point_shape]),position = jitter_ihc,size = dot_size) + 
        coord_fixed(ratio = plot_ratio) + theme_bw() + 
        scale_color_manual(values = c("#ec1c24","#2e3191")) +
        theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+
        facet_grid(~label)+
        theme(strip.text.x = element_text(size = 25)) +
        theme(legend.title = element_text(size = 30),legend.text = element_text(size = 25)) +
        labs(shape = "Cluster",color = "Library ID")
    }
  }
  
  if(point_color != "cluster" & point_color != "Experiment"){
    rotated_hc$label = point_color
    if(mid_value == "mean"){
      mid = mean(rotated_hc[,point_color],na.rm = T)
    }
    if(mid_value == "median"){
      mid = median(rotated_hc[,point_color],na.rm = T)
    }
    if(is.numeric(mid_value) == TRUE){
      mid = mid_value
    }
    
    if(plot_legend == FALSE){
      set.seed(seed)
      p1 = ggplot(rotated_hc) + geom_point(data = rotated_ohc,aes(shift_width,Rank,color = rotated_ohc[,point_color],shape = rotated_ohc[,point_shape]),position = jitter_ohc,size = dot_size) +
        geom_point(data = rotated_ihc,aes(0,Rank,color = rotated_ihc[,point_color],shape = rotated_ihc[,point_shape]),position = jitter_ihc,size = dot_size) + 
        coord_fixed(ratio = plot_ratio) + theme_bw() + 
        scale_color_gradient2(low = "white",high = "#E65948",mid = "#FEF6B1",midpoint = mid,na.value = "grey80") +
        theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+
        facet_grid(~label)+
        theme(strip.text.x = element_text(size = 25)) +
        theme(legend.title = element_text(size = 30),legend.text = element_text(size = 25),legend.position = "none") +
        labs(shape = "Cluster",color = point_color)
    }
    
    if(plot_legend == TRUE){
      set.seed(seed)
      p1 = ggplot(rotated_hc) + geom_point(data = rotated_ohc,aes(shift_width,Rank,color = rotated_ohc[,point_color],shape = rotated_ohc[,point_shape]),position = jitter_ohc,size = dot_size) +
        geom_point(data = rotated_ihc,aes(0,Rank,color = rotated_ihc[,point_color],shape = rotated_ihc[,point_shape]),position = jitter_ihc,size = dot_size) + 
        coord_fixed(ratio = plot_ratio) + theme_bw() + 
        scale_color_gradient2(low = "white",high = "#E65948",mid = "#FEF6B1",midpoint = mid,na.value = "grey80") +
        theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank())+
        facet_grid(~label)+
        theme(strip.text.x = element_text(size = 25)) +
        theme(legend.title = element_text(size = 30),legend.text = element_text(size = 25)) +
        labs(shape = "Cluster",color = point_color)
    }
  }
  
  return(p1)
}
