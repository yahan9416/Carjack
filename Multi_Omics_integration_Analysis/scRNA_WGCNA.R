#' In order to explore the association between genes, construct a gene co-expression network.
#' @param Seurat_obj Seurat object of single cell RNAseq
#' @param nfeatures How many feature will be used for analysis.
#' @param output_dir Path of output file
#' @author Ya Han

library(WGCNA)
library(Seurat)
library(tidyverse)
library(reshape2)
library(stringr)

scWGCNA<-function(Seurat_obj,nfeatures=3000,output_dir){
  datadf <- Seurat_obj@assays$RNA@data
  DefaultAssay(Seurat_obj)<-"RNA"
  Seurat_obj=FindVariableFeatures(Seurat_obj,nfeatures=nfeatures)
  HV_gene=VariableFeatures(Seurat_obj)

  Cell_minicluster_list<<-list()
  colname_source_cluster<<-NULL
  rownames(Seurat_obj@meta.data)=names(Idents(Seurat_obj))
  Seurat_obj@meta.data$sample_cluster=paste(Seurat_obj@meta.data$sample,Seurat_obj@meta.data$seurat_clusters,sep = "|")
  DefaultAssay(Seurat_obj)<-"integrated"
  
  cluster_sample_average<-function(cluster){
    
    sample_index=rownames(Seurat_obj@meta.data)[which( Seurat_obj@meta.data$sample_cluster ==cluster)]
    
    if(length(sample_index) > 20){
      temp_seurat=subset(Seurat_obj,cells=sample_index)
      
      cells <- length(sample_index)
      cluster.res=0.2
      dims.use=1:15
      temp_seurat <- FindNeighbors(temp_seurat, dims = dims.use)
      temp_seurat=FindClusters(temp_seurat,resolution = cluster.res)
      #print(cluster)
      #print(table(temp_seurat@meta.data$seurat_clusters))
      temp_seurat@meta.data$source_cluster=temp_seurat@meta.data$seurat_clusters
      
      Find_the_near_100Neig_incluster<-function(y){
        print(y)
        umap_df <<- temp_seurat@reductions$umap@cell.embeddings[which(temp_seurat@meta.data$source_cluster == y),]
        cell_dist<<-as.matrix(dist(umap_df))
        if(nrow(cell_dist) > 20){
          #if more than 50 cell will generate a new cluster or will not
          if(ceiling(length(which(temp_seurat@meta.data$source_cluster == y))/20) == round(length(which(temp_seurat@meta.data$source_cluster == y))/20) ){
            temp_cluster<<-ceiling(length(which(temp_seurat@meta.data$source_cluster == y))/20)
            #first generate the first cluster and remove less than 100 cells from the total, because there less than temp_cluster*100, so there must with overlap between clusters
            seed_cell=sample(rownames(cell_dist),1)
            cluster_cells=names(sort(cell_dist[match(seed_cell,rownames(cell_dist)),],decreasing=F)[1:20])
            Cell_minicluster_list<<-c(Cell_minicluster_list,list(cluster_cells))
            
            colname_source_cluster<<-c(colname_source_cluster,paste(cluster,y,sep="|"))
            less10_cell_number=(length(which(temp_seurat@meta.data$source_cluster == y)) - floor(length(which(temp_seurat@meta.data$source_cluster == y))/20)*20)
            
            
            index=match(cluster_cells,rownames(cell_dist))
            cell_dist<<-cell_dist[-1*index,-1*index]
            
            generate_50cell_expresion<-function(sub_cluster){
              
              if(sub_cluster != temp_cluster){
                seed_cell=sample(rownames(cell_dist),1)
                cluster_cells=names(sort(cell_dist[match(seed_cell,rownames(cell_dist)),],decreasing=F)[1:20])
                Cell_minicluster_list<<-c(Cell_minicluster_list,list(cluster_cells))
                index=match(cluster_cells,rownames(cell_dist))
                cell_dist<<-cell_dist[-1*index,-1*index]
              }else{
                Cell_minicluster_list<<-c(Cell_minicluster_list,list(rownames(cell_dist)))
              }
              colname_source_cluster<<-c(colname_source_cluster,paste(cluster,y,sub_cluster,sep="|"))
            } 
            result=apply(as.matrix(2:temp_cluster),1,generate_50cell_expresion)
            
          }else{
            temp_cluster<<-round(length(which(temp_seurat@meta.data$source_cluster == y))/15)
            generate_50cell_expresion<-function(sub_cluster){
              seed_cell=sample(rownames(cell_dist),1)
              if(sub_cluster==temp_cluster){
                cluster_cells=colnames(cell_dist)
                Cell_minicluster_list<<-c(Cell_minicluster_list,list(cluster_cells))
                colname_source_cluster<<-c(colname_source_cluster,paste(cluster,y,sub_cluster,sep="|"))
              }else{
                cluster_cells=names(sort(cell_dist[match(seed_cell,rownames(cell_dist)),],decreasing=F)[1:15])
                cell_index=match(cluster_cells,colnames(cell_dist))
                cell_dist<<-cell_dist[-1*cell_index,-1*cell_index]
                Cell_minicluster_list<<-c(Cell_minicluster_list,list(cluster_cells))
                colname_source_cluster<<-c(colname_source_cluster,paste(cluster,y,sub_cluster,sep="|"))
              }
              
            } 
            result=apply(as.matrix(1:temp_cluster),1,generate_50cell_expresion)
          }
          
        }else{
          #cluster including total less than 100 cells 
          
          Cell_minicluster_list<<-c(Cell_minicluster_list,list(colnames(cell_dist)))
          colname_source_cluster<<-c(colname_source_cluster,paste(cluster,y,sep="|"))
        }
        
      } 
      #source_cluster_moerthan100=names(which(table(temp_seurat@meta.data$source_cluster) >= 100))
      result=apply(as.matrix(unique(temp_seurat@meta.data$seurat_clusters)),1,Find_the_near_100Neig_incluster)
      
    }else{
      
      Cell_minicluster_list<<-c(Cell_minicluster_list,list(sample_index))
      
      colname_source_cluster<<-c(colname_source_cluster,paste(cluster,0,sep="|"))
      
    }
    
  }
  result=apply(as.matrix(unique(as.vector(unlist( Seurat_obj@meta.data$sample_cluster)))),1,cluster_sample_average)    
  
  if(length(output_dir) == 0){
    output_dir=getwd()
  }
  setwd(output_dir)
    save(Cell_minicluster_list,colname_source_cluster,file="SeuratObj_mini_cluster_20_15_cells_index.RData")

  result=unlist(lapply(Cell_minicluster_list,function(x) return(length(x))))
  less_than_5_index=which(result < 5)
  
  hv_gene_index=match(HV_gene,rownames(datadf))
  datadf=datadf[hv_gene_index,]
  datadf=as.matrix(datadf)
  Expmat_profile<<-NULL
  generate_expression_profile<-function(x){
    sam_ind=match(Cell_minicluster_list[[x]],colnames(datadf))
    Expmat_profile<<-cbind(Expmat_profile,rowMeans(datadf[,sam_ind]))
  }
  result=lapply(matrix(setdiff(1:length(Cell_minicluster_list),less_than_5_index)),generate_expression_profile)
  
  #Build WGCNA NETWORK
  {
    type = "unsigned" 
    corType = "pearson" 
    corFnc = ifelse(corType=="pearson", cor, bicor)
    corFnc
    maxPOutliers = ifelse(corType=="pearson",1,0.05) 
    robustY = ifelse(corType=="pearson",T,F)
    dataExpr  <- as.matrix(Expmat_profile)
    
    m.mad <- apply(dataExpr,1,mad)
    dataExprVar <- dataExpr[which(m.mad > 
                                    max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
    

    dataExpr <- as.data.frame(t(dataExprVar))
    
    gsg = goodSamplesGenes(dataExpr, verbose = 3)
    gsg$allOK
    gsg$goodSamples
    
    if (!gsg$allOK){
      # Optionally, print the gene and sample names that were removed:
      if (sum(!gsg$goodGenes)>0) 
        printFlush(paste("Removing genes:", 
                         paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
      if (sum(!gsg$goodSamples)>0) 
        printFlush(paste("Removing samples:", 
                         paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
      # Remove the offending genes and samples from the data:
      dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
    }
    
    nGenes = ncol(dataExpr)
    nSamples = nrow(dataExpr)
    
    sampleTree = hclust(dist(dataExpr), method = "average")
    plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
    powers = c(c(1:10), seq(from = 12, to=30, by=2))
    sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                            networkType="signed", verbose=5)
    
    par(mfrow = c(1,2))
    cex1 = 0.9
    # 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
    # 网络越符合无标度特征 (non-scale)
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",
         ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red")
    # 筛选标准。R-square=0.85
    abline(h=0.85,col="red")
    
    # Soft threshold与平均连通性
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
         cex=cex1, col="red")
    
    power = sft$powerEstimate
    softPower  = power
    softPower
    
    if (is.na(power)){
      power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                     ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                            ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                   ifelse(type == "unsigned", 6, 12))       
                     )
      )
    }
    
    cor <- WGCNA::cor
    
    net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,#nGenes
                           TOMType = "unsigned", minModuleSize = 10,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs=TRUE, corType = corType, 
                           maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                           saveTOMFileBase = paste0("dataExpr", ".tom"),
                           verbose = 3)
    
    moduleLabels = net$colors
    moduleColors = labels2colors(moduleLabels)
    moduleColors
    # Plot the dendrogram and the module colors underneath
    # 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
    pdf("plotDendroAndColors.pdf")
    plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    dev.off()
    dev.off()
    
    MEs = net$MEs

    MEs_col = MEs
    colnames(MEs_col) = paste0("ME", labels2colors(
      as.numeric(str_replace_all(colnames(MEs),"ME",""))))
    MEs_col = orderMEs(MEs_col)
    
    pdf("Correlation_between_module.pdf")
    plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                          marDendro = c(3,3,2,4),
                          marHeatmap = c(3,4,2,2),
                          plotDendrograms = T,
                          xLabelsAngle = 90)
    dev.off()
    
    dynamicColors=moduleColors
    MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
    #Module value of each sample
    MEs = MEList$eigengenes  
    saveRDS(MEs,"Seuratobj_WGCNA_Module_score.rds")
    saveRDS(net,file="Seuratobj_WGCNA_Module_network.rds")
  }
}