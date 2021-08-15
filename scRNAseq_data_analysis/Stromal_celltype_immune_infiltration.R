#Explore the relationship between stromal cell and immune cell infiltration
#First calculate the correatlion of cell type average gene expression and TCGA pation gene expression
#Then calcuate the Corratlion between CTL and correaltion of TCGA patient in the first step
#' @param Seurat         The Seurat object.
#' @param Cancer_exp     Path of the expression of cancer patient.
#' @param UpDEG_list     The upregulated gene list.
#' @param output_dir     The path of plot file store.

library(MAESRTO)
library(Seurat)
library(ggplot2)
Stromal_celltype_CTL<-function(Seurat,UpDEG_list,Cancer_exp,output_dir=NULL){
  
  cancer_expMat=readRDS(Cancer_exp)
  index=match(c("CD8A","CD8B","GZMA","GZMB","PRF1"),rownames(cancer_expMat))
  CTL_inf=colMeans(cancer_expMat[index,])
  
  expmat=GetAssayData(Seurat)
  expmat=as.matrix(expmat)
  cluster_DEG=UpDEG_list
  cluster_DEG=cluster_DEG[which(cluster_DEG$gene %in% rownames(expmat)),]
  if(length(output_dir) ==0){output_dir=getwd()}
  batch_do_ratio_survival<-function(x){
    cluster_DEG=cluster_DEG[order(cluster_DEG$avg_log2FC,decreasing = T),]
    cluster_gene=cluster_DEG$gene[which(cluster_DEG$cluster == x)][1:100]
    #ccell_index=which(MacrMono_seu@meta.data$curated_anno_gene == x)
    
    ccell_index=which(Seurat@meta.data$curated_anno == x)
    cgene_index=match(cluster_gene,rownames(expmat))
    cluster_avg=rowMeans(as.matrix(expmat[cgene_index,ccell_index]))
    background_avg=rowMeans(as.matrix(expmat[cgene_index,]))
    cluster_ratio=cluster_avg/background_avg
    
    
    
    tcga_gene_index=match(cluster_gene,rownames(cancer_expMat))
    temp_cancer_expMat=cancer_expMat[tcga_gene_index,]
    cor_result=apply(temp_cancer_expMat,2,function(y){unlist(cor.test(y,cluster_ratio))[c(4)]})
    data=data.frame(CTL=CTL_inf,Celltype_sig_cor=cor_result)
    
    p=ggplot(data,aes(x=Celltype_sig_cor,y=CTL))+geom_point(size=0.5,color="#3D7FBD")+theme_classic()
    ggsave(paste0(x,"_iCAF_2_dysfunction.pdf"),p,width = 3,height = 3)
    return(x,unlist(cor.test(data[,1],data[,2]))[4])
  }
  result=apply(matrix(unique(as.vector(unlist(cluster_DEG$cluster)))),1,batch_do_ratio_survival)
}
