#Calculating the signature score for each cell type, we could well describe the function of different cell types.
#In this function, all signarute gene are defined in your dataset.
#' @param Seurat_obj       The Seurat object.
#' @param gene_list        The cell type marker gene list
#' @param Celltype_order   The order of cell type
#' @param output_dir       The path of plot file store

library(Seurat)
library(ComplexHeatmap)
library(pheatmap)
library(circlize)

Celltype_annotation_heatmap<-function(Seurat_obj,gene_list,Celltype_order=NULL,output_dir=NULL){
  
  DefaultAssay(Seurat_obj)<-"RNA"
  expmat=GetAssayData(Seurat_obj,slot = "data",assay="RNA")
  #table(TNK_Seurat@meta.data$curated_anno)
  index=na.omit(match(gene_list,rownames(expmat)))
  Temp_exp=as.matrix(expmat[index ,])
  if(length(Celltype_order) ==0){
    Celltype_order=unique(as.vector(unlist(Seurat_obj@meta.data$curated_anno)))
  }
  
  
  Cell_type_marker<<-NULL
  get_celltype_AverExp<-function(x){
    cell_ind=match(rownames(Seurat_obj@meta.data)[which(Seurat_obj@meta.data$curated_anno == x)],colnames(expmat))
    celltype_aver=rowMeans(Temp_exp[,cell_ind])
    Cell_type_marker<<-cbind(Cell_type_marker,celltype_aver)
  }
  result=apply(as.matrix(Celltype_order),1,get_celltype_AverExp)
  colnames(Cell_type_marker)=cell_type
  
  if(length(output_dir) == 0){
    output_dir=getwd()
  }
  pdf("cell_type_annotation.pdf")
  col_fun = colorRamp2(c( 0, 3.5), c("white", "red"))
  Heatmap(Cell_type_marker,cluster_rows = FALSE,cluster_columns=FALSE,col = col_fun,column_split =c(rep("A",7) ,"B","B","C"))
  dev.off()
}