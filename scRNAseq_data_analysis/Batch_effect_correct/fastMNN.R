#' Correct batch effect by fastMNN and return a seurat object.
#' @param H5_files Path of expression matrix of h5 file
#' @param Seurat_obj If you have built a Seurat object, and you can use it as input
#' @param Batch_info A vector of batch information corresponding to cell barcode in expmat or seurat object
#' @param nfeatures How many feature will be used for find anchors and appeared in the integrated object.
#' @author Ya Han

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(batchelor)

Run_fastMNN<-function(H5_files=NULL,Seurat_obj=NULL,Batch_info,nfeatures=3000){
  if(length(H5_files) >0){
    # read data
    expr = Read10X_h5(H5_files)
    seuratObj <- CreateSeuratObject(counts = expr, project = "Seurat_obj", min.cells = 10, min.features = 500)
    seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
    seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
    seuratObj <- ScaleData(seuratObj, features = rownames(seuratObj$RNA))
    seuratObj <- RunPCA(seuratObj, features = VariableFeatures(object = seuratObj),npcs = 100)
    
    Seurat_obj=seuratObj
  }
  Seurat_obj@meta.data$Batch=Batch_info
  
  expMat=GetAssayData(Seurat_obj,assay="RNA")
  Seurat_obj <- FindVariableFeatures(Seurat_obj, selection.method = "vst",nfeatures=nfeatures)
  varible_gene=VariableFeatures(Seurat_obj)
  Corr_expmat <- fastMNN(expMat, batch=Batch_info,subset.row=varible_gene)
  seuratObj <- CreateSeuratObject(counts = Corr_expmat, project = "Seurat_obj_fastMNN", min.cells = 10, min.features = 500)
  seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
  seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
  seuratObj <- ScaleData(seuratObj, features = rownames(seuratObj$RNA))
  seuratObj <- RunPCA(seuratObj, features = VariableFeatures(object = seuratObj),npcs = 100)
  
  return(seuratObj)
  
}
