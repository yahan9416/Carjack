#' Integrate scATACseq and scRNAseq data by gene expression and RP matrix 
#' @param RNA Seurat object of single cell RNAseq
#' @param ATAC Seurat object of single cell ATAC
#' @param nfeatures How many feature will be used for find anchors and appeared in the integrated object.
#' @param project Name of project 
#' @param output_dir Path of output file
#' @author Ya Han

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)

Integrate_scRNA_scATAC<-function(RNA,ATAC,nfeatures=3000,project="Integrate_scRNA_scATAC",output_dir){
   ATAC$tech <- "ATAC"
   RNA$tech <- "RNA"

   DefaultAssay(ATAC) <- "ACTIVITY"
   RNA <- FindVariableFeatures(RNA, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
   ATAC <- FindVariableFeatures(ATAC, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
   
   transfer.anchors <- FindTransferAnchors(reference = RNA, query = ATAC, features = VariableFeatures(object = RNA), 
                                           reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
  
   setwd(output_dir)
   #transfer.anchors=readRDS("/fs/home/hanya/Project/Carcinogenesis/scATAC/Integrate_ATAC_RNA_together/Media_file/ATAC_RNA_allcell_transfer.anchors.rds")
   RNA@meta.data$assign.ident=RNA@meta.data$assign.ident
   celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = RNA$assign.ident, weight.reduction = ATAC[["lsi"]], dims = dims.use)
   ATAC@meta.data$assign.ident <- celltype.predictions$predicted.id
   ATAC@meta.data$prediction.score.max <- celltype.predictions$prediction.score.max
   saveRDS(ATAC@meta.data,file=paste0("RNAATAC_Integr_subcelltype_meta.rds"))
   saveRDS(celltype.predictions,file="ATAC_Predict_celltype_transfor.rds")
}
