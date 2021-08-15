#Running Monocle to infer the pseudotime.
##Creating an monocle object based on your three required inputs: the read count matrix, cell meta information, and the gene information:
#' @param Seurat_obj The input is a seurat obj
#' @param mean_cutoff The cut-off of average gene expression and used to select highly expressed gene to run monocle3
#' @param output_dir The path of plot file store

library(monocle)
library(Seurat)


Run_Monocle<-function(Seurat_obj,mean_cutoff=0.0125,output_dir=NULL){
  expmat=GetAssayData(Seurat_obj)
  cell_metadata=Seurat_obj@meta.data
  expmat=as.matrix(expmat)
  index=which(rowMeans(expmat) > mean_cutoff)
  expmat=expmat[index,]
  
  gene_metadata=rownames(expmat)
  gene_metadata=as.data.frame(gene_metadata)
  rownames(gene_metadata)=rownames(expmat)
  colnames(gene_metadata)="gene_short_name"
  
  colnames(expmat)=rownames(cell_metadata)
  
  pd <- new("AnnotatedDataFrame", data = cell_metadata)
  fd <- new("AnnotatedDataFrame", data = gene_metadata)
  
  #cds <- newCellDataSet(expmat, phenoData = pd, featureData = fd)
  colnames(expmat)=rownames(cell_metadata)
  
  cds <- newCellDataSet(as(expmat, "sparseMatrix"),
                        phenoData = pd,
                        featureData = fd)
  
  
  
  cds <- estimateSizeFactors(cds)
  cds=estimateDispersions(cds)
  #cds <- setOrderingFilter(cds, ordering_genes)
  cds <- reduceDimension(cds)
  cds <- orderCells(cds)
  
  
  p=plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size=0.5,show_tree = F,show_backbone=F,show_branch_points=F)
  ggsave("Monocle_Differention_Traj_pseudotime.pdf",p,width = 5,height = 4)
  
  retrun(cds)
}