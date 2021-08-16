#Running Monocle3 to infer the pseudotime.
##Creating an monocle3 object based on your three required inputs: the read count matrix, cell meta information, and the gene information:
#' @param Seurat_obj The input is a seurat obj
#' @param mean_cutoff The cut-off of average gene expression and used to select highly expressed gene to run monocle3
#' @param output_dir The path of plot file store

library(monocle3)
library(Seurat)


Run_Monocle3<-function(Seurat_obj,mean_cutoff=0.0125,output_dir=NULL){
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
  cds <- new_cell_data_set(expmat,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_metadata)
  
  ## Step 1: Normalize and pre-process the data
  cds <- preprocess_cds(cds, num_dim = 100)
  ## Step 3: Reduce the dimensions using UMAP
  cds <- reduce_dimension(cds,preprocess_method="PCA")
  ## Step 4: Cluster the cells
  cds <- cluster_cells(cds)
  ## Step 5: Learn a graph
  cds <- learn_graph(cds)
  
  ## Step 6: Order cells
  get_earliest_principal_node <- function(cds, time_bin="Naïve CD8 T"){
    cell_ids <- which(colData(cds)[, "curated_anno"] == time_bin)
    
    closest_vertex <-
      cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
      igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                (which.max(table(closest_vertex[cell_ids,]))))]
    
    root_pr_nodes
  }
  cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
  
  if(length(output_dir) == 0){
    output_dir=getwd()
  }
  setwd(output_dir)
  
  
  p=plot_cells(cds,
               color_cells_by = "pseudotime",
               label_cell_groups=FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE,
               graph_label_size=1.5)
  
  ggsave("Monocle3_Differention_Traj_pseudotime.pdf",p,width = 5,height = 4)
  
  
  retrun(cds)
}
