#Plot the signal of single cell level by Signac
#' @param ATAC                 Seurat object of single cell ATAC sequencing.
#' @param fragments_size       The path of fragments_size file Exhuasted_T_cells_fragment_size_sort.tsv.gz.
#' @param Genome_region        The region of genome to plot "chr20-51300000-51600000"
#' @param output_dir           The path of plot file store.


library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)




Signac_vis_signal<-function(ATAC,fragments_size,Genome_region="chr20-51300000-51600000",output_dir=NULL){
  DefaultAssay(ATAC)<-'ATAC'
  peak_count=GetAssayData(ATAC)
  metadata=ATAC@meta.data


  chrom_assay <- CreateChromatinAssay(
    counts = peak_count,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = fragments_size,
    min.cells = 10,
    min.features = 20,
    validate.fragments=FALSE)

  ATAC_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata)

  #ATAC_obj@meta.data$cluster=metadata$cluster[match(rownames(ATAC_obj@meta.data),metadata[,1])]
  #Idents(ATAC_obj)=ATAC_obj@meta.data$cluster

  p_sing=tile_plot <- TilePlot(
    object = ATAC_obj,
    region = Genome_region)

  if(length(output_dir) == 0){
    output_dir=getwd()
  }
  setwd(output_dir)
  ggsave("Cicero_example.pdf",p_sing)
  saveRDS(ATAC_obj,file="ATAC_new_version.rds")
}