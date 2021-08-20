#' Using Cicero to find the co-accessibility in the specific region
#' Return Cicero object and gene annotation result in Rata file
#' If no gene are specifc which generate a expamle of NFKB2 gene
#' @param ATAC                  Seurat object of single cell ATAC
#' @param human_hg38_size       Path of human_hg38_size file 
#' @param Genome_annotation_gtf Path of Homo_sapiens.GRCh38.104.gtf file
#' @param TF_chr                The chromatin "chr1"
#' @param TF_star_site          The star site of plot
#' @param TF_end_site           The end site of plot
#' @param output_dir            Path of output file
#' @author Ya Han
#' 
library(Seurat)
library(ggplot2)
library(future)
library(RColorBrewer)
library(cicero)


Run_cicero<-function(ATAC,human_hg38_size,Genome_annotation_gtf,TF_chr="chr20",TF_star_site=51300000,TF_end_site=51600000,output_dir=NULL){

  DefaultAssay(Transfor_ATAC)<-"ATAC"
  peak_mat=GetAssayData(Transfor_ATAC)

  peak_info=rownames(peak_mat)
  peak_info=gsub("-","_",peak_info)
  peak_info=data.frame(site_name=peak_info)
  rownames(peak_info)=rownames(peak_mat)

  cell_metadata=data.frame(cells=colnames(peak_mat),Celltype=Transfor_ATAC@meta.data$assign.ident)
  #which(rownames(cell_metadata) != rownames(exhuasted_score))
  rownames(cell_metadata)=colnames(peak_mat)

  pd <- new("AnnotatedDataFrame", data = cell_metadata)
  fd <- new("AnnotatedDataFrame", data = peak_info)
  peak_mat=as.matrix(peak_mat)
  input_cds <-  suppressWarnings(newCellDataSet(peak_mat,
                                              phenoData = pd,
                                              featureData = fd,
                                              expressionFamily=VGAM::binomialff(),
                                              lowerDetectionLimit=0))

  input_cds@expressionFamily@vfamily <- "binomialff"
  input_cds <- monocle::detectGenes(input_cds)

  #Ensure there are no peaks included with zero reads
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

  set.seed(2017)
  input_cds <- detectGenes(input_cds)
  input_cds <- estimateSizeFactors(input_cds)

# *** if you are using Monocle 3 alpha, you need to run the following line as well!
# NOTE: Cicero does not yet support the Monocle 3 beta release (monocle3 package). We hope
# to update soon!
#input_cds <- preprocessCDS(input_cds, norm_method = "none")
  input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                             reduction_method = 'tSNE', norm_method = "none")
  tsne_coords <- t(reducedDimA(input_cds))
  row.names(tsne_coords) <- row.names(pData(input_cds))
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)


  human_hg38_size=read.table(human_hg38_size,header = FALSE,sep="\t")
  #human_hg38_size=read.table("/fs/home/hanya/Reference/hg38.chrom.sizes",header = FALSE,sep="\t")
#sample_genome <- subset(human_hg38_size, V1 == "chr20")
  conns <- run_cicero(cicero_cds, human_hg38_size) # Takes a few minutes to run

  setwd(output_dir)
  gene_anno <- rtracklayer::readGFF(Genome_annotation_gtf)
# rename some columns to match plotting requirements
  gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
  gene_anno$gene <- gene_anno$gene_id
  gene_anno$transcript <- gene_anno$transcript_id
  gene_anno$symbol <- gene_anno$gene_name

  if(length(output_dir) == 0 ){
    output_dir=getwd()
  }
  setwd(output_dir)
  
  pdf("Rplot.pdf",height = 3,width = 6)
  gene_annotation_sample=rbind(gene_annotation_sample,gene_annotation_NFATC2)
  plot_connections(conns, TF_chr, TF_star_site, TF_end_site, 
                 gene_model = gene_anno, 
                 coaccess_cutoff = .25, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )
  dev.off()
  saveRDS(conns,gene_anno,"Cicero_object_gene_annotation.RData")
}
