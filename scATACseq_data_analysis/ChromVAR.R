#'Using ChromVAR to identify the drive TF of each cells 
#' @author Ya Han

library(BiocParallel)#BiocParallel to specify your preferred method of parallelization
register(MulticoreParam(10))
library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(MAESTRO)
library(Seurat)
library(ggplot2)


#Example
Peak_count=Read10X_h5("/mnt/Storage2/home/hanya/project/Carcinogenesis/scATAC/Maestro_process_ATAC/HNSC_Oral_scATAC_filtered_peak_count.h5")

#the peak count matrix min value should >= 0
#But after merge peak and calculate the peak count, there are 0 is the merge peak count
Chromvar_TF=RunchromVAR(Peak_count)

RunchromVAR <- function(inputMat, min.c = 50, min.p = 500, organism = "GRCh38")
{
  message("Reading File!")
  #   spam_input_matrix <- as.spam.dgCMatrix(inputMat)
  peaks <- data.frame(chr=unlist(strsplit(rownames(inputMat),'_'))[seq(1,nrow(inputMat)*3,3)], 
                      start=unlist(strsplit(rownames(inputMat),'_'))[seq(2,nrow(inputMat)*3,3)],
                      end=unlist(strsplit(rownames(inputMat),'_'))[seq(3,nrow(inputMat)*3,3)])  
  rownames(peaks) <- paste0(peaks[,1],'_',peaks[,2],'_',peaks[,3])
  peaks <- na.omit(peaks)
  message("Generating col.data!")
  #   col.Data <- DataFrame(celltype=colnames(inputMat), depth=apply.spam(spam_input_matrix,2,sum))
  col.Data <- DataFrame(celltype=colnames(inputMat), depth=Matrix::colSums(inputMat))
  col.Data <- col.Data[which(col.Data[,2]>min.p),]
  message("Generating row.data!")
  #   row.Data <- DataFrame(peaks=rownames(inputMat), depth=apply.spam(spam_input_matrix[,col.Data[,1]],1,sum))
  row.Data <- DataFrame(peaks=rownames(inputMat), depth=Matrix::rowSums(inputMat[,col.Data[,1]]))
  row.Data <- row.Data[which(row.Data[,2]>min.c),]
  #   finalPeak <- as.matrix(peaks[row.Data[,1],])
  #   finalCount <- as.matrix(inputMat[row.Data[,1],col.Data[,1]])
  finalPeak <- peaks[row.Data[,1],]
  finalCount <- inputMat[row.Data[,1],col.Data[,1]]
  rowRanges <- GRanges(as.character(finalPeak[,1]),IRanges(as.integer(finalPeak[,2]),as.integer(finalPeak[,3])),
                       strand="*",score=as.integer(5),qval=1,name=rownames(finalPeak))
  resizedRowRanges <- resize(rowRanges, width=500, fix="center")
  chromVAR_counts  = SummarizedExperiment(assays=SimpleList(counts=finalCount),rowRanges=resizedRowRanges, colData=col.Data)
  
  
  message("chromVAR analysis ...")
  if(organism == "GRCh38") {
    chromVAR_counts = addGCBias(chromVAR_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
    #jaspar的motif少一点
    motifs <- getJasparMotifs(species = "Homo sapiens")
    motif_ix = matchMotifs(motifs, chromVAR_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
    #data(human_pwms_v2) 是chromvar 里面的 database
    #motif_ix = matchMotifs(human_pwms_v2, chromVAR_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
  }
  if(organism == "GRCm38") {
    chromVAR_counts = addGCBias(chromVAR_counts, genome = BSgenome.Mmusculus.UCSC.mm10)
    #      motifs <- getJasparMotifs(species="Mus musculus")
    #      motif_ix = matchMotifs(motifs, chromVAR_counts, genome = BSgenome.Mmusculus.UCSC.mm10)
    data(mouse_pwms_v2)
    motif_ix = matchMotifs(mouse_pwms_v2, chromVAR_counts, genome = BSgenome.Mmusculus.UCSC.mm10)}  
  dev = computeDeviations(object = chromVAR_counts, annotations = motif_ix)
  variability = computeVariability(dev)
  zscore = assays(dev)$z
  rownames(zscore) = variability[,1]
  #   variability <- computeVariability(dev)
  
  message("Finished!")
  return(list(dev = dev, zscore = zscore, variability = variability))
}