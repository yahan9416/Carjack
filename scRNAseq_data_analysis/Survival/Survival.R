#Look at the survival effect of cell type upregulated gene on TCGA HNSC patients.
#There are three ways to look at the effect: 
#First, for each patient calculate the GSVA score of the cell type upregulated gene. 
#Second, calculate the correlation between the average expression of the upregulated gene in the cell type and HNSC patient expression level. 
#Third, calculate the correlation between the ratio of the cell type up-regulated gene expression compares with all cell types and HNSC patient expression level.

#' @param Seurat_obj    The input is a seurat obj, the meta information of seurat need including curated_anno.
#' @param upDEG_list    The cell-type/cluster up-regulated gene list
#' @param Cancer_exp    The expression of TCGA specific cancer type patient expression profile
#' @param Survival_info The file of clinical information
#' @param output_dir    The path of plot file store


library(Seurat)
library(survminer)
library(survival)
library(ggplot2)
library(GSVA)

#Third, calculate the correlation between the ratio of the cell type up-regulated gene expression compares with all cell types and HNSC patient expression level.
Ratio_correlation_survival<-function(Seurat_obj,upDEG_list=NULL,Cancer_exp,Survival_info,output_dir=NULL){
  
  cancer_expMat=readRDS(Cancer_exp)
  Survival_info=t(read.table(Survival_info,header = T,sep=",",row.names = 1))
  Survival_info<<-as.data.frame(Survival_info)
  time_index=match("_OS",colnames(Survival_info))
  statu_index=match("_OS_IND",colnames(Survival_info))
  
  if(length(upDEG_list) == 0){
    upDEG_list=FindAllMarkers(Seurat_obj,only.pos = TRUE)
  }
  cluster_DEG=upDEG_list
  cluster_DEG=cluster_DEG[which(cluster_DEG$gene %in% rownames(cancer_expMat)),]
  if(length(output_dir) == 0){
    output_dir=getwd()
  }
  
  survival_infor<<-NULL
  expmat=GetAssayData(Seurat_obj)
  batch_do_ratio_survival<-function(x){
    cluster_gene=cluster_DEG$gene[which(cluster_DEG$cluster == x)]
    ccell_index=which(Seurat_obj@meta.data$curated_anno == x)
    cgene_index=match(cluster_gene,rownames(expmat))
    cluster_avg=rowMeans(as.matrix(expmat[cgene_index,ccell_index]))
    background_avg=rowMeans(as.matrix(expmat[cgene_index,]))
    cluster_ratio=cluster_avg/background_avg
    
    survival_infor=data.frame(time=Survival_info[,time_index],statu=Survival_info[,statu_index])
    rownames(survival_infor)=gsub("\\.","-",rownames(Survival_info))
    tcga_gene_index=match(cluster_gene,rownames(cancer_expMat))
    temp_cancer_expMat=cancer_expMat[tcga_gene_index,]
    cor_result=apply(temp_cancer_expMat,2,function(y){unlist(cor.test(y,cluster_ratio))[c(4)]})
    sample_index=match(rownames(survival_infor),names(cor_result))
    cor_result=as.numeric(as.vector(unlist(cor_result)))
    survival_infor$Ratio_Corr=cor_result[sample_index]
    
    
    survival_infor$Cor_class=1
    survival_infor$Ratio_Corr=as.numeric(as.vector(survival_infor[,3]))
    survival_infor$Cor_class[which(survival_infor$Ratio_Corr < median(survival_infor$Ratio_Corr) )]=0
    #survival_infor=survival_infor[which(survival_infor$time < 365*5),]           
    
    fit<-survfit(Surv(time,statu)~Cor_class, data=survival_infor)
    surv_diff <- survdiff(Surv(time,statu)~Cor_class, data=survival_infor)
    p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
    #p.value
    
    setwd(output_dir)
    ggsurv <- ggsurvplot(fit,palette = "jco", pval = FALSE,risk.table = T)+ggtitle(paste0(x,"_survival_eff:",format(p.value,digits = 3)))
    p=ggarrange(ggsurv$plot, ggsurv$table, heights = c(1.5, 0.8),ncol = 1, nrow = 2)
    ggsave(paste0(x,"_ratio_correlation_TCGA_survival.pdf"),p,width = 7,height = 6)
    
    
  }
  result=apply(matrix(unique(cluster_DEG$cluster)),1,batch_do_ratio_survival)
}

#Second, calculate the correlation between the average expression of the upregulated gene in the cell type and HNSC patient expression level. 
Correlation_survival<-function(Seurat_obj,upDEG_list=NULL,Cancer_exp,Survival_info,output_dir=NULL){
  
  cancer_expMat=readRDS(Cancer_exp)
  Survival_info=t(read.table(Survival_info,header = T,sep=",",row.names = 1))
  Survival_info<<-as.data.frame(Survival_info)
  time_index=match("_OS",colnames(Survival_info))
  statu_index=match("_OS_IND",colnames(Survival_info))
  
  if(length(upDEG_list) == 0){
    upDEG_list=FindAllMarkers(Seurat_obj,only.pos = TRUE)
  }
  cluster_DEG=upDEG_list
  cluster_DEG=cluster_DEG[which(cluster_DEG$gene %in% rownames(cancer_expMat)),]
  if(length(output_dir) == 0){
    output_dir=getwd()
  }
  
  survival_infor<<-NULL
  expmat=GetAssayData(Seurat_obj)
  batch_do_ratio_survival<-function(x){
    cluster_gene=cluster_DEG$gene[which(cluster_DEG$cluster == x)]
    ccell_index=which(Seurat_obj@meta.data$curated_anno == x)
    cgene_index=match(cluster_gene,rownames(expmat))
    cluster_avg=rowMeans(as.matrix(expmat[cgene_index,ccell_index]))
    
    
    survival_infor=data.frame(time=Survival_info[,time_index],statu=Survival_info[,statu_index])
    rownames(survival_infor)=gsub("\\.","-",rownames(Survival_info))
    tcga_gene_index=match(cluster_gene,rownames(cancer_expMat))
    temp_cancer_expMat=cancer_expMat[tcga_gene_index,]
    cor_result=apply(temp_cancer_expMat,2,function(y){unlist(cor.test(y,cluster_avg))[c(4)]})
    sample_index=match(rownames(survival_infor),names(cor_result))
    cor_result=as.numeric(as.vector(unlist(cor_result)))
    survival_infor$Corr=cor_result[sample_index]
    
    
    survival_infor$Cor_class=1
    survival_infor$Corr=as.numeric(as.vector(survival_infor[,3]))
    survival_infor$Cor_class[which(survival_infor$Corr < median(survival_infor$Corr) )]=0
    #survival_infor=survival_infor[which(survival_infor$time < 365*5),]           
    
    fit<-survfit(Surv(time,statu)~Cor_class, data=survival_infor)
    surv_diff <- survdiff(Surv(time,statu)~Cor_class, data=survival_infor)
    p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
    #p.value
    
    setwd(output_dir)
    ggsurv <- ggsurvplot(fit,palette = "jco", pval = FALSE,risk.table = T)+ggtitle(paste0(x,"_survival_eff:",format(p.value,digits = 3)))
    
    p=ggarrange(ggsurv$plot, ggsurv$table, heights = c(1.5, 0.8),ncol = 1, nrow = 2)
    ggsave(paste0(x,"_ratio_correlation_TCGA_survival.pdf"),p,width = 7,height = 6)
    
    
  }
  result=apply(matrix(unique(cluster_DEG$cluster)),1,batch_do_ratio_survival)
}


#First, for each patient calculate the GSVA score of the cell type upregulated gene. 
DEG_GSVA_survival<-function(upDEG_list=NULL,Cancer_exp,Survival_info,output_dir=NULL){
  cancer_expMat=readRDS(Cancer_exp)
  Survival_info=t(read.table(Survival_info,header = T,sep=",",row.names = 1))
  Survival_info<<-as.data.frame(Survival_info)
  time_index=match("_OS",colnames(Survival_info))
  statu_index=match("_OS_IND",colnames(Survival_info))
  
  
  if(length(output_dir) == 0){
    output_dir=getwd()}
  
  batch_cluster_survival<-function(x){
    gsva_score=gsva(cancer_expMat,list(cluster_DEG$gene[which(cluster_DEG$cluster == x)]))
    gsva_score=t(gsva_score)
    colnames(gsva_score)=x
    survival_infor=data.frame(time=Survival_info[,time_index],statu=Survival_info[,statu_index])
    rownames(survival_infor)=gsub("\\.","-",rownames(Survival_info))
    index=match(rownames(survival_infor),rownames(gsva_score))
    survival_gsva<<-cbind(survival_infor,gsva_score[index])
    
    survival_gsva=cbind(survival_gsva,Cor_class=0)
    survival_gsva$Cor_class=1
    survival_gsva$gene_cor_ration=as.numeric(as.vector(survival_gsva[,5]))
    survival_gsva$Cor_class[which(survival_gsva$gene_cor_ration < median(survival_gsva$gene_cor_ration) )]=0
    
    fit<-survfit(Surv(time,statu)~Cor_class, data=survival_gsva)
    surv_diff <- survdiff(Surv(time,statu)~Cor_class, data=survival_gsva)
    p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
    setwd(output_dir)
    ggsurv <- ggsurvplot(fit,palette = "jco", pval = FALSE,risk.table = T)+ggtitle(paste0(x,"_survival_eff:",format(p.value,digits = 3)))
    p=ggarrange(ggsurv$plot, ggsurv$table, heights = c(1.5, 0.8),ncol = 1, nrow = 2)
    ggsave(paste0(x,"_gsva_TCGA_HNSC_survival.pdf"),p,width = 7,height = 6)
    
  }
  result=apply(matrix(unique(cluster_DEG$cluster)),1,batch_cluster_survival)
}
