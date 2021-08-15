#Calculating the signature score for each cell type, we could well describe the function of different cell types.
#In this function, all signarute gene are defined in your dataset.
#' @param T_seurat       The Seurat object only including CD4 and CD8 T cells.
#' @param upDEG_list     Up-regulated gene list
#' @param output_dir     The path of plot file store

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(ggpubr)

T_signature_score<-function(T_seurat,upDEG_list=NULL,output_dir=NULL){
  expmat=GetAssayData(T_seurat)
  if(length(upDEG_list) == 0){
    upDEG_list=FindAllMarkers(T_seurat,only.pos = TRUE)
  } 
  DEGs=upDEG_list
  #Find the cell type belong to CD8T cells
  
  DEGs=DEGs[grep("CD8",DEGs$cluster),]
  meta.data=T_seurat@meta.data
  
  #The first condition is DEGs and second is in CD8 Tcells
  index=grep("CD8",meta.data$curated_anno)
  expmat=expmat[,index]
  
  #cytotoxicity associated gene selected
  index=match("GZMK",rownames(expmat))
  cal_cor<-function(x){
    index2=match(x,rownames(expmat))
    cor=unlist(cor.test(as.vector(expmat[index2,]),as.vector(expmat[index,])))[3:4]
    return(c(x,cor))
  }
  result=apply(matrix(unique(DEGs$gene)),1,cal_cor)
  result=as.data.frame(t(result))
  result$estimate.cor=as.numeric(as.vector(result$estimate.cor))
  result=result[order(result$estimate.cor,decreasing = T),]
  result=result[which(!duplicated(result)),]
  cytotoxicity_fea=result[1:20,1]
  
  
  
  cd8_cell=rownames(meta.data)[grep("CD8",meta.data$curated_anno)]
  TNK_cd8t<-subset(x = T_seurat,cells=cd8_cell,invert = FALSE)
  rownames(TNK_cd8t@meta.data)=names(Idents(TNK_cd8t))
  TNK_cd8t <- AddModuleScore(
    object = TNK_cd8t,
    features = list(cytotoxicity_fea),
    assay	="RNA",
    ctrl = 5,
    name = 'Cytotoxic_GZMK'
  )
  
  index=match("HAVCR2",rownames(expmat))
  cal_cor<-function(x){
    index2=match(x,rownames(expmat))
    cor=unlist(cor.test(as.vector(expmat[index2,]),as.vector(expmat[index,])))[3:4]
    return(c(x,cor))
  }
  result1=apply(matrix(unique(DEGs$gene)),1,cal_cor)
  
  result1=as.data.frame(t(result1))
  result1$estimate.cor=as.numeric(as.vector(result1$estimate.cor))
  result1=result1[order(result1$estimate.cor,decreasing = T),]
  result1=result1[which(!duplicated(result1)),]
  cytotoxicity_fea=result1[1:20,1]
  
  TNK_cd8t <- AddModuleScore(
    object = TNK_cd8t,
    features = list(cytotoxicity_fea),
    assay	="RNA",
    ctrl = 5,
    name = 'Exhuastion_HAVCR2'
  )
  
  score$cancer_con=score$source
  score$cancer_con[which(score$cancer_con %in% c("OLK","OLP","OSF"))]="Precancerous"
  my_comparisons <- list( c("Normal", "OLK"), c("OSF",    "Normal"),c("Tumor","Normal"), c("OSF","OLK"),c("Tumor","OLK"),c("OSF","Tumor") )
  
  if(length(output_dir) == 0){
    output_dir=getwd()}
  
  score$cancer_con=score$source
  p=ggplot(score, aes(x=cancer_con, y=Cytotoxic_GZMK1,fill=cancer_con)) + 
    geom_boxplot()+
    theme_classic()+ scale_fill_brewer(palette="Oranges") + theme_classic()+ theme(legend.position="none")+theme(axis.text.x = element_text(angle = 45))+  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+ 
    stat_compare_means(label.y = 4.5) 
  ggsave("Exhausted_HAVCR2_Integrate_top20_source_corr.pdf",p,width = 3,height = 3.5)
  
  p=ggplot(score, aes(x=cancer_con, y=Exhuastion_HAVCR21,fill=cancer_con)) + 
    geom_boxplot()+
    theme_classic()+ scale_fill_brewer(palette="Oranges") + theme_classic()+ theme(legend.position="none")+theme(axis.text.x = element_text(angle = 45))+  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+ 
    stat_compare_means(label.y = 4.5) 
  ggsave("Cytotoxic_GZMK_Integrate_top20_source_corr.pdf",p,width = 3,height = 3.5)
}
