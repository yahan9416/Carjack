#' Calculate the influence of gene mutation on cell type proportion and visualization by heatmap.
#' @param Mut_ann_file Path of the ANNOVAR annotation output file
#' @param Gene A vector of gene will be show,if the number of gene mutated less than 3 sample, will be remove. If the parameter is null, will visualization the top frequent mutated 20 genes.
#' @param Meta_data_file An rds file could be a list including different cell type lineage meta-information.
#' @param Fig_wid Width of pdf file
#' @param Fig_Hei Height of pdf file
#' @param output_dir the path of pdf file save
#' The output is a PDF file
#' @author Ya Han


library(ggpubr)
library(rstatix)
library(ggplot2)
library(rstatix)
library(readr)

Mutation_Celltype_Proportion<-function(Mut_ann_file,Gene=NULL,Meta_data_file,output_dir=NULL,Fig_wid=10,Fig_Hei=5){
  
  file=list.files(Mut_ann_file)
  setwd(Mut_ann_file)
  sample_mutation<<-NULL
  batch_process_muitation<-function(x){
    temp=readr::read_csv(x)
    index=which(temp$'ExonicFunc.refGene' %in% c(".","synonymous SNV","unknown"))
    mut_gene=as.vector(unlist(unique(temp$'Gene.refGene'[-1*index])))
    sample=strsplit(x,"_filte")[[1]][1]
    print(length(mut_gene))
    #return(cbind(sample,mut_gene))
    temp_sampl_load=cbind(sample,mut_gene)
    print(dim(temp_sampl_load))
    sample_mutation<<-rbind(sample_mutation,temp_sampl_load)
    return(c(sample,length(mut_gene)))
  }
  result=apply(as.matrix(file),1,batch_process_muitation)
  
  Mutation_gene=names(which(table(sample_mutation[,2]) > 3))
  if(length(Gene) == 0){
    Gene=names(sort(table(sample_mutation[,2]),decreasing = T))
    if(length(intersect(Gene,Mutation_gene)) >20){
      Gene=intersect(Gene,Mutation_gene)[1:20]
    }else{Gene=intersect(Gene,Mutation_gene)}
    
  }else{
    Mutation_gene=intersect(Gene,Mutation_gene)
  }
  
  
  
  #check the correct meta information again and the rownames are right. (2021-01-19)
  meta.data=readRDS(Meta_data_file)
  
  
  Gene_muta_cell_type<<-NULL
  batch_celltype_dis<-function(x){
    seurat_obj=meta.data[[x]]
    temp_result=data.frame(Sample=as.vector(unlist(lapply(strsplit(rownames(seurat_obj),"@"),function(x) x[1]))),Cell_type=as.vector(unlist(seurat_obj$curated_anno_gene)))
    
    cluster_sample=table(temp_result$Sample,temp_result$Cell_type)
    source_pro=t(apply(cluster_sample,1,function(x) x/sum(x)))
    sample_pro=data.frame(Cell_type=rep(colnames(source_pro),each=dim(source_pro)[1]),Sample=rep(rownames(source_pro),dim(source_pro)[2]),Proportion=as.vector(source_pro))
    
    analysis_each_gene_mutation_state<-function(gene){
      
      TP53_mutsam=sample_mutation[which(sample_mutation[,2] == gene),1]
      sample_pro$Mut_Fre_class="Wild_type"
      sample_pro$Mut_Fre_class[which(sample_pro$Sample %in%   TP53_mutsam)]="Mutation"
      
      my_comparisons <- list( c("Mutation", "Wild_type"))
      stat.test=NULL
      stat.test <- sample_pro %>%
        group_by(Cell_type) %>%
        t_test(Proportion ~ Mut_Fre_class)
      
      temp_gene_effect=cbind(gene,stat.test)
      
      Gene_muta_cell_type<<-rbind(Gene_muta_cell_type,temp_gene_effect)
    }
    temp_result<-apply(as.matrix(Mutation_gene),1,analysis_each_gene_mutation_state)
  }
  result=apply(matrix(1:length(meta.data)),1,batch_celltype_dis)
  
  #return(Gene_muta_cell_type)
  #因为很多基因都是在很多的cell type 中是显著的，所以我们需要对基因和cell type 做筛选。 
  
  Gene_muta_cell_type$Proportion_of_cell="High"
  Gene_muta_cell_type$Proportion_of_cell[which(Gene_muta_cell_type$statistic < 0)]="Low"
  Gene_muta_cell_type$Proportion_of_cell[which(Gene_muta_cell_type$p > 0.05)]="ns"
  
  Sign_Gene_muta_cell_type=Gene_muta_cell_type
  Sign_Gene_muta_cell_type$P_value="ns"
  Sign_Gene_muta_cell_type$P_value[which(Sign_Gene_muta_cell_type$p <= 0.05)]="<= 0.05"
  Sign_Gene_muta_cell_type$P_value[which(Sign_Gene_muta_cell_type$p <= 0.01)]="<= 0.01"
  
  cell_linege=unlist(lapply(strsplit(Sign_Gene_muta_cell_type$Cell_type,"_"), function(x) x[1]))
  cell_linege[which(cell_linege %in% c("DCs","Macro","Mast","Mono","Neutrophils","pDCs"))]="Myeloid"
  cell_linege[which(cell_linege %in% c("Fibro","Myocytes","MyoFibro"))]="Fibroblast"
  cell_linege[which(cell_linege %in% c("B","Plasma"))]="B_plasma"
  cell_linege[which(cell_linege %in% c("CD4","CD8","NK","NKT","T"))]="T_NK"
  cell_linege[which(cell_linege %in% c("Malignant"))]="Mali"
  cell_linege[which(cell_linege %in% c("Epithelial"))]="Epi"
  
  if(length(output_dir) == 0){
    output_dir=getwd()
  }
  setwd(output_dir)
  
  
  S1<- ggplot(selecet_celltype_effect, aes(x= Cell_type, y=gene, size=P_value, color=Proportion_of_cell, group=Cell_type)) + geom_point() + theme_classic(base_size = 15)+
    theme(axis.text.x = element_text(angle=90,vjust=0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          strip.placement = "outside",
          strip.background = element_blank()
    )+
    facet_grid(~factor(cell_linege,levels=c("T_NK","B_plasma","Myeloid","Epith","Mali","Fibroblast","Endo")) , space="free_x", scales="free_x",switch = "x")+
    theme(strip.text.x = element_text(face="bold"))+
    scale_color_manual(values=c("#8675A9", "#70AF85", "#E8EAE6"))+
    scale_size_manual(values=c(6,4,1))
  ggsave(file="Mutation_effect_TME_celltype_Proportion.pdf",S1,width = Fig_wid,height = Fig_Hei)
  
  return(Gene_muta_cell_type)
}