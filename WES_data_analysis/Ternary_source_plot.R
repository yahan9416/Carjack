#' Visualization gene mutation in different tisssue type by Ternary 
#' @param Mut_ann_file Path of the ANNOVAR annotation output file
#' @param Gene a vector of gene will be show, and if the parameter is null, will visualization the top frequent mutated 20 genes.
#' @param output_dir the path of pdf file save
#' The output is a PDF file
#' @author Ya Han

#install.packages("ggtern") 
library("ggtern")
library("ggplot2")
library(RColorBrewer)
library(readr)

Ternary_source_plot<-function(Mut_ann_file,Gene=NULL,output_dir=NULL){
  
  file=list.files(Mut_ann_file)
  setwd(Mut_ann_file)
  if(length(Gene) == 0){
    sample_mutation<<-NULL
    batch_process_muitation<-function(x){
      temp=read.csv(x,header = T)
      index=which(temp$'ExonicFunc.refGene' %in% c(".","synonymous SNV","unknown"))
      mut_gene=unique(as.vector(unlist(unique(temp$'Gene.refGene'[-1*index]))))
      sample=strsplit(x,"_filte")[[1]][1]
      sample_mutation<<-rbind(sample_mutation,cbind(sample,mut_gene))
    }
    result=apply(as.matrix(file),1,batch_process_muitation)
    Gene=names(sort(table(sample_mutation[,2]),decreasing = T))[1:20]
    gene_list_length=length(Gene)
    
  }else{
    gene_list_length=length(Gene)
  }
  
  Alteration_mat<<-matrix(rep("  ",length(file)*gene_list_length),nrow=gene_list_length)
  batch_process_muitation<-function(x){
    temp=read.csv(file[x],header = T)
    index=which(temp$'ExonicFunc.refGene' %in% c(".","synonymous SNV","unknown"))
    temp=temp[-1*index,]
    intersect_gene=intersect(Gene,temp$'Gene.refGene')
    
    temp$'ExonicFunc.refGene'[which(temp$'ExonicFunc.refGene' == "frameshift deletion")]="Frame_Shift_Del"
    temp$'ExonicFunc.refGene'[which(temp$'ExonicFunc.refGene' == "nonsynonymous SNV")]="Missense_Mutation"
    temp$'ExonicFunc.refGene'[which(temp$'ExonicFunc.refGene' == "frameshift insertion")]="Frame_Shift_Ins"
    temp$'ExonicFunc.refGene'[which(temp$'ExonicFunc.refGene' == "stopgain")]="Nonsense_Mutation"
    temp$'ExonicFunc.refGene'[which(temp$'ExonicFunc.refGene' == "nonframeshift deletion")]="In_Frame_Del"
    temp$'ExonicFunc.refGene'[which(temp$'ExonicFunc.refGene' == "nonframeshift substitution")]="Inframe_INDEL"
    temp$'ExonicFunc.refGene'[which(temp$'ExonicFunc.refGene' == "nonframeshift insertion")]="In_Frame_Ins"
    
    generate_alter<-function(y){
      index=which(temp$'Gene.refGene' == intersect_gene[y])
      value=paste0(unique(temp$'ExonicFunc.refGene'[index]),collapse = "",sep=";")
      gene_index=match(intersect_gene[y],Gene)
      Alteration_mat[gene_index,x]<<-value
    }
    result=apply(matrix(1:length(intersect_gene)),1,generate_alter)
    
  }
  result=apply(as.matrix(1:length(file)),1,batch_process_muitation)
  
  
  #最后我们要展示的gene有
  rownames(Alteration_mat)=Gene
  colnames(Alteration_mat)=files
  
  
  Normal_index=grep("Normal",colnames(Alteration_mat))
  Tumor_index=grep("Tumor",colnames(Alteration_mat))
  Precancer_index=setdiff(1:26,c(Normal_index,Tumor_index))
  
  data_plot<<-NULL
  generat_ggtern_data<-function(x){
    index=match(x,rownames(Alteration_mat))
    ave_fre=length(which(Alteration_mat[index,] != "  "))/26
    norm_fre=length(which(Alteration_mat[index,Normal_index] != "  "))/3
    precan_fre=length(which(Alteration_mat[index,Precancer_index] != "  "))/length(Precancer_index)
    tumor_fre=length(which(Alteration_mat[index,Tumor_index] != "  "))/length(Tumor_index)
    enrichment=c("Normal","Precancer","Tumor")[which.max(c(norm_fre,precan_fre,tumor_fre))]
    data_plot<<-rbind(data_plot,c(norm_fre,precan_fre,tumor_fre,ave_fre))
    return(enrichment)
  }
  result=apply(as.matrix(Gene),1,generat_ggtern_data)
  
  data_plot<-as.data.frame(data_plot)
  data_plot<<-cbind(data_plot,result)
  colnames(data_plot)<-c("Normal","Precancer","Tumor","Average","Enrichment")
  rownames(data_plot)=Gene
  data_plot<-cbind(Gene,data_plot)
  colnames(data_plot)<-c("Gene","Normal","Precancer","Tumor","Average","Enrichment")
  nv=0.05
  pn=position_nudge_tern(y=nv,x=-nv/2,z=-nv/2)
  
  p<-ggtern(data=data_plot,aes(x=Normal,y=Precancer,z=Tumor))+     #X,Y,Z轴分别映射三分组表达值
    geom_point(aes(size=Average,color=Enrichment),alpha=0.8)+     #以点图形式呈现，大小映射的是平均值，颜色映射的是最大值对应组即enrichment）
    scale_color_brewer(palette="Dark2")+   #自定义设置颜色
    theme_rgbw(base_size = 12 )+   #设置背景样式和字体大小，背景样式可通过theme_bw查询
    labs(title = "Ternary plot")+   #设置标题
    geom_text(aes(label=Gene),size=3)+
    geom_mask()+
    theme(plot.title = element_text(size=15,hjust = 0.5)) ##标题大小和位置
  # + theme_legend_position(x="topright")   #更改图注位置,topleft, middleleft, bottomleft, topright, middleright and bottomright
  if(length(output_dir) == 0){
    output_dir=getwd()
  }
  ggsave("Ternary_plot_mutation.pdf",p)
}
