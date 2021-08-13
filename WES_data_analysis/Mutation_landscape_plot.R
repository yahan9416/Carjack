#' Visualization gene mutation by Complexheatmap 
#' @param Mut_ann_file Path of the ANNOVAR annotation output file
#' @param Gene a vector of gene will be show, and if the parameter is null, will visualization the top frequent mutated 20 genes.
#' @param Clinic_file a String of clinic information file and each columns represent a feature need to be show in heatmap
#' @param output_dir the path of pdf file save
#' The output is a complexheatmap
#' @author Ya Han

library(ComplexHeatmap)
library(ggplot2)
library(readr)

Mutation_landscape_plot<-function(Mut_ann_filev,Gene=NULL,Clinic_file=NULL,output_dir=NULL){
  
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
  #Sample=unlist(lapply(strsplit(file,"_"),function(x){ paste(x[1],x[2],sep  ="_")}))
  #index=match(sort(Sample),Sample)
  #Alteration_mat=Alteration_mat[,index]
  colnames(Alteration_mat)=file
  #Patient=unlist(lapply(strsplit(Sample[index],"_"),function(x){ x[1]}))
  
  col = c("Frame_Shift_Del" = "#FF7F00", "Missense_Mutation" = "#1F78B4", "Frame_Shift_Ins" = "#E31A1C","Nonsense_Mutation"="#B2DF8A","In_Frame_Del"="#A6CEE3","Inframe_INDEL"="#FDBF6F","In_Frame_Ins"="#FB9A99")
  column_title = "OncoPrint for Oral cancer gene alteration"
  heatmap_legend_param = list(title = "Alternations", at = c("Frame_Shift_Del", "Missense_Mutation", "Frame_Shift_Ins","Nonsense_Mutation","In_Frame_Del","Inframe_INDEL","In_Frame_Ins"), 
                              labels = c("Frame_Shift_Del", "Missense_Mutation", "Frame_Shift_Ins","Nonsense_Mutation","In_Frame_Del","Inframe_INDEL","In_Frame_Ins"))
  
  alter_fun = list(
    background = alter_graphic("rect", fill = "#F4F4F2"),	
    Frame_Shift_Del = alter_graphic("rect", fill = col["Frame_Shift_Del"]),
    Missense_Mutation = alter_graphic("rect", fill = col["Missense_Mutation"]),
    Frame_Shift_Ins = alter_graphic("rect", fill = col["Frame_Shift_Ins"]),
    Nonsense_Mutation = alter_graphic("rect", fill = col["Nonsense_Mutation"]),
    In_Frame_Del = alter_graphic("rect",width = 0.9, height = 0.4, fill = col["In_Frame_Del"]),
    Inframe_INDEL = alter_graphic("rect",width = 0.9, height = 0.4, fill = col["Inframe_INDEL"]),
    In_Frame_Ins = alter_graphic("rect",width = 0.9, height = 0.4, fill = col["In_Frame_Ins"])
  )
  
  
  #Process clinic information
  if(length(Clinic_file) >0){
    Clinic_info=read.table(Clinic_file,header = T,sep="\t")
    
    Gender=Clinic_info$Gender
    Age=Clinic_info$Age
    SmokingHistory=rep("No Smoking",26)
    SmokingHistory[which(Clinic_info$SmokingHistory > 0)]="Smoking"
    AlcoholHistory=rep("No Alcohol",26)
    AlcoholHistory[which(Clinic_info$FrequencyForAlcoholUse > 0)]="Alcohol intake"
    ArecanutHistory=rep("more than 20 years",26)
    ArecanutHistory[which(Clinic_info$ArecanutHistory  < 10)]="more than 10 years"
    ArecanutHistory[which(Clinic_info$ArecanutHistory  < 5)]="less than 5 years"
    
    ha = HeatmapAnnotation(age  = anno_points(Age,gp = gpar(col = "#CD5D7D")),
                           SmokingHistory = SmokingHistory,
                           AlcoholHistory = AlcoholHistory,
                           ArecanutHistory = ArecanutHistory,
                           col = list(SmokingHistory=c("No Smoking"="#F6F6F6","Smoking"="#CD5D7D"), AlcoholHistory=c("No Alcohol"="#F6F6F6","Alcohol intake"="#F1AE89"),ArecanutHistory=c("more than 20 years"="#C3AED6","more than 10 years"="#EFBBCF","less than 5 years"="#FFEADB")),
                           annotation_legend_param = list(SmokingHistory=list(title="SmokingHistory"),AlcoholHistory=list(title="AlcoholHistory"),ArecanutHistory=list(title="ArecanutHistory"))
    )
    
    
    
    if(length(output_dir) == 0 ){
      output_dir=getwd()
    }
    setwd(output_dir)
    
    pdf("Genome_alteration_add_clin_sample_Order.pdf",width = 10,height = 7)
    oncoPrint(Alteration_mat,
              alter_fun = alter_fun,  col=col,
              top_annotation =ha,
              column_title = column_title, heatmap_legend_param = heatmap_legend_param,show_column_names=TRUE,column_order=colnames(Alteration_mat))
    dev.off()
  }else{
    
    
    if(length(output_dir) == 0 ){
      output_dir=getwd()
    }
    setwd(output_dir)
    
    pdf("Genome_alteration_add_clin_sample_Order.pdf",width = 10,height = 7)
    oncoPrint(Alteration_mat,
              alter_fun = alter_fun,  col=col,
              column_title = column_title, heatmap_legend_param = heatmap_legend_param,show_column_names=TRUE,column_order=colnames(Alteration_mat))
    dev.off()
  }
  
}