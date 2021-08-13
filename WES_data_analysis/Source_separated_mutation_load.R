#' Calculate the number of gene mutation in each sample and using dotplot to visualize the source seperated mutation load
#' @param Mut_ann_file Path of the ANNOVAR annotation output file
#' @param Fig_wid Width of pdf file
#' @param Fig_Hei Height of pdf file
#' @param output_dir the path of pdf file save
#' The output is a PDF file
#' @author Ya Han

library(readr)
library(ggsci)
library(ggplot2)

Source_separated_mutation_load<-function(Mut_ann_file,Fig_wid=3,Fig_Hei=6,output_dir=NULL){
  
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
  
  mutation_load=t(result)
  mutation_load=as.data.frame(mutation_load)
  colnames(mutation_load)[1:2]=c("Sample","Mutation_load")
  mutation_load$Mutation_load=as.numeric(as.vector(unlist(mutation_load$Mutation_load)))
  mutation_load$Patient=unlist(lapply(strsplit(mutation_load[,1],"_"),function(x) x[[1]]))
  mutation_load$Source=unlist(lapply(strsplit(mutation_load[,1],"_"),function(x) x[[2]]))
  
  mutation_load$Source_can="Precancer"
  mutation_load$Source_can[which(mutation_load$Source == "Tumor")]="Tumor"
  mutation_load$Source_can[which(mutation_load$Source == "Normal")]="Normal"
  
  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }
  
  p<-ggplot(mutation_load, aes(x=Source_can, y=Mutation_load, fill=Source_can)) +
    geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.8)+ theme_bw() +scale_fill_brewer(palette="Dark2")+scale_colour_brewer(palette="Dark2")+ stat_summary(fun.data=data_summary, color="blue")
  if(length(output_dir) == 0){
    output_dir=getwd()
  }
  setwd(output_dir)
  ggsave("Patient_Source_mut_gene_load.pdf",p,width = Fig_wid,height = Fig_Hei)  

}