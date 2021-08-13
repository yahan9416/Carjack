#' Visualization of the mutation signature load at one clinic factor
#' @param Mut_Signature Path of the significant mutation signature output file example:./HNSC_WES_Mutation_Sig_New/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt
#' @param Clinic_file A string of the clinic information file.
#' @param Fig_wid Width of pdf file
#' @param Fig_Hei Height of pdf file
#' @param output_dir the path of pdf file save
#' The output is a PDF file
#' @author Ya Han

library(ggplot2)
library("ggsci")
library(ggpubr)

Clinic_Mutationsignature<-function(Mut_Signature,Clinic_file,output_dir=NULL,Fig_wid=3,Fig_Hei=5){
  
  mutation_sign=read.table(Mut_Signature,header = T,sep="\t")
  Clinic_info=read.table(Clinic_file,header = T,sep="\t")
  
  Patient=Clinic_info$SampleID[which(Clinic_info$ArecanutHistory  == 0)]
  mutation_sign$Patien=unlist(lapply(strsplit(mutation_sign[,1],"_"),funcrion(x) x[1]))
  mutation_sign$Arecanut="Arecanut"
  mutation_sign$Arecanut[which(mutation_sign$Patien %in% Patient)]="no Arecanut"
  
  batch_separately_draw_plot<-function(x){
    p<-ggplot(mutation_sign, aes(x=Arecanut, y=mutation_sign[,x], fill=Arecanut,colors=Arecanut)) +
      geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.8)+ theme_bw() + scale_fill_npg()+ scale_color_npg()+ stat_summary(fun.data=data_summary, color="blue")
    if(length(output_dir) == 0){
      output_dir=getwd()
    }
    setwd(output_dir)
    ggsave(paste0("Patient_",x,"_load1.pdf"),p,width = Fig_wid,height = Fig_Hei)
  }
  result=apply(matrix(setdiff(colnames(mutation_sign),c("Patien","Samples","Arecanut"))),1,batch_separately_draw_plot)
  
}


