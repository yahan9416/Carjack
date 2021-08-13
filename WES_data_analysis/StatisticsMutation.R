#' Generate the matix of which gene mutation at specific sample
#' @param filepath Path of the ANNOVAR annotation output file
#' The output is a matix, row is gene and columns is sample
#' @author Ya Han

library(readr)

StatisticsMutation<-function(file.path){
  file=list.files(fold)
  #read in annovar annotation result
  sample_mutation<<-NULL
  batch_process_muitation<-function(x){
    setwd(fold)
    temp=read.csv(x,header = T)
    index=which(temp$'ExonicFunc.refGene' %in% c(".","synonymous SNV","unknown"))
    mut_gene=as.vector(unlist(unique(temp$'Gene.refGene'[-1*index])))
    sample=strsplit(x,"_filte")[[1]][1]
    print(length(mut_gene))
    #return(cbind(sample,mut_gene))
    sample_mutation<<-rbind(sample_mutation,cbind(sample,mut_gene))
    return(c(sample,length(mut_gene)))
  }
  result=apply(as.matrix(file),1,batch_process_muitation)
  
  #Convert the gene mutation into a table, 行是sample，列是基因
  unique_gene=unique(unlist(sample_mutation[,2]))
  mut_sample=unique(unlist(sample_mutation[,1]))
  Mut_matrix<<-matrix(rep("Wild_type",length(unique_gene)*length(mut_sample)),ncol=length(mut_sample))
  batch_sample_mutation_gene<-function(temp_samp_index){
    temp_gene=sample_mutation[which(sample_mutation[,1] == mut_sample[temp_samp_index]),2]
    Mut_matrix[match(temp_gene,unique_gene),temp_samp_index]]<<-"Mutation"
  }
  result=apply(matrix(1:length(mut_sample)),1,batch_sample_mutation_gene)
  
  return(Mut_matrix)
} 
