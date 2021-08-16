#' Correct batch effect by Conos
#' @param Seuobj_files path of seurat object
#' @param Output_dir Path of save the output file

#' @author Ya Han

library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(conos)


Run_Conos<-function(Seuobj_files,Output_dir=NULL){
  seratobj_files=list.files(Seuobj_files)
  sample_name=unlist(lapply(strsplit(seratobj_files,"_"), function(x) {paste(x[3],x[4],sep="_")} ))
  seratobj_files=paste0(Seuobj_files,seratobj_files)
  
  seurat_obj_list<<-list()
  batch_readin_seuratobj<-function(x){
    seuratobj=readRDS(x)
    seuratobj$RNA=RunTSNE(seuratobj$RNA)
    seurat_obj_list<<-c(seurat_obj_list,seuratobj$RNA)
    
  }
  result=apply(as.matrix(seratobj_files),1,batch_readin_seuratobj)
  names(seurat_obj_list)=sample_name
  
  con <- Conos$new(seurat_obj_list)
  if(length(Output_dir) == 0){
    Output_dir=getwd()
  }
  setwd(Output_dir)
  
  con$buildGraph()
  # find communities
  con$findCommunities()
  # plot joint graph
  con$plotGraph()
  # plot panel with joint clustering results
  retrun(cos)
}

