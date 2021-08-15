#Do functional enrichment analysis and then visualization the result by heatmap
#' @param UpDEG_list     The upregulated gene list.
#' @param gmt_path       The path of gmt file.
#' @param output_dir     The path of plot file store.


library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)
library(org.Hs.eg.db)

Enrichment_heatmp<-function(UpDEG_list,gmt_path,output_dir=NULL){
  
  geneset_name=c("h.all.v7.0.entrez.gmt","c1.all.v7.0.entrez.gmt","c2.all.v7.0.entrez.gmt","c3.all.v7.0.entrez.gmt","c4.all.v7.0.entrez.gmt","c5.all.v7.0.entrez.gmt","c6.all.v7.0.entrez.gmt","c7.all.v7.0.entrez.gmt")
  description=c("Hallmark gene sets","Positional gene sets","Curated gene sets","Motif gene sets","Computational gene sets","GO gene sets","Oncogenic signatures","Immunologic signatures")
  genesets=matrix(c(geneset_name,description),ncol=2,byrow=F)
  #DEG_gene=seuratObj.markers
  DEG_gene=UpDEG_list
  
  allcluster_annootation<<-NULL
  cluster_batch_enrich<-function(temp_clu){
    gene_set=DEG_gene$gene[which(DEG_gene$cluster == temp_clu)]
    gene_set <- bitr(gene_set, fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
    
    all_enrichment_result<<-NULL
    enrichment_analysis<-function(x){
      current_gene_set=read.gmt(paste0(gmt_path,genesets[x,1]))
      enrichment_result=enricher(gene_set$ENTREZID, TERM2GENE=current_gene_set)
      if ( dim(enrichment_result)[1] >0 ) {
        enrichment_result=cbind(as.data.frame(enrichment_result)[,c(1,3,5,7)],collection=genesets[x,2])
        all_enrichment_result<<-rbind(all_enrichment_result,enrichment_result)}
    }
    result=apply(matrix(c(1,3,6)),1,enrichment_analysis)
    allcluster_annootation<<-rbind(allcluster_annootation,cbind(all_enrichment_result,temp_clu))
  }
  annotation_result=apply(matrix(as.vector(unlist(DEG_gene$cluster))),1,cluster_batch_enrich)
  
  
  
  if(length(output_dir) ==0){output_dir=getwd()}
  setwd(output_dir)
  
  allcluster_annootation=as.data.frame(allcluster_annootation)
  
  function_term=unique(allcluster_annootation$ID[which(allcluster_annootation$collection == "Hallmark gene sets")])
  temp_ano=allcluster_annootation[which(allcluster_annootation[,1] %in% function_term),]
  colnames(temp_ano)[6]="Cluster" 
  temp_ano=as.data.frame(temp_ano)
  
  heatmap_annotation<<-matrix(rep(0,length(unique(function_term))*length(as.vector(unlist(DEG_gene$cluster)))),ncol=length(as.vector(unlist(DEG_gene$cluster))))
  
  
  Generate_heatmap_data<-function(temp_clu){
    heatmap_annotation[match(temp_ano$ID[which(temp_ano$Cluster == temp_clu)],function_term),match(temp_clu,as.vector(unlist(DEG_gene$cluster)))]<<-as.numeric(as.vector(unlist(temp_ano$qvalue[which(temp_ano$Cluster == temp_clu)])))
  }
  result=apply(matrix(as.vector(unlist(DEG_gene$cluster))),1,Generate_heatmap_data)
  
  heatmap_annotation=-log10(heatmap_annotation)
  heatmap_annotation[which(is.infinite(heatmap_annotation))]=0
  colnames(heatmap_annotation)=as.vector(unlist(DEG_gene$cluster))
  rownames(heatmap_annotation)=function_term
  rownames(heatmap_annotation)=tolower(rownames(heatmap_annotation))
  rownames(heatmap_annotation)=gsub("_"," ",rownames(heatmap_annotation))
  heatmap_annotation[which(heatmap_annotation > 10)]=5
  heatmap_annotation[which(heatmap_annotation > 5)]=5
  
  pdf("upDEG_enrichment_hallmark.pdf",height = 15,width = 15)
  a=pheatmap(heatmap_annotation, cluster_rows = F,cluster_cols = F,color=colorRampPalette(c("white","#86659A"))(10),cellwidth=40,cellheight =20)
  dev.off()
}

