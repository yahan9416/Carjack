# @Author: Ya Han
# @Dat: 2020-09-28
# @Last Modified by: Ya


# Note:
# The SCVI is a tool used for correct the batch effect 


#conda activate py36
import sys
import scvi
import scanpy as sc

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch


hnsccobj=sc.read_10x_h5("/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/Data/Carcinogenesis_merge_QCfiltered_gene_count.h5", genome="GRCh38",gex_only=False)
patient=hnsccobj.obs.index
#pay attention to we need to coorect the batch by patient in stead of sample.
#if we using the sample as the batch, this will remove some import bology different 
batch=[]
for i in patient:
	batch.append(i.split("_")[0])


hnsccobj.obs["batch"]=batch

hnsccobj.layers["counts"] = hnsccobj.X.copy()

sc.pp.normalize_total(hnsccobj, target_sum=10e4)
sc.pp.log1p(hnsccobj)
hnsccobj.raw = hnsccobj
sc.pp.highly_variable_genes(
    hnsccobj,
    n_top_genes=3000,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
)


scvi.data.setup_anndata(hnsccobj, layer="counts", batch_key="batch")

model = scvi.model.SCVI(hnsccobj)
model.train()
model.trainer.history # check training history
latent = model.get_latent_representation() # get latent

os.chdir("/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/SCVI_remove_batch")
model.save("my_model_patient_batch/")
#model1 = scvi.model.SCVI.load(hnsccobj, "my_model_all_cells/", use_cuda=True)

hnsccobj.obsm["X_scVI"] = model.get_latent_representation()
#hnsccobj.obsm["X_normalized_scVI"] = model.get_normalized_expression(library_size=10e4)
#aa=model.get_normalized_expression(library_size=10e4)
#aa.to_csv(r'scvi_process_top3000_feature_allCells.txt', header=True, index=True, sep='\t', mode='a')





#after SCVI correct the batch 
#reload the model

model = scvi.model.SCVI.load(hnsccobj, "/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/SCVI_remove_batch/my_model_patient_batch", use_cuda=True)
#hnsccobj.obsm["X_scVI"] = model.get_latent_representation()
hnsccobj.obsm["X_normalized_scVI"] = model.get_normalized_expression(library_size=10e4)
hnsccobj.obsm["X_scVI"] = model.get_latent_representation()
sc.pp.neighbors(hnsccobj, use_rep="X_scVI",n_pcs=100)
sc.tl.umap(hnsccobj, min_dist=0.2)
os.chdir("/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/SCVI_remove_batch")
sc.pl.umap(
    hnsccobj,
    color="batch",
    frameon=False,
    save="SCVI_coorect_scanpy_batch_umap.png"
)

sc.pl.umap(hnsccobj,
		 color=['PTPRC', 'CD3D', 'KLRD1','CD79A','IGKC','S100A9','CPA3','EPCAM','MYLK','KRT14','COL1A2','MCAM'],
		 frameon=False,
		 save="SCVI_coorect_scanpy_Marker_umap.png")