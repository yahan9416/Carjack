# @Author: Ya Han
# @Dat: 2021-04-20
# @Last Modified by: Ya
# @Last Modified time: 2021-04-20

# Note:
# Convert scATACseq data from sce to anndata format



import scanpy as sc
import numpy as np
import anndata2ri
from rpy2.robjects import r
anndata2ri.activate()
rscript = 'readRDS("/mnt/Storage2/home/hanya/project/Carcinogenesis/scATAC/Maestro_process_ATAC/Integrate_RNAseq_scATAC/Seerated_celltype/Myeloid_subcelltype_transfor/DC/DC_ATAC_sce_object.rds")'
#rscript = 'readRDS("{RDS_file_path}")'.format(RDS_file_path = input_RDS)
adata = r(rscript)
adata.var.columns = [str(i) for i in adata.var.columns]
adata.obs.columns = [str(i) for i in adata.obs.columns]
adata.write("/mnt/Storage2/home/hanya/project/Carcinogenesis/scATAC/Maestro_process_ATAC/Integrate_RNAseq_scATAC/Seerated_celltype/Myeloid_subcelltype_transfor/DC/DC_ATAC_sce_object.h5ad")


#SCRIPT enrich -e ./DC_ATAC_sce_object.h5ad -i ./DC_ATAC_count_peak_count.h5 -s hs -p Oral_DC_scATAC
