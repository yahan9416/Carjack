# @Author: Ya Han
# @Dat: 2021-03-20
# @Last Modified by: Ya
# @Last Modified time: 2021-03-20

# Note:
# The analysis starts from peak count matrix , and then identifying the doublets from all cells 




import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
from MAESTRO.scATAC_H5Process import *


os.chdir('/mnt/Storage2/home/hanya/project/Carcinogenesis/scATAC/Maestro_process_ATAC/Peak_count_merge_old_Meas/scrublet_Doublets')
for file in os.listdir('/mnt/Storage2/home/hanya/project/Carcinogenesis/scATAC/Maestro_process_ATAC/Peak_count_merge_old_Meas/Sample_peak_count'):
	print(file)
	scrna_count = read_10X_h5('/mnt/Storage2/home/hanya/project/Carcinogenesis/scATAC/Maestro_process_ATAC/Peak_count_merge_old_Meas/Sample_peak_count/'+file)
	counts_matrix = scrna_count.matrix
	features = scrna_count.names.tolist()
	barcodes = scrna_count.barcodes.tolist()
	scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
	doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
	np.savetxt("dbl_" + file + ".txt", scrub.predicted_doublets_, fmt='%d')


