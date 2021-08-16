# @Author: Ya Han
# @Dat: 2020-08_12
# @Last Modified by: Ya
# @Last Modified time: 2020-08_25

# Note:
# The analysis starts from count matrix , and then identifying the doublets from all cells 




import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os


os.chdir('/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/scrublet_Doublets')
for file in os.listdir('/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/Data'):
	if(os.path.isfile('/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/Data/'+file) == False):
		counts_matrix = scipy.io.mmread('/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/Data/'+file+'/filtered_feature_bc_matrix/matrix.mtx').T.tocsc()
		scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
		doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
		np.savetxt("dbl_" + file + ".txt", scrub.predicted_doublets_, fmt='%d')


