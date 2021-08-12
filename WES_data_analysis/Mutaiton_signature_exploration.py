# @Author: Ya Han
# @Dat: 2020-01-28
# @Last Modified by: Ya
# @Last Modified time: 2020-02_28

# Note:
# The analysis starts from a vcf file , and then exploring patterns of small mutational events 


#conda activate py36
#from SigProfilerMatrixGenerator import install as genInstall
#genInstall.install('GRCh38')

import os
import pandas as pd
import numpy as np
from SigProfilerExtractor import sigpro as sig

os.chdir("/mnt/Storage2/home/hanya/project/Carcinogenesis/WES/7filtered_vcf/Mutect_output/VCF_FILE/input")

#This step is find the signature
data = "/mnt/Storage2/home/hanya/project/Carcinogenesis/WES/7filtered_vcf/Mutect_output/VCF_FILE/input"
sig.sigProfilerExtractor("vcf", "HNSC_WES_MutationSig", data, minimum_signatures=1, maximum_signatures=25,reference_genome="GRCh38",opportunity_genome="GRCh38",cpu=20)
