# @Author: Ya Han
# @Dat: 2021-04-28
# @Last Modified by: Ya
# @Last Modified time: 2021-04-28

# Note:
# Lisa is used to identify the drive TF.



import os
os.chdir('/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/Without_P12OLP/CCA_result/Nfeatures5000_includingInteract_CCA/Seperate_cell_lineage/B_Plasma/Gene_regulatory')
for file in os.listdir('/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/Without_P12OLP/CCA_result/Nfeatures5000_includingInteract_CCA/Seperate_cell_lineage/B_Plasma/Gene_regulatory'):
	if(file.find("DEGs") > -1):
		print(file)
		cmd1='lisa oneshot hg38 {input_file} -b 400 -c 5 --seed=2556 --save_metadata> {output_file}'.format(input_file=file,output_file=file.replace(".txt","_lisa.tsv"))
		print(cmd1)
		os.system(cmd1)
