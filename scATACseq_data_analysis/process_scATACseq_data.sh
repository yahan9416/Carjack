# @Author: Ya Han
# @Date:   2021-02-26
# @Last Modified by:   ya
# @Last Modified time: 2021-04-20 
# Celltype: all 
# Genome build: GRCh38
# Meta data: no
# run analysis pipeline
# Sample type: HNSCC


# Acivate conda environment
source activate /fs/home/hanya/miniconda3/envs/MAESTRO
# run the following command for detailed usage of tools.
MAESTRO scatac-init --help
macs2 callpeal --help
snakemake --help
samtools sort  --help



# Examples
MAESTRO scatac-init --platform 10x-genomics --format bam  --bam ./P11_Tumor/outs/possorted_bam.bam \
--species GRCh38 --cores 12 --directory /mnt/Storage2/home/hanya/project/Carcinogenesis/scATAC/Maestro_process_ATAC \
--outprefix Maestro --giggleannotation /mnt/Storage2/home/hanya/project/Carcinogenesis/scATAC/giggle.all/giggle.GRCh38 \
--fasta /mnt/Storage/home/hanya/Reference_data/Reference_Samtools/GRCh38.d1.vd1.fa \
--whitelist ./10X_ATAC_data/P11_Tumor/outs/raw_peak_bc_matrix/barcodes.tsv

snakemake -j 12


#Call peak
macs2 callpeak -t HW3_Q2_tumor.bed -f BED -g hs -q 0.05 -n HW3_Q2_tumor_peak



#samtools
nohup samtools sort  --threads 16 possorted_bam.bam -o possorted_bam_sort.bam &

