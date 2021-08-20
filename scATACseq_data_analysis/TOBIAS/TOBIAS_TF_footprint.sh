# @Author: Ya Han
# @Dat: 2021-04-14
# @Last Modified by: Ya
# @Last Modified time: 2021-07-04

# Note:
# Pseudo-bulk ATAC-seq profiles are generated by pooling together cells within each cell type.
# After Pseudo-bulk ATAC-seq, generate the footprint figure of specific TF in cell-type bulk ATACseq data.

# Acivate conda environment
source activate /fs/home/hanya/miniconda3/envs/MAESTRO


# run the following command for detailed usage of tools.
TOBIAS ATACorrect      --help
TOBIAS BINDetect       --help
TOBIAS PlotAggregate   --help

# the following commands can be used for format conversion
TOBIAS ATACorrect          # Corrects this bias to yield a corrected signal
TOBIAS BINDetect           # generate TF_all.bed for each cell type
TOBIAS PlotAggregate       # Plot the TF bindsit signal deviation

# Examples:
cd /mnt/Storage2/home/hanya/project/Carcinogenesis/scATAC/Maestro_process_ATAC/Maestro_seperate_process_sample/CD8T_precancer/
TOBIAS ATACorrect --bam ./CD8_CCR7_Precancer.bam \
                  --genome /mnt/Storage/home/hanya/Reference_data/Reference_Samtools/GRCh38.d1.vd1.fa \
                  --peaks ./CD8_CCR7_Precancer_MA/Result/Analysis/CD8_CCR7_Precancer_final_peaks.bed \
                  --outdir CD8_CCR7_Precancer --cores 12


TOBIAS BINDetect --motifs ./MA1634.1.jaspar \
                 --signals ./CD8T_precancer/CD8_HAVCR_Precancer/CD8_HAVCR_Precancer_uncorrected.bw ./CD8T_Tumor/CD8_HAVCR_Tumor/CD8_HAVCR_Tumor_uncorrected.bw \
                 --genome  /mnt/Storage/home/hanya/Reference_data/Reference_Samtools/GRCh38.d1.vd1.fa \
                 --peaks ./CD8T_cellbarcode_seperate_analysis/CD8_HAVCR2/Result/Analysis/CD8_HAVCR2_final_peaks.bed \
                 --outdir CD8_HAVCR2_BATF_BINDetect_uncoor \
                 --cond_names CD8_HAVCR2_Precancer CD8_HAVCR2_Tumor --cores 8



cd /mnt/Storage2/home/hanya/project/Carcinogenesis/scATAC/Maestro_process_ATAC/Maestro_seperate_process_sample/
TOBIAS PlotAggregate --TFBS ./CD8_GZMK_BATF_BINDetect/BATF_MA1634.1/beds/BATF_MA1634.1_all.bed  \
                     --signals ./CD8T_Tumor/CD8_CCR7_Tumor/CD8_CCR7_Tumor_corrected.bw  ./CD8T_precancer/CD8_CCR7_Precancer/CD8_CCR7_Precancer_corrected.bw ./CD8T_Tumor/CD8_GZMK_Tumor/CD8_GZMK_Tumor_corrected.bw ./CD8T_precancer/CD8_GZMK_Precancer/CD8_GZMK_Precancer_corrected.bw ./CD8T_precancer/CD8_GZMB_Precancer/CD8_GZMB_Precancer_uncorrected.bw ./CD8T_Tumor/CD8_GZMB_Tumor/CD8_GZMB_Tumor_uncorrected.bw ./CD8T_precancer/CD8_HAVCR_Precancer/CD8_HAVCR_Precancer_uncorrected.bw ./CD8T_Tumor/CD8_HAVCR_Tumor/CD8_HAVCR_Tumor_uncorrected.bw \
                     --output ./BATF_CD8T_comparison_source.pdf --share_y both --plot_boundaries --signal-on-x
