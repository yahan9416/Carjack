# @Author: Ya Han
# @Dat: 2021-07-01
# @Last Modified by: Ya


# Note:
# The analysis pipeline starts from a bam file and return a loom file


# Acivate conda environment
source activate /fs/home/hanya/miniconda3/envs/MAESTRO


# run the following command for detailed usage of tools.
velocyto run --help


# Examples:
cd /fs/home/hanya/Project/Carcinogenesis/scRNAseq/Data/P03_OSF/
nohup velocyto run -b ./barcodes.tsv -o ./  ./possorted_genome_bam.bam /fs/home/hanya/Reference/Homo_sapiens.GRCh38.104.gtf >log.txt &

