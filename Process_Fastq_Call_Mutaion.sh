# @Author: Ya Han
# @Dat: 2020-08-11
# @Last Modified by: Ya
# @Last Modified time: 2020-09_07 15:07

# Note:
# The analysis starts from a fastq file, and then perform quality control for fastq file, mapping, quality control for bam file, mutation calling and mutation filtering.


#acivate conda environment
source activate /fs/home/hanya/miniconda3/envs/MAESTRO


# run the following command for detailed usage of tools
fastqc --help
bwa --help
samtools --help
gatk --help


# the following commands can be used for format conversion
fastqc                           #FastQC process to see the quanlity
bwa mem                          #Mapping.
samtools                         #Convert sam file to bam file.
gatk AddOrReplaceReadGroups      #Add appropriately label read group (@RG) fields, coordinate sort and index a BAM file.
gatk FixMateInformation          #Verify mate-pair information between mates and fix if needed.
gatk BaseRecalibrator            #Builds a model of covariation based on the input data and a set of known variants, producing a recalibration file.
gatk ApplyBQSR                   #Adjusts the base quality scores in the data based on the model.
gatk Mutect2                     #Call somatic short variants.
gatk FilterMutectCalls           #Remove false positives variants.


# Note:
# If there're multiple samples (fastq file) in your analysis, usually we need process sample seperately. 

# Example
# FastQC v0.11  https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc
# Run FastQC process to see the quanlity
fastqc -o /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/ -t 4 /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF.fq



# BWA Version: 0.7.12-r1039 https://github.com/lh3/bwa
# Reference /Reference_data/GRCh38.d1.vd1.fa.tar.gz Download from https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
# Run mapping
bwa mem -t 10 /mnt/Storage/home/hanya/Reference_data/Reference_Samtools/GRCh38.d1.vd1.fa /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF.fq.gz > /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF.sam



#Samtools Version: samtools-1.5  https://sourceforge.net/projects/samtools/files/samtools/1.10/
#Samtools used to convert sam file to bam file. The index file using the GRCh38.d1.vd1.fa.tar.gz Download from https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files
samtools import /mnt/Storage/home/hanya/Reference_data/Reference_Samtools/GRCh38.d1.vd1.fa.fai  /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF.sam /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF.bam



#Generate a output report includes mapping statistics such as Mapped reads.
#For the mapping result unique mapped reads is not suit for pair end sequence to select the high quality mapping reads.
#considering the insert and deletion and gap, the mapping uniqueness is not used to select
#Then each align tool give a score used to access the mapping quality.
#MapQ = -10 log10(P), cutoff 15 means 0.05% occurs non-unique mapped, cutoff 0.01%, But many people choose 30 as the cutoff.  
samtools view -bhS -q 30 /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF.bam > /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_MapQ30.bam
	
#Mapping quality processing will remove some unaligned read, and this will make the error in GATK "Read is marked as paired, but its pair was not found".
#Before GATK call varitation should use gatk ValidateSamFile to check the format of bam.
#After Mapping quality control will appear ERROR:MATE_NOT_FOUND. this could use gatk FixMateInformation command to repair.
#gatk FixMateInformation： Verify mate-pair information between mates and fix if needed. This tool ensures that all mate-pair information is in sync between each read and its mate pair





#The input bam file of GATK must combined with a .bai(after samtools index, The samtool sort is a necessay index of samtool index)
#Before run gatk any other tool the ValidateSamFile should be run to see the sam file format.
#The error MISSING_READ_GROUP meaning the GATK input bam file need Read group or @RG/SM
#MarkDuplicates need more than 12 hours for one sample, and without option for threads selection.
#AddOrReplaceReadGroups to appropriately label read group (@RG) fields, coordinate sort and index a BAM file.
/mnt/Storage/home/hanya/gatk-4.1.8.1/gatk AddOrReplaceReadGroups -I /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_MapQ30.bam -O  /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_MapQ30_RG.bam -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM P01_OSF --SORT_ORDER coordinate --CREATE_INDEX true
/mnt/Storage/home/hanya/gatk-4.1.8.1/gatk FixMateInformation -I /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_MapQ30_RG.bam -O /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_MapQ30_RG_Fix_RP.bam



#second step is to improve the quality of gapped alignment via indel realignment.
#Some aligners (such as Novoalign) and variant callers (such as GATK HaplotypeCaller Mutect2) involve indel alignment improvement. After indel realignment, BQSR (BaseRecalibrator from the GATK suite) is recommended to improve the accuracy of base quality scores prior to variant calling.
#Note that indel realignment is no longer necessary for variant discovery if you plan to use a variant caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2. However it is still required when using legacy callers such as UnifiedGenotyper or the original MuTect.
#BQSR: it is a data pre-processing step that detects systematic errors made by the sequencing machine when it estimates the accuracy of each base call.
/mnt/Storage/home/hanya/gatk-4.1.8.1/gatk BaseRecalibrator -I /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_MapQ30_RG_Fix_RP.bam -R /mnt/Storage/home/hanya/Reference_data/Reference_Samtools/GRCh38.d1.vd1.fa --known-sites /mnt/Storage/home/hanya/Reference_data/GATK_ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --known-sites /mnt/Storage/home/hanya/Reference_data/GATK_ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz -O /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_MapQ30_.table
gatk ApplyBQSR -R mnt/Storage/home/hanya/Reference_data/Reference_Samtools/GRCh38.d1.vd1.fa -I /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_MapQ30_RG_Fix_RP.bam --bqsr-recal-file /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_MapQ30_.table -O /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_MapQ30_RG_Fix_RP_BQSR.bam





#gatk-4.1.8.1   https://github.com/broadinstitute/gatk/releases/download/4.1.8.1/gatk-4.1.8.1.zip
#GATK reference using same with BWA ans Samtools and need to create a .dict file
#gatk CreateSequenceDictionary -R /mnt/Storage/home/hanya/Reference_data/Reference_Samtools/GRCh38.d1.vd1.fa -O /mnt/Storage/home/hanya/Reference_data/Reference_Samtools/GRCh38.d1.vd1.dict
#need some human data to filter SNP and find really SNV. https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=
#GATK call SNV with two step: Mutect2 was used to Call somatic short variants ;FilterMutectCalls to remove false positives.
#FilterMutectCalls now requires a reference, which should be the same fasta file input to Mutect2.

/mnt/Storage/home/hanya/gatk-4.1.8.1/gatk Mutect2 -R /mnt/Storage/home/hanya/Reference_data/Reference_Samtools/GRCh38.d1.vd1.fa -I {non_normalfile} -I {normal_file} -tumor /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_MapQ30_RG_Fix_RP.bam -normal /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_Blood_MapQ30_RG_Fix_RP.bam --germline-resource /mnt/Storage/home/hanya/Reference_data/GATK_ref/somatic-hg38_af-only-gnomad.hg38.vcf.gz  --panel-of-normals /mnt/Storage/home/hanya/Reference_data/GATK_ref/somatic-hg38_1000g_pon.hg38.vcf.gz -O /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_vcf.gz
/mnt/Storage/home/hanya/gatk-4.1.8.1/gatk FilterMutectCalls -V /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_vcf.gz -R /mnt/Storage/home/hanya/Reference_data/Reference_Samtools/GRCh38.d1.vd1.fa -O /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/P01_OSF_filtered.vcf.gz


