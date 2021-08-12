# @Author: Ya Han
# @Dat: 2020-01-07
# @Last Modified by: Ya
# @Last Modified time: 2020-02_22

# Note:
# The analysis starts from a vcf file , and then convert the VCF file into a avinput format and then using annotate_variation.pl to do annotation
# The ANNOVAR don't need install
# The reference data used in ANNOVAR, could be download by the code

# annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
# annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/ 
# annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/
# The genome reference are store in /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/annovar/humandb

import os
os.chdir('/mnt/Storage2/home/hanya/project/Carcinogenesis/WES/7filtered_vcf/')
for file in os.listdir('/mnt/Storage2/home/hanya/project/Carcinogenesis/WES/7filtered_vcf/'):
	if(file.find("filtered.vcf") > -1):
		print(file)
		cmd1='/mnt/Storage2/home/hanya/project/Carcinogenesis/WES/annovar/convert2annovar.pl -format vcf4 --allsample  -filter pass {input_file} --outfile {output_file}'.format(input_file=file,output_file=file+".filter_avinput")
		print(cmd1)
		os.system(cmd1)




os.chdir('/mnt/Storage2/home/hanya/project/Carcinogenesis/WES/7filtered_vcf/')
for file in os.listdir('/mnt/Storage2/home/hanya/project/Carcinogenesis/WES/7filtered_vcf/'):
	if(file.find("avinput") > -1):
		print(file)
		cmd1='/mnt/Storage2/home/hanya/project/Carcinogenesis/WES/annovar/table_annovar.pl {input_file} /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/annovar/humandb/ -buildver hg38 -out {output_file}  -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -csvout -polish -xref /mnt/Storage2/home/hanya/project/Carcinogenesis/WES/annovar/example/gene_fullxref.txt'.format(input_file=file,output_file=file+".annovar")
		print(cmd1)
		os.system(cmd1)
