#!/bin/bash
#SBATCH --job-name=svaba
#SBATCH --partition=64c512g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --output=%j_svaba.out

eval "$(conda shell.bash hook)"
conda activate wgs

sample=P10

sample_array=(S3_Acinar1 S3_CA1 S3_Duct1 S3_LG1 S3_HG3 S3_HG2 S3_HG4 S3_HG1 S3_LG2 S3_LG3 S43_CA1 S43_Duct1 S43_LG1 S43_HG3 S43_HG1 S43_LG2 S43_HG2 S43_LG3 S43_LG4 S43_LG5)

for ROI in "${sample_array[@]}"; do
	
	gunzip -c ./${ROI}/${sample}_${ROI}.sv.vcf.gz > ./${ROI}/${sample}_${ROI}.sv.vcf
	python /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/filter_svaba_vaf.py --input /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/${sample}/${ROI}/${sample}_${ROI}.sv.vcf \
		--output /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/${sample}/${ROI}/${sample}_${ROI}.somatic.vaf0.10.ad2.sv.vcf \
		--tumor-sample t001 \
		--min-vaf 0.10 \
		--min-ad 2 \
		--force

done
