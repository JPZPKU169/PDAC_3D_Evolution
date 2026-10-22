#!/bin/bash
#SBATCH --job-name=svaba
#SBATCH --partition=64c512g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --output=%j_svaba.out

eval "$(conda shell.bash hook)"
conda activate svaba

if [ ! -d "./svaba/" ];then
    mkdir ./svaba
else
    echo "文件夹svaba已存在"
fi

sample=P10

sample_array=(S3_Acinar1 S3_CA1 S3_Duct1 S3_LG1 S3_HG3 S3_HG2 S3_HG4 S3_HG1 S3_LG2 S3_LG3 S43_CA1 S43_Duct1 S43_LG1 S43_HG3 S43_HG1 S43_LG2 S43_HG2 S43_LG3 S43_LG4 S43_LG5)


for ROI in "${sample_array[@]}"; do
	cd ./${ROI}

	/dssg/home/acct-guowei_med/anzhijie/pdac3d/software/svaba/build/svaba run -t ${sample}_${ROI}.rmdup.bqsr.bam -n ../NOR_ref/${sample}_NOR_ref.rmdup.bqsr.bam -G /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/hg38_index/hg38.fa -a ${sample}_${ROI} -p 16 --blacklist /dssg/home/acct-guowei_med/anzhijie/pdac3d/software/svaba/tracks/hg38.combined_blacklist.bed

	/dssg/home/acct-guowei_med/anzhijie/pdac3d/software/svaba/build/svaba postprocess -i ${sample}_${ROI} -t 16 -m 12G
	
	/dssg/home/acct-guowei_med/anzhijie/pdac3d/software/svaba/build/svaba tovcf -i ${sample}_${ROI}.bps.sorted.dedup.txt.gz -b ${sample}_${ROI}.rmdup.bqsr.bam -a ${sample}_${ROI}

	cd ../
done
