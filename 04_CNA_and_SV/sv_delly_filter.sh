#!/bin/bash
#SBATCH --job-name=delly
#SBATCH --partition=64c512g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --output=%j_delly.out

eval "$(conda shell.bash hook)"
conda activate bcftools

if [ ! -d "./delly/" ];then
    mkdir ./delly
else
    echo "delly present"
fi

sample=P10

sample_array=(S3_Acinar1 S3_CA1 S3_Duct1 S3_LG1 S3_HG3 S3_HG2 S3_HG4 S3_HG1 S3_LG2 S3_LG3 S43_CA1 S43_Duct1 S43_LG1 S43_HG3 S43_HG1 S43_LG2 S43_HG2 S43_LG3 S43_LG4 S43_LG5)

for ROI in "${sample_array[@]}"; do

	echo -e "${sample}_${ROI}\ttumor\n${sample}_NOR_ref\tcontrol" > "./delly/${ROI}_info.tsv"

	singularity exec --bind /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/${sample}/:/data --bind /dssg/home/acct-guowei_med/anzhijie/pdac3d/references:/references /dssg/home/acct-guowei_med/anzhijie/pdac3d/singularity/delly-v2.5.1.sif delly filter -p -f somatic -o /data/delly/${sample}_${ROI}.somatic.bcf -a 0.10 -s /data/delly/${ROI}_info.tsv /data/delly/${sample}_${ROI}.bcf


done
conda deactivate

conda activate wgs
for ROI in "${sample_array[@]}"; do
	python /dssg/home/acct-guowei_med/anzhijie/pdac3d/software/delly2bnd.py -v ./delly/${sample}_${ROI}.somatic.bcf -r /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/hg38_index/hg38.fa -o ./delly/${sample}_${ROI}.somatic.bnd.vcf
done
conda deactivate

conda activate bcftools
for ROI in "${sample_array[@]}"; do
	bcftools view -O v -o ./delly/${sample}_${ROI}.somatic.vcf ./delly/${sample}_${ROI}.somatic.bcf
	bcftools sort -O v -o ./delly/${sample}_${ROI}.somatic.bnd.sorted.vcf ./delly/${sample}_${ROI}.somatic.bnd.vcf
done

conda deactivate
