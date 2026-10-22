#!/bin/bash
#SBATCH --job-name=delly
#SBATCH --partition=64c512g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --output=%j_delly.out

eval "$(conda shell.bash hook)"
conda activate wgs

if [ ! -d "./delly/" ];then
    mkdir ./delly
else
    echo "delly present"
fi

sample=z1420838

sample_array=(S10Acinar S10CA1 S10Duct1 S10HG1 S10HG2 S10HG4 S10HG5 S10LG2 S10LG3 S10LG4 S90CA1 S90Duct1 S90HG1 S90HG2 S90HG3 S90LG1 S90LG2 S90LG3 S90LG4 S90LG5)

for ROI in "${sample_array[@]}"; do

	echo -e "P1_${ROI}\ttumor\nP1_NOR_ref\tcontrol" > "./delly/${ROI}_info.tsv"

        singularity exec --bind /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/${sample}/:/data --bind /dssg/home/acct-guowei_med/anzhijie/pdac3d/references:/references /dssg/home/acct-guowei_med/anzhijie/pdac3d/singularity/delly-v2.5.1.sif delly sr -q 25 -o /data/delly/${sample}_${ROI}.bcf -g /references/hg38_index/hg38.fa /data/${ROI}/${sample}_${ROI}.rmdup.bqsr.bam /data/NOR_ref/${sample}_NOR_ref.rmdup.bqsr.bam

	singularity exec --bind /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/${sample}/:/data --bind /dssg/home/acct-guowei_med/anzhijie/pdac3d/references:/references /dssg/home/acct-guowei_med/anzhijie/pdac3d/singularity/delly-v2.5.1.sif delly filter -f somatic -o /data/delly/${sample}_${ROI}.somatic.bcf -s /data/delly/${ROI}_info.tsv /data/delly/${sample}_${ROI}.bcf

done

conda deactivate
