#!/bin/bash
#SBATCH --job-name=ascat
#SBATCH --partition=64c512g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH -o wgs_ascat.%j.out

path=$(pwd)
sample=$1
ROI=$2
gd=$3

cd ./${ROI}

cp /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/codes/ascat.R ./ascat/

eval "$(conda shell.bash hook)"
conda activate wgs

##### ASCAT #####
Rscript ./ascat/ascat.R $sample $ROI $gd

conda deactivate
