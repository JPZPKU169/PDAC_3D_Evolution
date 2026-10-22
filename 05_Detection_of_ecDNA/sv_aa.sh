#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH --get-user-env
#SBATCH --partition=64c512g
#SBATCH -J aa
#SBATCH -n 20
#SBATCH --ntasks-per-node=20

export AA_SRC=/dssg/home/acct-guowei_med/anzhijie/pdac3d/software/AmpliconArchitect/src
export AC_SRC=/dssg/home/acct-guowei_med/anzhijie/pdac3d/software/AmpliconClassifier-main
export AA_DATA_REPO=/dssg/home/acct-guowei_med/anzhijie/pdac3d/software/data_repo

python /dssg/home/acct-guowei_med/anzhijie/pdac3d/software/AmpliconSuite-pipeline/GroupedAnalysisAmpSuite.py -i /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/P10/aa_input.txt -o /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/P10/AA_downsample50 -t 40 --ref GRCh38 --downsample 50
