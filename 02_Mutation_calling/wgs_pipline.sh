#!/bin/bash
#SBATCH --job-name=wgs_tumor
#SBATCH --partition=64c512g
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH -o job_pipline.%j.out 

path=$(pwd)
sample=$1
ROI=$2

cd ./${ROI}

if [ ! -d "./qualimap/" ];then
    mkdir ./qualimap
else
    echo "qualimap present"
fi
if [ ! -d "./mutect2/" ];then
    mkdir ./mutect2
else
    echo "mutect2 present"
fi
if [ ! -d "./strelka/" ];then
    mkdir ./strelka
else
    echo "strelka present"
fi
if [ ! -d "./varscan2/" ];then
    mkdir ./varscan2
else
    echo "varscan2 present"
fi

if [ ! -d "./ascat/" ];then
    mkdir ./ascat
    cp /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/codes/ascat.R ./ascat/
else
    echo "ascat present"
    cp /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/codes/ascat.R ./ascat/
fi


eval "$(conda shell.bash hook)"
conda activate wgs

if [[ -f "${sample}_${ROI}_1_fastp_tiled.fq.gz" ]]; then
    fastq_file1="${sample}_${ROI}_1_fastp_tiled.fq.gz"
    fastq_file2="${sample}_${ROI}_2_fastp_tiled.fq.gz"
else
    fastq_file1="${sample}_${ROI}_1_fastp.fq.gz"
    fastq_file2="${sample}_${ROI}_2_fastp.fq.gz"
fi


# run bwa
bwa mem -t 24 -R '@RG\tID:'$sample'_'$ROI'_'$lane'\tPL:illumina\tSM:'$sample'_'$ROI'' /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/hg38_index/hg38.fa ${fastq_file1} ${fastq_file2}  | samtools view -bS -q 20  | \
 samtools sort -m 20G -T ${sample}_${ROI}.tmp -o ${sample}_${ROI}.sort.bam
# run sambamba
sambamba markdup -r -t 24 -p --overflow-list-size=8000000 --sort-buffer-size=10240 --tmpdir=./tmp ${sample}_${ROI}.sort.bam ${sample}_${ROI}.rmdup.bam
rm ${sample}_${ROI}.sort.bam

##### BQSR #####
gatk BaseRecalibrator \
   -I ${sample}_${ROI}.rmdup.bam \
   -R /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/hg38_index/hg38.fa \
   --known-sites /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/BQSR_known_vcf/Homo_sapiens_assembly38.dbsnp138.vcf \
   --known-sites /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/BQSR_known_vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O ./recal_data.table
gatk ApplyBQSR \
   -R /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/hg38_index/hg38.fa \
   -I ${sample}_${ROI}.rmdup.bam \
   --bqsr-recal-file ./recal_data.table \
   -O ${sample}_${ROI}.rmdup.bqsr.bam

rm ${sample}_${ROI}.rmdup.bam
# samtools index
samtools index ${sample}_${ROI}.rmdup.bqsr.bam


##### Mutect2 #####
gatk --java-options "-Xmx30g" Mutect2 -R /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/hg38_index/hg38.fa -I ${sample}_${ROI}.rmdup.bqsr.bam -tumor ${sample}_${ROI} -I ../NOR_ref/${sample}_NOR_ref.rmdup.bqsr.bam -normal ${sample}_NOR_ref \
 --panel-of-normals /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/1000g_pon.hg38.vcf.gz -O ./mutect2/${sample}_${ROI}_mutect2_raw.vcf.gz --native-pair-hmm-threads 24 --f1r2-tar-gz ./mutect2/f1r2.tar.gz --callable-depth 8
# LearnReadOrientationModel
gatk LearnReadOrientationModel -I ./mutect2/f1r2.tar.gz -O ./mutect2/read-orientation-model.tar.gz
# GetPileupSummaries
gatk GetPileupSummaries \
 -I ${sample}_${ROI}.rmdup.bqsr.bam \
 -V /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/gatk_GetPileupSummaries/small_exac_common_3.hg38.vcf.gz \
 -L /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/gatk_GetPileupSummaries/small_exac_common_3.hg38.vcf.gz \
 -O ./mutect2/pileups.table
# Estimate contamination
gatk CalculateContamination -I ./mutect2/pileups.table -matched ../NOR_ref/nor_pileups.table  -tumor-segmentation ./mutect2/segments.table -O ./mutect2/contamination.table
# FilterMutectCalls
gatk --java-options "-Xmx8g" FilterMutectCalls -R /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/hg38_index/hg38.fa -V ./mutect2/${sample}_${ROI}_mutect2_raw.vcf.gz -O ./mutect2/${sample}_${ROI}_mutect2.vcf.gz \
 -XL /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/DAC_Exclusion_List_Regions/ENCFF356LFX.bed --tumor-segmentation ./mutect2/segments.table --contamination-table ./mutect2/contamination.table --ob-priors ./mutect2/read-orientation-model.tar.gz

##### Strelka2 #####
/dssg/home/acct-guowei_med/anzhijie/pdac3d/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam ../NOR_ref/${sample}_NOR_ref.rmdup.bqsr.bam --tumorBam ${sample}_${ROI}.rmdup.bqsr.bam \
 --referenceFasta /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/hg38_index/hg38.fa --runDir ./strelka
./strelka/runWorkflow.py -m local -j 24


##### Varscan2 #####
samtools mpileup -q 1 -f /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/hg38_index/hg38.fa ${sample}_${ROI}.rmdup.bqsr.bam > ./varscan2/${sample}_${ROI}.mpileup

# tumor normal comparison
java -jar /dssg/home/acct-guowei_med/anzhijie/pdac3d/software/VarScan2/VarScan.v2.4.6.jar somatic ../NOR_ref/${sample}_NOR_ref.mpileup ./varscan2/${sample}_${ROI}.mpileup  ./varscan2/varscan2 \
 --min-coverage 5 --min-var-freq 0.1 --output-vcf 1

# processSomatic
java -jar /dssg/home/acct-guowei_med/anzhijie/pdac3d/software/VarScan2/VarScan.v2.4.6.jar processSomatic ./varscan2/varscan2.snp.vcf
java -jar /dssg/home/acct-guowei_med/anzhijie/pdac3d/software/VarScan2/VarScan.v2.4.6.jar processSomatic ./varscan2/varscan2.indel.vcf

# false positive filter
egrep -hv "^#" ./varscan2/varscan2.snp.Somatic.hc.vcf | awk 'OFS="\t" {print $1, $2-1, $2+1}' | sort -k1,1 -k2,2n | /dssg/home/acct-guowei_med/anzhijie/pdac3d/software/bedtools merge -i - > ./varscan2/somatic_hc.bed
egrep -hv "^#" ./varscan2/varscan2.indel.Somatic.hc.vcf | awk 'OFS="\t" {print $1, $2 - ((length($4) > length($5)) ? length($4) : length($5)) , $2 + ((length($4) > length($5)) ? length($4) : length($5)) }' | sort -k1,1 -k2,2n | \
 /dssg/home/acct-guowei_med/anzhijie/pdac3d/software/bedtools merge -i -  >  ./varscan2/indel.bed

rm ./varscan2/${sample}_${ROI}.mpileup
