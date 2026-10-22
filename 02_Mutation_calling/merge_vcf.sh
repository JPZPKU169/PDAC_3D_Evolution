#!/bin/bash
#SBATCH --job-name=merge_vcf
#SBATCH --partition=small
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output=%j_qc.out
#SBATCH --get-user-env

sample=$1

echo "CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR" > ./${sample}_mutation_all_roi.vcf
sample_array=(S3_Acinar1 S3_CA1 S3_Duct1 S3_LG1 S3_HG3 S3_HG2 S3_HG4 S3_HG1 S3_LG2 S3_LG3 S43_CA1 S43_Duct1 S43_LG1 S43_HG3 S43_HG1 S43_LG2 S43_HG2 S43_LG3 S43_LG4 S43_LG5)
for roi in "${sample_array[@]}"; do

    cat ./"$roi"/${sample}_${roi}_merged_uniq_all_filtered2.vcf | egrep -hv "^#" >> ./${sample}_mutation_all_roi.vcf

done

cat ${sample}_mutation_all_roi.vcf | awk -F'\t' '!seen[$1 FS $2 FS $4 FS $5]++' > ${sample}_mutation_all_roi_uniq.vcf
rm ${sample}_mutation_all_roi.vcf

sed '1d' ${sample}_mutation_all_roi_uniq.vcf | awk 'OFS="\t" {print $1, $2 - ((length($4) > length($5)) ? length($4) : length($5)) , $2 + ((length($4) > length($5)) ? length($4) : length($5)) }' | sort -k1,1 -k2,2n | \
 /dssg/home/acct-guowei_med/anzhijie/pdac3d/software/bedtools merge -i -  >  ./${sample}_mutation_all_uniq.bed


for roi in "${sample_array[@]}"; do

    singularity run --bind /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/${sample}/:/data --bind /dssg/home/acct-guowei_med/anzhijie/pdac3d/references:/references /dssg/home/acct-guowei_med/anzhijie/pdac3d/singularity/bam-readcount.sif \
     -f /references/hg38_index/hg38.fa -q 20 -b 20 -d 1000 -l /data/${sample}_mutation_all_uniq.bed -w 5 /data/${roi}/${sample}_${roi}.rmdup.bqsr.bam > /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/${sample}/${roi}/${sample}_${roi}_all_mutation.bamRC

done

for roi in "${sample_array[@]}"; do

    perl ./all_mutation_count.pl -var-file ./${sample}_mutation_all_roi_uniq.vcf -readcount-file ./${roi}/${sample}_${roi}_all_mutation.bamRC -output-file ./${roi}/${sample}_${roi}_all_mutation_site_readcount.txt -roi $roi -sample $sample

done

output_file=${sample}_all_mutation_site_all_roi_readcount.txt
cp ${sample}_mutation_all_roi_uniq.vcf $output_file

for roi in "${sample_array[@]}"; do
    paste $output_file ./${roi}/${sample}_${roi}_all_mutation_site_readcount.txt > temp && mv temp "$output_file"
done
