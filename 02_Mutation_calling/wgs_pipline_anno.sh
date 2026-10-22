#!/bin/bash
#SBATCH --job-name=wgs_anno
#SBATCH --partition=small
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --output=%j_anno.out
#SBATCH --get-user-env

#set -xv
sample=$1
ROI=$2

cd ./$ROI/


# false positive filter
java -jar /dssg/home/acct-guowei_med/anzhijie/pdac3d/software/VarScan2/VarScan.v2.4.6.jar fpfilter ./varscan2/varscan2.snp.Somatic.hc.vcf ./varscan2/snp_hc.bamRC --output-file ./varscan2/varscan2.snp.Somatic.pass.vcf --filtered-file ./varscan2/varscan2.snp.Somatic.vaf.filtered.vcf \
 --min-var-count 5 --min-var-freq 0.15 --min-ref-readpos 0.20 --min-var-readpos 0.15 --min-ref-dist3 0.20 --min-var-dist3 0.15 --min-strandedness 0 --min-strand-reads 5 --min-ref-basequal 15 --min-var-basequal 20  --min-ref-dist3 0.20 --min-var-dist3 0.15 --max-rl-diff 0.25 \
 --max-var-mmqs 150 --max-ref-mmqs 100 --max-mmqs-diff 150 --min-ref-mapqual 10 --min-var-mapqual 10 --max-mapqual-diff 15
java -jar /dssg/home/acct-guowei_med/anzhijie/pdac3d/software/VarScan2/VarScan.v2.4.6.jar fpfilter ./varscan2/varscan2.indel.Somatic.hc.vcf ./varscan2/indel_hc.bamRC --output-file ./varscan2/varscan2.indel.Somatic.pass.vcf --filtered-file ./varscan2/varscan2.indel.Somatic.vaf.filtered.vcf \
 --min-var-count 5 --min-var-freq 0.15 --min-ref-readpos 0.20 --min-var-readpos 0.15 --min-ref-dist3 0.20 --min-var-dist3 0.15 --min-strandedness 0 --min-strand-reads 5 --min-ref-basequal 15 --min-var-basequal 20  --min-ref-dist3 0.20 --min-var-dist3 0.15 --max-rl-diff 0.05 \
 --max-var-mmqs 100 --max-ref-mmqs 100 --max-mmqs-diff 50 --min-ref-mapqual 10 --min-var-mapqual 10 --max-mapqual-diff 10


zcat ./mutect2/${sample}_${ROI}_mutect2.vcf.gz | grep -E '#|PASS' > ./mutect2/mutect2_filtered_once.vcf

perl /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/codes/mutect2_filter_snv_v4.pl -s ./mutect2/mutect2_filtered_once.vcf -total_depth 8 -var_depth 5 -vaf 0.15 > ./mutect2/mutect2_filtered_last_snv.vcf
perl /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/codes/mutect2_filter_indel_v4.pl -i ./mutect2/mutect2_filtered_once.vcf -total_depth 8 -var_depth 5 -vaf 0.15  > ./mutect2/mutect2_filtered_last_indel.vcf

perl /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/codes/strelka2_filter_v4.pl -s ./strelka/results/variants/somatic.snvs.vcf.gz -total_depth_snv 8 -var_depth_snv 5 -vaf_snv 0.15 > ./strelka/strelka.snv.filtered.vcf
perl /dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/codes/strelka2_filter_v4.pl -i ./strelka/results/variants/somatic.indels.vcf.gz -total_depth_indel 8 -var_depth_indel 5 -vaf_indel 0.15 > ./strelka/strelka.indel.filtered.vcf



echo "##fileformat=VCFv4.2" > ${sample}_${ROI}_merged_snv.vcf
echo "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR" >> ${sample}_${ROI}_merged_snv.vcf
awk 'NR==FNR {A[$1"\t"$2"\t"$4"\t"$5]; next}(($1"\t"$2"\t"$4"\t"$5) in A)' ./strelka/strelka.snv.filtered.vcf ./mutect2/mutect2_filtered_last_snv.vcf | egrep -hv "^#"  >> ${sample}_${ROI}_merged_snv.vcf
awk 'NR==FNR {A[$1"\t"$2"\t"$4"\t"$5]; next}(($1"\t"$2"\t"$4"\t"$5) in A)' ./strelka/strelka.snv.filtered.vcf ./varscan2/varscan2.snp.Somatic.pass.vcf | egrep -hv "^#" >> ${sample}_${ROI}_merged_snv.vcf
awk 'NR==FNR {A[$1"\t"$2"\t"$4"\t"$5]; next}(($1"\t"$2"\t"$4"\t"$5) in A)' ./mutect2/mutect2_filtered_last_snv.vcf ./varscan2/varscan2.snp.Somatic.pass.vcf | egrep -hv "^#" >> ${sample}_${ROI}_merged_snv.vcf


awk '!seen[$1 FS $2 FS $4 FS $5]++' ${sample}_${ROI}_merged_snv.vcf > ${sample}_${ROI}_merged_uniq_snv.vcf
rm ${sample}_${ROI}_merged_snv.vcf

echo "##fileformat=VCFv4.2" > ${sample}_${ROI}_merged_indel.vcf
echo "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR" >> ${sample}_${ROI}_merged_indel.vcf
awk 'NR==FNR {A[$1"\t"$2"\t"$4"\t"$5]; next}(($1"\t"$2"\t"$4"\t"$5) in A)' ./strelka/strelka.indel.filtered.vcf ./mutect2/mutect2_filtered_last_indel.vcf | egrep -hv "^#"  >> ${sample}_${ROI}_merged_indel.vcf
awk 'NR==FNR {A[$1"\t"$2"\t"$4"\t"$5]; next}(($1"\t"$2"\t"$4"\t"$5) in A)' ./strelka/strelka.indel.filtered.vcf ./varscan2/varscan2.indel.Somatic.pass.vcf | egrep -hv "^#" >> ${sample}_${ROI}_merged_indel.vcf
awk 'NR==FNR {A[$1"\t"$2"\t"$4"\t"$5]; next}(($1"\t"$2"\t"$4"\t"$5) in A)' ./mutect2/mutect2_filtered_last_indel.vcf ./varscan2/varscan2.indel.Somatic.pass.vcf | egrep -hv "^#" >> ${sample}_${ROI}_merged_indel.vcf


awk '!seen[$1 FS $2 FS $4 FS $5]++' ${sample}_${ROI}_merged_indel.vcf > ${sample}_${ROI}_merged_uniq_indel.vcf
rm ${sample}_${ROI}_merged_indel.vcf


cp ${sample}_${ROI}_merged_uniq_snv.vcf  ${sample}_${ROI}_merged_uniq_all.vcf
cat ${sample}_${ROI}_merged_uniq_indel.vcf | egrep -hv "^#" >> ${sample}_${ROI}_merged_uniq_all.vcf

awk -v RS='\r?\n' '
NR==FNR {
    chrom=$2;
    start=$3;
    end=$4;
    ranges[chrom, ++count[chrom]] = start FS end;
    next;
}

{
    if ($0 ~ /^#/) {
        print $0;
        next;
    }

    chrom=$1;
    pos=$2;

    is_in_range = 0;
    for (i = 1; i <= count[chrom]; i++) {
        split(ranges[chrom, i], range, FS);
        if (pos >= range[1] && pos <= range[2]) {
            is_in_range = 1;
            break;
        }
    }

    if (!is_in_range) {
        print $0;
    }
}
' /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/UCSC/microsat.txt ${sample}_${ROI}_merged_uniq_all.vcf > ${sample}_${ROI}_merged_uniq_all_filtered1.vcf

awk -v RS='\r?\n' '
NR==FNR {
    chrom=$2;
    start=$3;
    end=$4;
    ranges[chrom, ++count[chrom]] = start FS end;
    next;
}

{
    if ($0 ~ /^#/) {
        print $0;
        next;
    }

    chrom=$1;
    pos=$2;

    is_in_range = 0;
    for (i = 1; i <= count[chrom]; i++) {
        split(ranges[chrom, i], range, FS);
        if (pos >= range[1] && pos <= range[2]) {
            is_in_range = 1;
            break;
        }
    }

    if (!is_in_range) {
        print $0;
    }
}' /dssg/home/acct-guowei_med/anzhijie/pdac3d/references/UCSC/simpleRepeat.txt ${sample}_${ROI}_merged_uniq_all_filtered1.vcf > ${sample}_${ROI}_merged_uniq_all_filtered2.vcf

/dssg/home/acct-guowei_med/anzhijie/pdac3d/software/annovar/convert2annovar.pl -format vcf4old ${sample}_${ROI}_merged_uniq_all_filtered2.vcf > ${sample}_${ROI}.avinput
/dssg/home/acct-guowei_med/anzhijie/pdac3d/software/annovar/table_annovar.pl ${sample}_${ROI}.avinput /dssg/home/acct-guowei_med/anzhijie/pdac3d/software/annovar/humandb/ -buildver hg38 -out ${sample}_${ROI} -remove -protocol refGene,knownGene,gnomad_genome,cosmic70,clinvar_20250721 -operation g,g,f,f,f -nastring . -polish
