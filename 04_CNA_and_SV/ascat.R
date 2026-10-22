library(ASCAT)
library(tidyverse)

args<-commandArgs(trailingOnly = T)
sampleID = args[1]
roi = args[2]
gd = args[3]

print("before setting wd")
getwd()
setwd("./ascat/")
print("after setting wd:")
getwd()
##### extract 
print("sample is:")
print(sampleID)

print("roi is:")
print(roi)

print(paste0("../",sampleID,"_",roi,".rmdup.bqsr.bam"))
print(paste0("../../NOR_ref/",sampleID, "_NOR_ref.rmdup.bqsr.bam"))
print(gd)
print(gd %in% c("XX","XY"))

# extract logR and BAF from high throughput sequencing (HTS, WGS) data and running ASCAT
ascat.prepareHTS(
  tumourseqfile = paste0("/dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/", sampleID, "/", roi, "/", sampleID, "_", roi, ".rmdup.bqsr.bam"),
  normalseqfile = paste0("/dssg/home/acct-guowei_med/anzhijie/pdac3d/wgs/", sampleID, "/NOR_ref/", sampleID, "_NOR_ref.rmdup.bqsr.bam"),
  tumourname = roi,
  normalname = "NOR_ref",
  allelecounter_exe = "~/anaconda3/envs/wgs/bin/alleleCounter",
  alleles.prefix = "/dssg/home/acct-guowei_med/anzhijie/pdac3d/references/ascat/G1000_alleles/G1000_alleles_hg38_chr",
  loci.prefix = "/dssg/home/acct-guowei_med/anzhijie/pdac3d/references/ascat/G1000_loci/G1000_loci_hg38_chr",
  gender = gd,
  genomeVersion = "hg38",
  nthreads = 24,
  tumourLogR_file = "./Tumor_LogR.txt",
  tumourBAF_file = "./Tumor_BAF.txt",
  normalLogR_file = "./Germline_LogR.txt",
  normalBAF_file = "./Germline_BAF.txt")

ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", Germline_LogR_file = "Germline_LogR.txt", Germline_BAF_file = "Germline_BAF.txt", gender = "XY", genomeVersion = "hg38")
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "/dssg/home/acct-guowei_med/anzhijie/pdac3d/references/ascat/GC_G1000/GC_G1000_hg38.txt", replictimingfile = "/dssg/home/acct-guowei_med/anzhijie/pdac3d/references/ascat/RT_G1000/RT_G1000_hg38.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = TRUE)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')
