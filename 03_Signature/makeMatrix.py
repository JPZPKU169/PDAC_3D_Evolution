import argparse
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

#parser = argparse.ArgumentParser(description='make SBS matrix')
#parser.add_argument("--sample", type=str)
#args = parser.parse_args()
#sample=args.sample
vcf_path = "/io01dat03/PROJECT/zoujiping/pdac/wgs/sigprofiler/vcfs_filtered4"
sample = "PDAC"

matrices = matGen.SigProfilerMatrixGeneratorFunc(sample, "GRCh38", vcf_path, plot=True, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
