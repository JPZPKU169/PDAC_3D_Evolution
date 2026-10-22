from SigProfilerExtractor import sigpro as sig

if __name__ == '__main__':
    # This probably means that you are not using fork to start your child processes and you have forgotten to use the proper idiom
    sig.sigProfilerExtractor(input_type="matrix", 
                             output="results_SBS96_max10_filtered4", 
                             input_data="/io01dat03/PROJECT/zoujiping/pdac/wgs/sigprofiler/vcfs_filtered4/output/SBS/PDAC.SBS96.all", 
                             reference_genome="GRCh38", 
                             context_type="96",
                             exome=False,
                             minimum_signatures=2, 
                             maximum_signatures=10, 
                             nmf_replicates=800, cpu=16)
