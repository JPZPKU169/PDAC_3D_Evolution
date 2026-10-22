
##### Meskit bootstrap 1000 #####
case_all = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10")

library(MesKit)
for(spID in case_all){
  print(spID)
  library(MesKit)
  roi_info = read.csv(paste0("E:/Projects/PDAC/results/WGS/Tree_tune/", spID, "_ROIs_info.csv"), row.names = 1)
  rois = roi_info$rois
  
  
  maf = read.table(file = paste0("E:/Projects/PDAC/results/WGS/", spID, "/", spID, "_mutation_list.txt"), sep = "\t")
  readcount = read.table(file = paste0("E:/Projects/PDAC/results/WGS/", spID, "/", spID, "_readcount.txt"), check.names = FALSE)
  maf = filter(maf, sample %in% rois)
  
  # 筛选出KRAS突变ROI
  kras_exonic = filter(maf, Gene.refGene == "KRAS" & Func.refGene == "exonic")
  kras_exonic = kras_exonic$sample
  
  
  library(MesKit)
  print("Construct phylogenic tree using rescued mutation:")
  if(0.2 * length(unique(maf_unified$sample)) > 2){
    min_check_count = round(0.2 * length(unique(maf_unified$sample)))
  }else{
    min_check_count = 2
  }
  
  maf_rescued = supplement_maf(maf, readcount_unified, min_check_samplecount = min_check_count, 
                               min_alt_reads = 3, min_vaf = 0.05)
  
  maf_meskit_rescued = maf_rescued
  maf_meskit_rescued = filter(maf_meskit_rescued, Variant_Type == "SNP")
  maf_meskit_rescued = filter(maf_meskit_rescued, sample %in% rois)
  maf_meskit_rescued = maf_meskit_rescued[,c("case","Hugo_Symbol","Chr","Start","End","Variant_Classification","Variant_Type","Ref",
                                             "Alt","s_Ref","s_Alt","s_VAF","sample")]
  colnames(maf_meskit_rescued) = c("Patient_ID","Hugo_Symbol","Chromosome", "Start_Position", "End_Position", "Variant_Classification",
                                   "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2", "Ref_allele_depth",
                                   "Alt_allele_depth",   "VAF", "Tumor_Sample_Barcode")
  maf_meskit_rescued = unique(maf_meskit_rescued)
  write.table(maf_meskit_rescued, file = paste0("E:/Projects/PDAC/results/WGS/Tree_tune/Meskit_bootstrap1000/", spID, "_maf_rescued_meskit.maf"), 
              sep = "\t", quote = F, row.names = F)
  
  meskit_obj = readMaf(mafFile = paste0("E:/Projects/PDAC/results/WGS/Tree_tune/Meskit_bootstrap1000/", spID, "_maf_rescued_meskit.maf"),
                       clinicalFile = paste0("E:/Projects/PDAC/results/WGS/Tree_tune/Meskit_bootstrap1000/", spID, "_meskit_clinical.txt"),
                       refBuild = 'hg38')
  
  meth = "MP"
  print(meth)
  library(MesKit)
  phylo_res = getPhyloTree(meskit_obj, method = meth, bootstrap.rep.num = 1000)
  saveRDS(phylo_res, file = paste0("E:/Projects/PDAC/results/WGS/Tree_tune/Meskit_bootstrap1000/rescued_", meth, "/", spID, "_rescued_", meth,
                                   "_phylo_res.rds"))
  detach("package:MesKit", unload = TRUE)
  library(ape)
  tree = phylo_res@tree
  rooted_tree <- root(tree, outgroup = "NORMAL", resolve.root = TRUE)
  rooted_tree$node.label = phylo_res@bootstrap.value
  
  
  library(ggtree)
  ggtree(rooted_tree, 
         layout = "rectangular", 
         size = 0.8,
         color = "steelblue") +
    
    coord_flip(clip = "off") + 
    scale_x_reverse() +
    
    geom_tiplab(aes(color = label %in% kras_exonic), angle = 90, hjust = 1, offset = 0.05, size = 3) +
    scale_color_manual(values = c("FALSE"="black","TRUE" = "red")) +
    
    geom_nodelab(aes(subset = !isTip, label = label), 
                 size = 3,
                 color = "black",
                 vjust = -0.5,
                 hjust = -0.2) +
    
    theme_tree() +
    theme(plot.margin = margin(10, 10, 30, 10),
          legend.position = "none")
  ggsave(filename = paste0("E:/Projects/PDAC/results/WGS/Tree_tune/Meskit_bootstrap1000/rescued_", meth, "/", spID, "_minAlt3_rescued_minVAF0.05_", meth,
                           "_tree.png"), width = length(unique(maf_unified$sample)) * 0.6, height = 12, units = "cm", dpi = 300)
  

}


supplement_maf <- function(maf_unified, readcount_unified, min_check_samplecount, min_alt_reads = 3, min_vaf = 0.05) {
  require(data.table)
  require(dplyr)
  
  maf_unified_dt = as.data.table(maf_unified)
  readcount_unified_dt = as.data.table(readcount_unified)

  mutation_counts <- maf_unified_dt[, .(sample_count = uniqueN(sample)), by = mutation_id]
  
  frequent_mutations <- mutation_counts[sample_count > min_check_samplecount, ]
  frequent_mutations <- frequent_mutations[sample_count < length(unique(maf_unified_dt$sample)), mutation_id] # 不考虑在所有样本都存在的突变

  maf_frequent <- maf_unified_dt[mutation_id %in% frequent_mutations]

  all_samples <- unique(maf_unified_dt$sample)

  supplemental_rows <- list()

  for (mut in frequent_mutations) {
    mutated_samples <- maf_frequent[mutation_id == mut, sample]
    
    non_mutated_samples <- setdiff(all_samples, mutated_samples)
    
    ref_cols <- paste0(non_mutated_samples, "_Ref")
    ref_values <- readcount_unified_dt[mutation_id == mut, ..ref_cols]
    alt_cols <- paste0(non_mutated_samples, "_Alt")
    alt_values <- readcount_unified_dt[mutation_id == mut, ..alt_cols]
    vaf_cols = paste0(non_mutated_samples, "_VAF")
    vaf_values = readcount_unified_dt[mutation_id == mut, ..vaf_cols]
    
    alt_df <- data.table(
      sample = non_mutated_samples,
      ref_value = as.numeric(ref_values[1, ]),
      alt_value = as.numeric(alt_values[1, ]),  # 假设每个突变在readcount中只有一行
      vaf_value = as.numeric(vaf_values[1, ])
    )
    
    alt_df <- alt_df[alt_value >= min_alt_reads & vaf_value >= min_vaf]
    
    if (nrow(alt_df) > 0) {
      template_row <- maf_frequent[mutation_id == mut][1]
      
      new_rows <- lapply(alt_df$sample, function(s) {
        new_row <- copy(template_row)
        new_row[, sample := s]
        new_row[, s_Ref := alt_df$ref_value[match(s, alt_df$sample)]]
        new_row[, s_Alt := alt_df$alt_value[match(s, alt_df$sample)]]
        new_row[, s_VAF := alt_df$vaf_value[match(s, alt_df$sample)]]
        return(new_row)
      })
      
      supplemental_rows <- c(supplemental_rows, new_rows)
    }
  }
  
  if (length(supplemental_rows) > 0) {
    supplemental_maf <- rbindlist(supplemental_rows)
    result_maf <- rbindlist(list(maf_unified_dt, supplemental_maf), use.names = TRUE)
  } else {
    result_maf <- maf_unified_dt
  }
  
  setDF(result_maf)
  
  return(result_maf)
}
