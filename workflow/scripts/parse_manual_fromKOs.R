parse_manual_fromKOs = function(KEGG_KO_tab, 
                                KEGG_manual_tab) {
  #! library
  suppressMessages(suppressWarnings(library(tidyverse)))
  
  #! import raw data
  ko_list = read.delim(KEGG_KO_tab)
  ko_list = ko_list[, c("filename","locus_tag","gene_ann")]
  
  #! create container for outputs
  trait_manual = data.frame(matrix(NA, 0, 4))
  colnames(trait_manual) = c("filename","locus_tag","gene_ann","trait")
  trait_manual = list(trait_manual)
  
  
  ##-- add DHPS annotation 
  #! seek for the key KO:
  #! "K15509" (sulfopropanediol 3-dehydrogenase, gene hpsN)
  dhps_ann1 = ko_list[ko_list$gene_ann %in% "K15509", ]
  if (nrow(dhps_ann1) > 0) {
    dhps_ann1$trait = "DHPS_catabolism"
    trait_manual = c(trait_manual, list(dhps_ann1))
  }
  ##--
  
  
  ##-- Taurine 
  #! seek for the key KOs:
  
  #! 1) "K07255" + "K07256" (taurine dehydrogenase, tauXY genes) & K03852 (sulfoacetaldehyde acetyltransferase, xsc gene)
  gene_set_A = c("K07255","K07256","K03852")
  tau_ann_A = ko_list %>% 
    filter(gene_ann %in% gene_set_A)
  if (nrow(tau_ann_A) > 0) {
    tau_ann_A_compl = tau_ann_A %>% 
      mutate(value = 1) %>%
      group_by(filename, gene_ann) %>%
      summarise(value = sum(value)) %>%
      pivot_wider(id_cols = filename, names_from = gene_ann, values_from = value, values_fill = 0) %>%
      column_to_rownames("filename")
    #! use completeness rule as in 10.1038/s42003-022-03184-4
    tau_ann_A_compl = rownames(tau_ann_A_compl)[apply(tau_ann_A_compl, 1, function(i) sum(i > 0)) > 2] 
    tau_ann_A$trait = NA
    tau_ann_A$trait[tau_ann_A$filename %in% tau_ann_A_compl & tau_ann_A$gene_ann %in% gene_set_A] = "Taurine_in_TCA"
    trait_manual = c(trait_manual, list(tau_ann_A))
  }
  
  #! 2) "K03851" (taurine-pyruvate aminotransferase, tpa gene) & K03852 (sulfoacetaldehyde acetyltransferase, xsc gene)
  gene_set_B = c("K03851","K03852")
  tau_ann_B = ko_list %>% 
    filter(gene_ann %in% gene_set_B)
  if (nrow(tau_ann_B) > 0) {
    tau_ann_B_compl = tau_ann_B %>% 
      mutate(value = 1) %>%
      group_by(filename, gene_ann) %>%
      summarise(value = sum(value)) %>%
      pivot_wider(id_cols = filename, names_from = gene_ann, values_from = value, values_fill = 0) %>%
      column_to_rownames("filename")
    #! use completeness rule as in 10.1038/s42003-022-03184-4
    tau_ann_B_compl = rownames(tau_ann_B_compl)[apply(tau_ann_B_compl, 1, function(i) sum(i > 0)) > 2] 
    tau_ann_B$trait = NA
    tau_ann_B$trait[tau_ann_B$filename %in% tau_ann_B_compl & tau_ann_B$gene_ann %in% gene_set_B] = "Taurine_in_TCA"
    trait_manual = c(trait_manual, list(tau_ann_B))
  }
  ##--
  
  
  #! merge tables
  trait_manual = do.call(rbind, trait_manual)
  trait_manual = trait_manual %>% 
    filter(!is.na(trait)) %>%
    group_by(across(c(-trait))) %>%
    summarise(trait = paste0(trait, collapse = ";"))
  
  #! save output
  write.table(trait_manual, file = KEGG_manual_tab, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}


parse_manual_fromKOs(KEGG_KO_tab = snakemake@input[["KEGG_KO_tab"]],
                     KEGG_manual_tab = snakemake@output[["KEGG_manual_tab"]]
)