#! wrapper function
KM_reconstruction_wrapper = function(KEGG_KO_tab, 
                                     KEGG_KM_tab,
                                     KM_reco_script, KM_str, ncore) {
  #! load packages
  suppressMessages(suppressWarnings(library(tidyverse)))
  source(KM_reco_script)
  
  #! check if KEGG module structure file exist, otherwise fetch it
  if (file.exists(KM_str)) {
    KM_structures <- readRDS(KM_str)
  } else {
    # it takes a few minutes
    KMdiagram_fetcher(ncore = ncore, path = dirname(KM_str))
    KM_structures <- readRDS(KM_str)
  }
  
  #! parse KO annotation
  KEGG_KM_tab_long <- read.delim(KEGG_KO_tab, header = T)
  KEGG_KM_tab_wide <- KEGG_KM_tab_long %>% 
    mutate(value = 1) %>% 
    pivot_wider(id_cols = filename, names_from = gene_ann, names_sort = T,
                values_from = value, values_fn = sum, values_fill = 0) %>%
    column_to_rownames("filename")
  
  #! run KEGG module reconstruction
  #! gap-filling params from https://doi.org/10.1038/s42003-022-03184-4
  KMreco <- KMreco(indata = KEGG_KM_tab_wide, KM_str = KM_structures,  
                  len_breaks = c(3), allowed_gaps = c(0,1), ncore = ncore)
  
  #! filter KMs completeness (!'complete_perc' includes already gap-filling bonus)
  KMreco_filt <- KMreco[KMreco$complete_perc >= 1, ]
  
  #! parse output
  KMreco_filt2 <- KMreco_filt %>% group_by(across(c(-trait))) %>%
    summarise(trait = paste0(unique(trait), collapse = ";"))
  KMreco_filt2$KM_long <- sapply(KMreco_filt2$trait, function(i) {
    i2 <- KM_structures[names(KM_structures) %in% unlist(strsplit(i, ";"))]
    i3 <- paste0(sapply(i2, "[[", "KM_desc"), collapse = ";")
    return(i3)
  })

  
  #! save output
  write.table(KMreco_filt2, file = KEGG_KM_tab, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}

KM_reconstruction_wrapper(KEGG_KO_tab = snakemake@input[["KEGG_KO_tab"]],
                          KEGG_KM_tab = snakemake@output[["KEGG_KM_tab"]],
                          KM_reco_script = snakemake@params[["KM_reco_script"]],
                          KM_str = snakemake@params[["KM_str"]],
                          ncore = snakemake@params[["ncore"]]
)