# #___ for testing
# setwd("~/VSC5_data/metags/B03_M0003_aquatic_MAGs/")
# mapping_dir <- "2_assembly/mapping_assembly"
# out_repo <- "99_reports/Mapping_contig.report"
# mapping_dir <- "2_assembly/mapping_mags"
# out_repo <- "99_reports/Mapping_mag.report"
# #___ for testing



mapping_report <- function(out_repo, 
                           mapping_dir, use_filtered="TRUE") {
  
  #! load packages & functions
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressMessages(suppressWarnings(library(Biostrings)))
  
  
  #! import data
  mapping_list <- list.files(mapping_dir, pattern = ".stats$", full.names = T, recursive = T)
  if (use_filtered == "TRUE") {
    mapping_list <- mapping_list[grepl("_filt.stats", mapping_list)]
  } else {
    mapping_list <- mapping_list[!grepl("_filt.stats", mapping_list)]
  }
  mapping_stats <- lapply(mapping_list, function(i) read.delim(i, h=F))
  names(mapping_stats) <- gsub(".stats", "", basename(mapping_list))
  
  
  #! parse data
  mapping_repo <- do.call(rbind, mapping_stats)
  mapping_repo <- data.frame(sample=gsub("\\.[0-9]+$", "", row.names(mapping_repo)), mapping_repo)
  mapping_repo <- mapping_repo %>%
    group_by(sample) %>%
    summarise(mapped=sum(V3), 
              unmapped=sum(V4),
              n_contigs=sum(V2 > 0),
              N50=N50(V2[V2 > 0]),
              longest=max(V2[V2 > 0]),
              mean_len=round(mean(V2[V2 > 0]), 0)
    ) %>%
    rowwise() %>% mutate(map_eff=mapped/sum(mapped+unmapped)) %>% 
    relocate(map_eff, .after = unmapped)
  
  
  #! outputs
  write.table(mapping_repo, file = out_repo, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}



mapping_report(out_repo = snakemake@output[["out_repo"]],
               mapping_dir = snakemake@params[["mapping_dir"]],
               use_filtered = snakemake@params[["use_filtered"]]
)