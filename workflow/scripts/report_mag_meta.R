# #___ for testing
# setwd("~/VSC_data/metags/B03_M0001_Chicken_gut_microbiome/")
# mag_dir <- "3_MAG_analysis/drep_fasta/"
# checkm_dir <- "3_MAG_analysis/checkm/"
# gtdbtk_dir <- "3_MAG_analysis/gtdbtk/"
# out_tab <- "99_reports/"
# GFF_table <- "4_MAG_annotation/prokka_tab.tsv"
# #___ for testing



mags_report <- function(mag_dir, checkm_dir, gtdbtk_dir,
                        out_tab,
                        GFF_table=NULL) {
  
  #! load packages & functions
  suppressMessages(suppressWarnings(library(Biostrings)))
  
  
  #! import data
  mag_list <- list.files(mag_dir, pattern = ".fa$", full.names = T, recursive = T)
  mag_files <- lapply(mag_list, function(i) readDNAStringSet(i))
  names(mag_files) <- basename(mag_list)
  mag_checkm <- read.delim(file.path(checkm_dir, "results.tsv"), h=T)
  mag_gtdbtk <- read.delim(file.path(gtdbtk_dir, "gtdbtk.bac120.summary.tsv"), h=T)
  
  
  #! compute metrics
  mag_meta <- data.frame(MAG=gsub(".fa", "", mag_checkm$MAG_ids), mag_checkm)
  mag_meta$Taxonomy <- mag_gtdbtk$classification[match(mag_meta$MAG, mag_gtdbtk$user_genome)]
  mag_meta$len <- sapply(mag_meta$MAG_ids, function(i) sum(mag_files[[i]]@ranges@width))
  mag_meta$n_contigs <- sapply(mag_meta$MAG_ids, function(i) length(mag_files[[i]]))
  mag_meta$N50 <- sapply(mag_meta$MAG_ids, function(i) N50(mag_files[[i]]@ranges@width))
  mag_meta$longest <- sapply(mag_meta$MAG_ids, function(i) max(mag_files[[i]]@ranges@width))
  mag_meta$mean_len <- sapply(mag_meta$MAG_ids, function(i) round(mean(mag_files[[i]]@ranges@width), 0))
  
  if (file.exists(GFF_table)) {
    suppressMessages(suppressWarnings(library(tidyverse)))
    
    GFF_ann <- read.delim(GFF_table) %>%
      group_by(filename, type) %>%
      summarise(n_gene=n())
    
    mag_meta$ngene <- sapply(mag_meta$MAG, function(i) sum(GFF_ann$n_gene[GFF_ann$filename == i]))
    mag_meta$tRNA <- sapply(mag_meta$MAG, function(i) sum(GFF_ann$n_gene[GFF_ann$filename == i &
                                                                           GFF_ann$type == "tRNA"]))
    mag_meta$rRNA <- sapply(mag_meta$MAG, function(i) sum(GFF_ann$n_gene[GFF_ann$filename == i &
                                                                           GFF_ann$type == "rRNA"]))
  }
  
  
  #! outputs
  write.table(mag_meta, file = out_tab, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}



mags_report(mag_dir = snakemake@input[["mag_dir"]],
            checkm_dir = snakemake@input[["checkm_dir"]],
            gtdbtk_dir = snakemake@input[["gtdbtk_dir"]],
            GFF_table = snakemake@input[["GFF_table"]],
            out_tab = snakemake@output[["out_tab"]]
)