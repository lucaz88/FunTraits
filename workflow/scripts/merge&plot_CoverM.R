# #___ for testing
# setwd("~/VSC_data/lucaz/metags/B03_M0001_Chicken_gut_microbiome/")
# coverm_out <- c("_MAG_analysis/coverm/D14A--D14A_coverm","_MAG_analysis/coverm/D21A--D21A_coverm",
#                "_MAG_analysis/coverm/D35B--D35B_coverm","_MAG_analysis/coverm/D35C--D35C_coverm",
#                "_MAG_analysis/coverm/POOL--POOL_coverm")
# gtdbtk_dir <- "_MAG_analysis/gtdbtk"
# merged_tpm <- "merged_TPM.tsv"
# merged_cover <- "merged_cov.tsv"
# hm_tpm <- "hm_TPM.html"
# hm_cover <- "hm_cov.html"
# mapper_meth <- "mean relative_abundance covered_fraction count rpkm tpm"
# taxa_col <- "../../_DBs/MY_taxa_cols.tsv"
# dist_mt <- "bray"
# aggl_mt <- "ward.D2"
# #___ for testing



merge_CoverM <- function(coverm_out, gtdbtk_dir, 
                        merged_tpm, merged_cover, hm_tpm, hm_cover,
                        mapper_meth, taxa_col=NULL, smpl_col=NULL, dist_mt="bray", aggl_mt="ward.D2") {
  
  
  #! load packages
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressMessages(suppressWarnings(library(heatmaply)))
  suppressMessages(suppressWarnings(library(vegan)))
  suppressMessages(suppressWarnings(library(ggplot2)))
  suppressMessages(suppressWarnings(library(scales)))
  
  
  #! set variables
  cols <- c("Genome", strsplit(mapper_meth, " ")[[1]])
  gtdb_summary <- file.path(gtdbtk_dir, "gtdbtk.bac120.summary.tsv")
  
  
  #! import data
  covm_files <- lapply(coverm_out, function(i) read.delim(file.path(i, "MAGs_mapping.txt"), h=T, col.names = cols))
  names(covm_files) <- gsub(".txt", "", basename(coverm_out))
  covm_files <- do.call(rbind, covm_files)
  covm_files <- cbind(Sample=gsub("\\.[0-9]+$", "", rownames(covm_files)), covm_files)
  gtdb_tax <- read.delim(gtdb_summary, h=T) %>%
    mutate(Taxa = sapply(lapply(strsplit(classification, ";"), 
                                function(i) gsub("[a-z]__", "", i)), 
                         function(j) ifelse(any(grepl("Proteobacteria|Firmicutes", j)), j[3], j[2]))) %>%
    select(user_genome, classification, Taxa) %>%
    replace_na(list(Taxa = "unknown"))

  
  
  #! parse data
  # covm_files <- covm_files %>%
  #   group_by(Sample) %>% 
  #   mutate(count2 = sum(count, na.rm = T)/sum(relative_abundance[Genome != "unmapped"])*sum(relative_abundance[Genome == "unmapped"])) # try to infer number unmapped reads

  covm_tpm <- covm_files %>%
    filter(Genome != "unmapped") %>%
    pivot_wider(id_cols = Genome, names_from = Sample, values_from = tpm, values_fill = 0)
  covm_tpm <- merge(covm_tpm, gtdb_tax, by.x = "Genome", by.y = "user_genome", all.x = T)
  
  covm_cover <- covm_files %>%
    filter(Genome != "unmapped") %>%
    pivot_wider(id_cols = Genome, names_from = Sample, values_from = covered_fraction, values_fill = 0)
  covm_cover <- merge(covm_cover, gtdb_tax, by.x = "Genome", by.y = "user_genome", all.x = T)
  
  if (!is.null(taxa_col)) {
    taxa_col <- read.delim(taxa_col, h=T)
    taxa_col <- taxa_col[taxa_col$Taxa %in% covm_tpm$Taxa, ]
    taxa_col <- structure(taxa_col$HEX, names=taxa_col$Taxa)
  }
  # if (!is.null(smpl_col)) {
  #   smpl_col <- read.delim(smpl_col, h=T)
  #   smpl_col <- smpl_col[smpl_col$Sample %in% covm_tpm$Sample, ]
  #   smpl_col <- structure(smpl_col$HEX, names=smpl_col$Sample)
  # }
  
  
  #! plot output
  cell_col <- function(x) {
    scale_fill_gradientn(
      colours = c("white", "#330066", "#e60b0b"),
      values = c(0, rescale(mean(unlist(x), na.rm = T), to = c(0,1), from = c(min(unlist(x), na.rm = T), max(unlist(x), na.rm = T))), 1)
    )}
  
  mat_tpm <- covm_tpm %>% 
    column_to_rownames("Genome") %>%
    select(-c(classification, Taxa)) #%>%
    # log() %>% 
    # mutate(across(everything(), ~replace(., . ==  0 , NA)))
  over_text <- matrix(data = "", nrow = nrow(mat_tpm), ncol = ncol(mat_tpm))
  over_text <- apply(mat_tpm, 2, function(i) {
    paste0(i, gtdb_tax$classification[match(rownames(mat_tpm), gtdb_tax$user_genome)])
  })
  row_hc = hclust(vegdist(mat_tpm, dist_mt, na.rm = T), method = aggl_mt)
  col_hc = hclust(vegdist(t(mat_tpm), dist_mt, na.rm = T), method = aggl_mt)
  
  heatmaply(mat_tpm,
            scale = "none",
            scale_fill_gradient_fun = cell_col(mat_tpm),
            Rowv = row_hc, Colv = col_hc,
            # showticklabels = c(F, F), # remove x and y labels, essential with large matrices
            # col_side_colors = data.frame(Taxa=covm_tpm$Taxa),
            # col_side_palette = smpl_col,
            row_side_colors = data.frame(Taxa=covm_tpm$Taxa),
            row_side_palette = taxa_col,
            custom_hovertext = over_text,
            label_names = c("x", "y", "TPM"),
            file = hm_tpm
            )
  
  mat_cover <- covm_cover %>% 
    column_to_rownames("Genome") %>%
    select(-c(classification, Taxa))
  heatmaply(mat_cover,
            scale = "none",
            scale_fill_gradient_fun = cell_col(mat_cover),
            Rowv = row_hc, Colv = col_hc,
            # showticklabels = c(F, F), # remove x and y labels, essential with large matrices
            # col_side_colors = data.frame(Taxa=covm_tpm$Taxa),
            # col_side_palette = smpl_col,
            row_side_colors = data.frame(Taxa=covm_cover$Taxa),
            row_side_palette = taxa_col,
            custom_hovertext = over_text,
            label_names = c("x", "y", "MAG coverage"),
            file = hm_cover
  )
  
  
  #! save output
  write.table(covm_tpm, file = merged_tpm, 
              col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(covm_cover, file = merged_cover, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}



merge_CoverM(coverm_out = snakemake@input[["coverm_out"]],
             gtdbtk_dir = snakemake@input[["gtdbtk_dir"]],
             merged_tpm = snakemake@output[["merged_tpm"]],
             merged_cover = snakemake@output[["merged_cover"]],
             hm_tpm = snakemake@output[["hm_tpm"]],
             hm_cover = snakemake@output[["hm_cover"]],
             mapper_meth = snakemake@params[["mapper_meth"]],
             taxa_col = snakemake@params[["taxa_col"]],
             # smpl_col = snakemake@params[["smpl_col"]],
             dist_mt = snakemake@params[["dist_mt"]],
             aggl_mt = snakemake@params[["aggl_mt"]]
)