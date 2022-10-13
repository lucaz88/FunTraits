# #___ for testing
# setwd("~/VSC_data/lucaz/metags/B03_M0001_Chicken_gut_microbiome")
# MASTER_table <- "4_MAG_annotation/MASTER_table.tsv"
# gtdbtk_dir <- "3_MAG_analysis/gtdbtk"
# ann_modules <- c("prokka","KEGG_KO","KEGG_KM","KEGG_manual","antiSMASH","BioV_transp","blast_phytohormones","blast_vibrioferrin","blast_DMSP","dbCAN_CAZy")
# # ann_modules <- ann_modules[-c(7,8)]
# min_trait_occur <- 3
# dist_mt <- "jaccard"
# aggl_mt <- "ward.D2"
# taxa_col <- "~/VSC_data/lucaz/_DBs/MY_taxa_cols.tsv"
# outfile <- "4_MAG_annotation/hm_MASTERtraits_jacc.html"
# 
# plot_FunLuca_hm(MASTER_table, gtdbtk_dir, ann_modules, dist_mt, aggl_mt, outfile)
# #___ for testing



plot_FunLuca_hm <- function(MASTER_table, gtdbtk_dir, 
                            outfile,
                            ann_modules, min_trait_occur, dist_mt, aggl_mt, taxa_col=NULL) {
  #! libraries
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressMessages(suppressWarnings(library(vegan)))
  suppressMessages(suppressWarnings(library(heatmaply)))
  suppressMessages(suppressWarnings(library(ggplot2)))
  
  
  #! parse inputs
  ann_modules_filt <- ann_modules[!ann_modules %in% c("prokka","KEGG_KO")] # remove ann that do not produce traits
  hm_matrix <- read.delim(MASTER_table) %>%
    separate_rows(matches(paste0(ann_modules_filt, collapse = "|")), sep = ";") %>% # some traits shared same gene/s
    pivot_longer(cols = matches(paste0(ann_modules_filt, collapse = "|")), 
                 names_to = "ann_modules", values_to = "trait") %>%
    filter(!is.na(trait)) %>%
    mutate(value = 1) %>%
    pivot_wider(id_cols = filename, 
                names_from = trait, names_sort = T,
                values_from = value, values_fn = sum, values_fill = 0) %>%
    column_to_rownames("filename")
  gtdb_summary <- read.delim(file.path(gtdbtk_dir, "gtdbtk.bac120.summary.tsv"))
  
  
  #! filter
  hm_matrix <- hm_matrix[, apply(hm_matrix, 2, function(i) sum(i > 0)) >= min_trait_occur] # a trait must be present in at least 3 genomes
  
  
  #! tweaks for testing
  hm_matrix[hm_matrix > 1] <- 1
  
  
  #! row annotations (from GTDBtk)
  row_ann <- gtdb_summary  %>%
    mutate(Taxa = sapply(lapply(strsplit(classification, ";"), 
                                function(i) gsub("[a-z]__", "", i)), 
                         function(j) ifelse(any(grepl("Proteobacteria|Firmicutes", j)), j[3], j[2]))) %>%
    right_join(data.frame(user_genome=row.names(hm_matrix))) %>%
    replace_na(list(Taxa = "unknown")) %>%
    arrange(match(user_genome, row.names(hm_matrix))) %>%
    select(Taxa)
  row_ann$Taxa[row_ann$Taxa %in% names(table(row_ann$Taxa))[table(row_ann$Taxa) < 5]] = "Others"
  if (!is.null(taxa_col)) {
    taxa_col <- read.delim(taxa_col, h=T)
    taxa_col <- taxa_col[taxa_col$Taxa %in% row_ann$Taxa, ]
    taxa_col <- structure(taxa_col$HEX, names=taxa_col$Taxa)
    if (!all(unique(row_ann$Taxa) %in% names(taxa_col))) { # reset color in case some taxa is not mapped in the provided palette
      taxa_col <- NULL
    }
  }
  
  #! col annotations
  col_ann <- read.delim(MASTER_table) %>%
    separate_rows(matches(paste0(ann_modules_filt, collapse = "|")), sep = ";") %>% # some traits shared genes
    pivot_longer(cols = matches(paste0(ann_modules_filt, collapse = "|")), 
                 names_to = "ann_modules", values_to = "trait") %>%
    filter(!is.na(trait))
  col_ann <- col_ann[match(colnames(hm_matrix), col_ann$trait), "ann_modules", drop=F]
  
  
  #! cell cols
  cell_col <- scale_fill_gradientn(
    colours = c("#f6f2f0", "#330066"),
    values = c(0, 1)
  )
  
  
  #! cell extra text
  over_text <- matrix(data = "", nrow = nrow(hm_matrix), ncol = ncol(hm_matrix))
  over_text <- apply(over_text, 2, function(i) {
    paste0(i, gtdb_summary$classification[match(rownames(hm_matrix), gtdb_summary$user_genome)])
  })
  
  
  #! calculate dendrograms
  row_hc <- hclust(vegdist(hm_matrix, dist_mt, binary = T), method = aggl_mt)
  col_hc <- hclust(vegdist(t(hm_matrix), dist_mt, binary = T), method = aggl_mt)
  
  
  #! plotting & saving
  dir.create(dirname(outfile), recursive = T, showWarnings = F)
  heatmaply(as.matrix(hm_matrix), scale = "none", 
            xlab = "Genomic traits", ylab = "Genomes",
            showticklabels = c(F, F), # remove x and y labels, essential with large matrices
            margins = c(70,70,20,20),
            scale_fill_gradient_fun = cell_col,
            Rowv = row_hc, Colv = col_hc,
            col_side_colors = col_ann,
            # col_side_palette = smpl_col,
            row_side_colors = row_ann,
            row_side_palette = taxa_col,
            # colorbar_xanchor='left',colorbar_yanchor='top', colorbar_xpos=1.1, colorbar_ypos=-0.1,
            custom_hovertext = over_text,
            file = outfile
  )
  
  
}



plot_FunLuca_hm(MASTER_table = snakemake@input[["MASTER_table"]],
                gtdbtk_dir = snakemake@input[["gtdbtk_dir"]],
                outfile = snakemake@output[["outfile"]],
                ann_modules = snakemake@params[["ann_modules"]],
                min_trait_occur = snakemake@params[["min_trait_occur"]],
                dist_mt = snakemake@params[["dist_mt"]],
                aggl_mt = snakemake@params[["aggl_mt"]],
                taxa_col = snakemake@params[["taxa_col"]]
)