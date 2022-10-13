#___ for testing
setwd("~/VSC_data/lucaz/metags/B03_M0001_Chicken_gut_microbiome/")
merged_tpm <- "3_MAG_analysis/coverm/merged_TPM.tsv"
MASTER_table <- "4_MAG_annotation/MASTER_table.tsv"
out_dir <- "5_plots"
ann_modules <- c("prokka","KEGG_KO","KEGG_KM","KEGG_manual","antiSMASH","BioV_transp","blast_phytohormones","blast_vibrioferrin","blast_DMSP","dbCAN_CAZy")
taxa_col <- "../../_DBs/MY_taxa_cols.tsv"
min_trait_occur <- 3
min_tpm <- 5
dist_mt <- "bray"
aggl_mt <- "ward.D2"
#___ for testing



merge_FunLuca <- function(merged_tpm, MASTER_table, 
                          out_dir,
                          ann_modules, taxa_col=NULL, smpl_col=NULL, min_trait_occur, 
                          min_tpm, dist_mt="bray", aggl_mt="ward.D2") {
  
  #! load packages & functions
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressMessages(suppressWarnings(library(apcluster)))
  suppressMessages(suppressWarnings(library(vegan)))
  suppressMessages(suppressWarnings(library(heatmaply)))
  suppressMessages(suppressWarnings(library(scales)))
  cell_col <- function(x) {
    scale_fill_gradientn(
      colours = c("white", "#330066", "#e60b0b"),
      values = c(0, rescale(mean(unlist(x), na.rm = T), to = c(0,1), 
                            from = c(min(unlist(x), na.rm = T), max(unlist(x), na.rm = T))), 1)
    )}
  aff_prop_hc <- function(x) {
    # trait_cor <- cor(trait_smpl_mat, method = )
    # trait_cor[trait_cor < 0] <- 0
    trait_cor <- 1-as.matrix(vegdist(t(x), dist_mt, na.rm = T))
    
    apcl_trait = apcluster(s = trait_cor, details=T, q=0.5, lam=0.5, seed=1234, maxits=1000, convits=500)
    trait_cl_tmp = do.call(rbind, lapply(1:length(apcl_trait@clusters), function(i) data.frame(i, apcl_trait@clusters[[i]])))
    trait_cl = trait_cl_tmp$i[order(trait_cl_tmp$apcl_trait.clusters..i.., decreasing = F)]
    
    trait_cl_color = structure(gray.colors(length(unique(trait_cl))), names=unique(trait_cl))
    
    trait_dist = as.matrix(cophenetic(as.dendrogram(aggExCluster(s = trait_cor, x = apcl_trait))))
    trait_dist = trait_dist[match(row.names(trait_cor), row.names(trait_dist)), 
                            match(colnames(trait_cor), colnames(trait_dist))]
    
    trait_hc = hclust(dist(trait_dist), aggl_mt)
    trait_hc = as.hclust(reorder(as.dendrogram(trait_hc),
                                 wts = colSums(trait_cor), agglo.FUN = "mean")) # improve dendro sorting
    trait_hc$order = as.integer(trait_hc$order) # otherwise it rises an issue when plotting with iheatmapr
    
    return(list(hc = trait_hc,
                clusters = trait_cl,
                clusters_col = trait_cl_color))
  }
  
  
  #! import data
  merged_tpm <- read.delim(merged_tpm, h = T)
  ann_modules_filt <- ann_modules[!ann_modules %in% c("prokka","KEGG_KO")] # remove ann that do not produce traits
  FunLuca_mat <- read.delim(MASTER_table) 
  dir.create(out_dir, recursive = T, showWarnings = F)
  
  
  #! parse data
  MAG_ann_mat <- FunLuca_mat %>%
    separate_rows(matches(paste0(ann_modules_filt, collapse = "|")), sep = ";") %>% # some traits shared same gene/s
    pivot_longer(cols = matches(paste0(ann_modules_filt, collapse = "|")), 
                 names_to = "ann_modules", values_to = "trait") %>%
    filter(!is.na(trait)) %>%
    mutate(value = 1) %>%
    group_by(filename, trait) %>%
    summarise(ngene = sum(value))
  MAG_tpm <- merged_tpm %>%
    pivot_longer(cols = ends_with("_coverm"), 
                 names_to = "Samples", values_to = "tpm") %>%
    filter(tpm > 0) 
  MAG_ann_tpm_mat <- merge(MAG_tpm, MAG_ann_mat, by.x = "Genome", by.y = "filename")

  
  #! filtering
  MAG_ann_tpm_mat_filt <- MAG_ann_tpm_mat
  trait_freq <- table(MAG_ann_mat$trait)
  MAG_ann_tpm_mat_filt <- MAG_ann_tpm_mat_filt[MAG_ann_tpm_mat_filt$trait %in% names(trait_freq)[trait_freq > min_trait_occur], ] # filter trait based on frequency
  MAG_ann_tpm_mat_filt <- MAG_ann_tpm_mat_filt[MAG_ann_tpm_mat_filt$trait > min_tpm, ] # filter trait based min TPM counts
  
  trait_smpl_mat <- MAG_ann_tpm_mat_filt %>%
    group_by(Samples, trait) %>%
    summarise(tpm = sum(tpm)) %>% 
    pivot_wider(id_cols = Samples, 
                names_from = trait, names_sort = T,
                values_from = tpm, values_fn = sum, values_fill = NA) %>%
    column_to_rownames("Samples")
  
  taxa.trait_smpl_mat <- MAG_ann_tpm_mat_filt %>%
    unite(taxa.trait, c("Taxa", "trait")) %>%
    group_by(Samples, taxa.trait) %>%
    summarise(tpm = sum(tpm)) %>% 
    pivot_wider(id_cols = Samples, 
                names_from = taxa.trait, names_sort = T,
                values_from = tpm, values_fn = sum, values_fill = NA) %>%
    column_to_rownames("Samples")
  
  
  #! output - Traits VS Samples
  trait_analysis <- aff_prop_hc(trait_smpl_mat)
  trait_hc <- trait_analysis$hc
  smpl_hc <- hclust(vegdist(trait_smpl_mat, dist_mt, na.rm = T), method = aggl_mt)
  
  # if (!is.null(smpl_col)) {
  #   smpl_col <- read.delim(smpl_col, h=T)
  #   smpl_col <- smpl_col[smpl_col$Sample %in% covm_tpm$Sample, ]
  #   smpl_col <- structure(smpl_col$HEX, names=smpl_col$Sample)
  # }
  trait_smpl_hm <- file.path(out_dir, "trait_smpl_hm.html")
  
  heatmaply(trait_smpl_mat,
            scale = "none",
            scale_fill_gradient_fun = cell_col(trait_smpl_mat),
            Rowv = smpl_hc, Colv = trait_hc,
            # showticklabels = c(F, F), # remove x and y labels, essential with large matrices
            # custom_hovertext = over_text,
            col_side_colors = data.frame(Trait_clusters=trait_analysis$clusters),
            col_side_palette = trait_analysis$clusters_col,
            label_names = c("x", "y", "TPM"),
            file = trait_smpl_hm
  )
  write.table(t(cbind(Samples = c("Trait_clusters", row.names(trait_smpl_mat)), 
                    rbind(trait_analysis$clusters, trait_smpl_mat))), 
              file = file.path(out_dir, "trait_smpl_hm.tsv"), 
              col.names = F, row.names = T, quote = F, sep = "\t")
  
  
  #! output - Traits divided per taxa VS Samples
  trait_analysis <- aff_prop_hc(taxa.trait_smpl_mat)
  trait_hc <- trait_analysis$hc
  smpl_hc <- hclust(vegdist(taxa.trait_smpl_mat, dist_mt, na.rm = T), method = aggl_mt)
  
  if (!is.null(taxa_col)) {
    taxa_col_df <- read.delim(taxa_col, h=T)
    taxa_col_df <- taxa_col_df[taxa_col_df$Taxa %in% unique(gsub("_.*", "", colnames(taxa.trait_smpl_mat))), ]
    taxa_col_df <- structure(taxa_col_df$HEX, names=taxa_col_df$Taxa)
  }
  # if (!is.null(smpl_col)) {
  #   smpl_col <- read.delim(smpl_col, h=T)
  #   smpl_col <- smpl_col[smpl_col$Sample %in% covm_tpm$Sample, ]
  #   smpl_col <- structure(smpl_col$HEX, names=smpl_col$Sample)
  # }
  taxa.trait_smpl_hm <- file.path(out_dir, "taxa.trait_smpl_hm.html")
  
  heatmaply(taxa.trait_smpl_mat,
            scale = "none",
            scale_fill_gradient_fun = cell_col(taxa.trait_smpl_mat),
            Rowv = smpl_hc, Colv = trait_hc,
            # showticklabels = c(F, F), # remove x and y labels, essential with large matrices
            # row_side_colors = data.frame(Samples=row.names(taxa.trait_smpl_mat)),
            # row_side_palette = smpl_col,
            col_side_colors = data.frame(Taxa=gsub("_.*", "", colnames(taxa.trait_smpl_mat)),
                                         Trait_clusters=trait_analysis$clusters),
            col_side_palette = c(taxa_col_df, trait_analysis$clusters_col),
            # custom_hovertext = over_text,
            label_names = c("x", "y", "TPM"),
            file = taxa.trait_smpl_hm
  )
  write.table(t(cbind(Samples = c("Trait_clusters", row.names(taxa.trait_smpl_mat)), 
                      rbind(trait_analysis$clusters, taxa.trait_smpl_mat))), 
              file = file.path(out_dir, "taxa.trait_smpl_hm.tsv"), 
              col.names = F, row.names = T, quote = F, sep = "\t")
  
  
  #! other output
  write.table(MAG_ann_tpm_mat, file = file.path(out_dir, "MAG_ann_tpm_mat.tsv"), 
              col.names = T, row.names = F, quote = F, sep = "\t")
}



merge_FunLuca(merged_tpm = snakemake@input[["merged_tpm"]],
              MASTER_table = snakemake@input[["MASTER_table"]],
              out_dir = snakemake@output[["out_dir"]],
              ann_modules = snakemake@params[["ann_modules"]],
              taxa_col = snakemake@params[["taxa_col"]],
              # smpl_col = snakemake@params[["smpl_col"]],
              min_trait_occur = snakemake@params[["min_trait_occur"]],
              min_tpm = snakemake@params[["min_tpm"]],
              dist_mt = snakemake@params[["dist_mt"]],
              aggl_mt = snakemake@params[["aggl_mt"]]
)