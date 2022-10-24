make_master_table <- function(ann_tables, 
                             MASTER_table,
                             KM_long_names) {
  
  #! inport all annotaion tables
  all_ann <- lapply(ann_tables, function(i) read.delim(i, h=T))
  names(all_ann) <- gsub("_tab.tsv", "", basename(ann_tables))
  
  #! use Prokka GFF table as scaffold
  if (!any(grepl("prokka", names(all_ann)))) {
    stop("Somenthing must have gone wrong...there is no Prokka table in your results!\n")
  }
  ann_master_TAB <- all_ann[["prokka"]]
  all_ann <- all_ann[names(all_ann) != "prokka"]
  
  #! parse 'KEGG_KO' differently because the related traits are in a different file ('KEGG_KM')
  if (any(grepl("KEGG_KO", names(all_ann)))) {
    i_ann <- all_ann[[which(names(all_ann) == "KEGG_KO")]]
    i_ann <- i_ann[, c("filename","locus_tag","gene_ann")]
    colnames(i_ann)[colnames(i_ann) == "gene_ann"] <- "KO"
    
    ann_master_TAB <- merge(ann_master_TAB, i_ann, by = c("filename", "locus_tag"), all.x=T)
    all_ann <- all_ann[names(all_ann) != "KEGG_KO"]
  }
  
  #! parse 'KEGG_KM' differently because it needs KO annotation for mapping (as there are no gene names)
  if (any(grepl("KEGG_KM", names(all_ann)))) {
    i_ann <- all_ann[[which(names(all_ann) == "KEGG_KM")]]
    if (KM_long_names) {
      i_ann <- i_ann[, c("filename","gene_ann","KM_long")]
      colnames(i_ann)[colnames(i_ann) == "KM_long"] <- "trait"
    } else {
      i_ann <- i_ann[, c("filename","gene_ann","trait")]
    }
    colnames(i_ann)[colnames(i_ann) == "gene_ann"] <- "KO"
    colnames(i_ann)[colnames(i_ann) == "trait"] <- "KEGG_KM"
    
    ann_master_TAB <- merge(ann_master_TAB, i_ann, by = c("filename", "KO"), all.x=T)
    ann_master_TAB <- ann_master_TAB[, c(colnames(ann_master_TAB)[!colnames(ann_master_TAB) %in% c("KO","KEGG_KM")], "KO","KEGG_KM")]
    all_ann <- all_ann[names(all_ann) != "KEGG_KM"]
  }
  
  #! add all other annotated traits
  for (i in seq_along(all_ann)) {
    i_name <- names(all_ann)[i]
    i_ann <- all_ann[[i]]
    i_ann <- i_ann[, c("filename","locus_tag","trait")]
    colnames(i_ann)[colnames(i_ann) == "trait"] <- i_name
      
    ann_master_TAB <- merge(ann_master_TAB, i_ann, by = c("filename", "locus_tag"), all.x=T)

  }
  
  #! save output
  write.table(ann_master_TAB, file = MASTER_table, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}

make_master_table(ann_tables = snakemake@input[["ann_tables"]],
                  MASTER_table = snakemake@output[["MASTER_table"]],
                  KM_long_names = snakemake@params[["KM_long_names"]]
)