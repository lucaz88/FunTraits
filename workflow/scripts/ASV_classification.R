ASV_classification <- function(ASV_fasta, ASV_tab, 
                               dir_parsed, dir_ASV, class_method, ncore, 
                               classif_out, ASV_tab_wtax, ASV_tab_wtax_merged) {
  cat("\n### CLASSIFYING ASVs\n")
  
  
  
  #! load packages
  suppressMessages(suppressWarnings(library(seqinr)))
  
  
  
  #! set variables
  ASV_tab <- read.delim(ASV_tab, h = T)
  
  
  
  #! ASV classification
  implemented_methods <- c("sina_SILVA_138.1_SSU", "sina_SILVA_138.1_LSU", 
                           "blast_SILVA_138.1_SSU", "blast_SILVA_138.1_LSU", "blast_PR2_4.13_SSU")
  
  if (!class_method %in% implemented_methods) {
    stop(paste0("\n\nProvided 'class_method' is not supported. Valid options are:\n",
                paste0("- ", implemented_methods, "\n"), "\n\n"))
  }
  
  db_path <- "/gpfs/data/fs71579/lucaz/_DBs"
  ref_db <- ifelse(class_method == "sina_SILVA_138.1_SSU", file.path(db_path, "silva/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb"), 
                   ifelse(class_method == "sina_SILVA_138.1_LSU", file.path(db_path, "silva/SILVA_138.1_LSURef_NR99_30_06_20_opt.arb"), 
                          ifelse(class_method == "blast_SILVA_138.1_SSU", file.path(db_path, "silva/SILVA_138.1_SSURef_NR99_tax_silva.fasta"), 
                                 ifelse(class_method == "blast_SILVA_138.1_LSU", file.path(db_path, "silva/SILVA_138.1_LSURef_NR99_tax_silva.fasta"), 
                                        ifelse(class_method == "blast_PR2_4.13_SSU", file.path(db_path, "silva/SILVA_138.1_LSURef_NR99_30_06_20_opt.arb"),
                                               stop("\n\nSomething is wrong with the ref DB selection\n\n"))))))
  
  if (grepl("sina_", class_method)) {
    system2(command = "sina", 
            args = c(paste("-i", ASV_fasta),
                     paste("-o", file.path(dir_ASV, paste0(class_method, ".align.gz"))),
                     paste("-o", classif_out),
                     paste("--db", ref_db),
                     "--search", "--lca-fields tax_embl,tax_slv,tax_ltp", 
                     "--lca-quorum 0.5", "--search-max-result 10", 
                     paste("--threads", ncore), "--turn"),
            stdout = F, stderr = T)
    ASV_tax <- read.delim(classif_out, sep = ",", h=T)
    ASV_tax2 <- ASV_tax$lca_tax_slv[match(ASV_tab$ASVs, ASV_tax$name)]
    ASV_tax2[is.na(ASV_tax2)] <- "uncl."
    
  } else if (grepl("blast_", class_method)) {
    if (!file.exists(paste0(ref_db, ".ndb"))) {
      system2(command = "makeblastdb",
              args = c(paste("-in", ref_db),
                       "-dbtype nucl"))
    }
    system2(command = "blastn", 
            args = c(paste("-query", ASV_fasta),
                     paste("-db", ref_db),
                     paste("-out", classif_out),
                     "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'",
                     "-evalue 1e-5", paste("-num_threads", ncore), "-max_target_seqs 10"
            ))
    ASV_tax <- read.delim(classif_out, sep = "\t", h=F)
    colnames(ASV_tax) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                          "qstart","qend","sstart","send","evalue","bitscore","annotation")
    ASV_tax$annotation <- sapply(strsplit(ASV_tax$annotation, split = " "), 
                                 function(i) paste0(i[-1], collapse = "_"))
    # ASV_tax$Taxonomy[ASV_tax$evalue < 1e-80, ] <- "uncl."
    lca_quorum <- 0.5 # fraction of hits with same taxonomy to accept LCA
    ASV_tax <- do.call(rbind, lapply(split(ASV_tax, f = ASV_tax$qseqid), function(i) {
      #! find best hit
      i2 <- i[i[, "evalue"] == min(i[ ,"evalue"]) & 
               i[, "bitscore"] == max(i[, "bitscore"]), , drop=F]
      #! solve promiscuity
      if (nrow(i2) > 1) {
        #! check for quorum
        i_lca <- do.call(rbind, strsplit(i2[, "annotation"], ";"))
        good_ranks <- which(apply(i_lca, 2, function(j) max(table(j))/nrow(i_lca) ) >= lca_quorum)
        good_hits <- i_lca[, max(good_ranks)] == names(table(i_lca[, max(good_ranks)]))[which.max(table(i_lca[, max(good_ranks)]))]
        i_lca <- paste(i_lca[good_hits, good_ranks, drop=F][1,], collapse = ";")
        i2 <- i2[1, , drop=F]
        i2[1, "annotation"] <- i_lca
      }
      return(i2)
    }))
    ASV_tax2 <- ASV_tax$annotation[match(ASV_tab$ASVs, ASV_tax$qseqid)]
    ASV_tax2[is.na(ASV_tax2)] <- "uncl."
    
  } else {
    (stop("\n\nSomething is wrong with the method selected\n\n"))
    
  }
  
  
  
  
  
  #! save output
  ##! save ASVs table with Taxonomy
  ASV_tab.wtax <- data.frame(ASV_tab, Taxonomy=ASV_tax2)
  write.table(ASV_tab.wtax, file = ASV_tab_wtax,
              col.names = T, row.names = F, quote = F, sep = "\t")
  
  
  ##! save ASVs table MERGED BY taxonomy
  ASV_tab.wtax.merged <- ASV_tab.wtax[, -c(1,2,ncol(ASV_tab.wtax))]
  ASV_tab.wtax.merged <- rowsum(ASV_tab.wtax.merged, group = ASV_tax2)
  row.names(ASV_tab.wtax.merged)[is.na(row.names(ASV_tab.wtax.merged))] <- "Unclassified"
  ASV_tab.wtax.merged <- data.frame(Taxonomy=row.names(ASV_tab.wtax.merged), ASV_tab.wtax.merged)
  write.table(ASV_tab.wtax.merged, file = ASV_tab_wtax_merged,
              col.names = T, row.names = F, quote = F, sep = "\t")
  
}


ASV_classification(ASV_fasta = snakemake@input[["ASV_fasta"]],
                   ASV_tab = snakemake@input[["ASV_tab"]],
                   dir_parsed = snakemake@params[["dir_parsed"]],
                   dir_ASV = snakemake@params[["dir_ASV"]],
                   class_method = snakemake@params[["class_method"]],
                   ncore = snakemake@params[["ncore"]],
                   classif_out = snakemake@output[["classif_out"]],
                   ASV_tab_wtax = snakemake@output[["ASV_tab_wtax"]],
                   ASV_tab_wtax_merged = snakemake@output[["ASV_tab_wtax_merged"]]
)