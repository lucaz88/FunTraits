parse_phytohormones <- function(vibferr_files, 
                               blast_vibrioferrin_tab) {
  #! load packages
  suppressMessages(suppressWarnings(library(stringr)))
  suppressMessages(suppressWarnings(library(tidyverse)))
  
  #! import data
  vibferr_ann1 <- lapply(vibferr_files, function(i) 
    tryCatch(read.delim(i, h=F, quote="", fill=F, check.names = F), 
             error = function(x) cat("No vibrioferrin related genes in", i,"\n")))
  names(vibferr_ann1) <- gsub(".tsv", "", basename(vibferr_files))
  vibferr_ann1 <- do.call(rbind, vibferr_ann1)
  colnames(vibferr_ann1) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                             "qstart","qend","sstart","send","evalue","bitscore","annotation") # standard BLAST output format 6
  colnames(vibferr_ann1)[colnames(vibferr_ann1) == "qseqid"] <- "locus_tag"
  vibferr_ann1 <- cbind(filename=gsub("\\.[0-9]+$", "", row.names(vibferr_ann1)), vibferr_ann1)
  
  #! filtering
  vibferr_ann2 <- vibferr_ann1
  vibferr_ann2 <- vibferr_ann2[vibferr_ann2$evalue < 1e-5, ] # filter based on min e-value, Amin et al, 2012
  vibferr_ann2 <- vibferr_ann2[vibferr_ann2$pident > 30, ] # filter based on min similarity, Amin et al, 2012
  
  #! get coherent gene names
  vibferr_ann2$gene_ann <- sapply(vibferr_ann2$annotation, function(i) {
    ifelse(grepl("pvu", i, ignore.case = T), paste0("pvu", substr(gsub(".*pvu", "", i, ignore.case = T), 1, 1)),
           ifelse(grepl("pvs", i, ignore.case = T), paste0("pvs", substr(gsub(".*pvs", "", i, ignore.case = T), 1, 1)), NA
           ))
  }) 
  
  #! parse blast output
  vibferr_ann3 <- do.call(rbind, lapply(split(vibferr_ann2, f = vibferr_ann2$locus_tag), function(i) {
    #! find best hit
    i2 <- i[i$evalue == min(i$evalue) & i$bitscore == max(i$bitscore), , drop=F]
    #! solve promiscuity (i.e. label as unclear if ties have different gene names)
    if (nrow(i2) > 1) {
      i2 <- i2[1, , drop=F]
      if (length(unique(i$gene_ann)) == 1) {
        i2[1, "gene_ann"] <- unique(i$gene_ann)
      } else { i2[1, "gene_ann"] <- "unclear_ann" }
    }
    return(i2)
  }))
  
  #! handle exception of no hits passing the filtering step 
  if (is.null(vibferr_ann3)) {
    vibferr_ann3 <- vibferr_ann2
  }
  
  #! check for completeness
  vibferr_ann4 <- vibferr_ann3 %>% 
    mutate(value = 1) %>% 
    pivot_wider(id_cols = filename, names_from = gene_ann, names_sort = T,
                values_from = value, values_fn = sum, values_fill = 0)
  min_genes <- 4 # out of 5 total genes in both operons
  vibferr_ann5 <- data.frame(filename=vibferr_ann4$filename,
                            pvsABCDE=apply(vibferr_ann4[, grepl("pvs", colnames(vibferr_ann4))], 1, function(x) ifelse(sum(x>0) >= min_genes, 1, 0)),
                            pvuABCDE=apply(vibferr_ann4[, grepl("pvu", colnames(vibferr_ann4))], 1, function(x) ifelse(sum(x>0) >= min_genes, 1, 0)))
  
  #! select genes
  vibferr_ann6 <- vibferr_ann3[(vibferr_ann3$filename %in% vibferr_ann5$filename[vibferr_ann5$pvsABCDE > 0] &
                                 grepl("pvs", vibferr_ann3$gene_ann)) |
                                (vibferr_ann3$filename %in% vibferr_ann5$filename[vibferr_ann5$pvuABCDE > 0] 
                                 & grepl("pvu", vibferr_ann3$gene_ann)), ]
  vibferr_ann6$trait <- ifelse(grepl("pvs", vibferr_ann6$gene_ann), "pvsABCDE",
                              ifelse(grepl("pvu", vibferr_ann6$gene_ann), "pvuABCDE", NA))
  
  #! save output
  write.table(vibferr_ann6, file = blast_vibrioferrin_tab, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}


parse_phytohormones(vibferr_files = snakemake@input[["vibferr_files"]],
                    blast_vibrioferrin_tab = snakemake@output[["blast_vibrioferrin_tab"]]
)