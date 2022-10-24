parse_DMSP <- function(DMSP_files, 
                      blast_DMSP_tab) {
  #! load packages
  suppressMessages(suppressWarnings(library(stringr)))
  suppressMessages(suppressWarnings(library(tidyverse)))
  
  #! import data
  DMSP_ann1 <- lapply(DMSP_files, function(i) 
    tryCatch(read.delim(i, h=F, quote="", fill=F, check.names = F), 
             error = function(x) cat("No DMSP related genes in", i,"\n")))
  names(DMSP_ann1) <- gsub(".tsv", "", basename(DMSP_files))
  DMSP_ann1 <- do.call(rbind, DMSP_ann1)
  colnames(DMSP_ann1) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                          "qstart","qend","sstart","send","evalue","bitscore","annotation") # standard BLAST output format 6
  colnames(DMSP_ann1)[colnames(DMSP_ann1) == "qseqid"] <- "locus_tag"
  DMSP_ann1 <- cbind(filename=gsub("\\.[0-9]+$", "", row.names(DMSP_ann1)), DMSP_ann1)
  
  #! filtering
  DMSP_ann2 <- DMSP_ann1
  DMSP_ann2 <- DMSP_ann2[DMSP_ann2$evalue < 1e-70, ] # filter based on min e-value - Gärdes et al., 2013
  # DMSP_ann2 <- DMSP_ann2[DMSP_ann2$V12 > 60, ] # filter based on min bit score
  
  #! get coherent gene names
  DMSP_ann2$gene_ann <- sapply(DMSP_ann2$annotation, function(i) {
    ifelse(grepl("ddd", i, ignore.case = T), paste0("ddd", toupper(substr(gsub(".* ddd|.*=ddd", "", i, ignore.case = T), 1, 1))),
           ifelse(grepl("dmd", i, ignore.case = T), paste0("dmd", toupper(substr(gsub(".* dmd|.*=dmd", "", i, ignore.case = T), 1, 1))),
                  ifelse(grepl("acu", i, ignore.case = T), paste0("acu", toupper(substr(gsub(".* acu|.*=acu", "", i, ignore.case = T), 1, 1))), "NA"
                  )))
  })
  
  #! parse blast output
  DMSP_ann3 <- do.call(rbind, lapply(split(DMSP_ann2, f = DMSP_ann2$locus_tag), function(i) {
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
  if (is.null(DMSP_ann3)) {
    DMSP_ann3 <- DMSP_ann2
  }
  
  #! check for completeness
  DMSP_ann4 <- DMSP_ann3 %>% 
    mutate(value = 1) %>% 
    pivot_wider(id_cols = filename, names_from = gene_ann, names_sort = T,
                values_from = value, values_fn = sum, values_fill = 0)
  min_genes_1 <- 1 # out of 1 as ddd|ALMA genes are all redundant 
  min_genes_2 <- 3 # out of 4 for the dmd patway 
  # See 10.1146/annurev-marine-120710-100827 and 10.3389/fmicb.2017.00637
  DMSP_ann5 <- data.frame(filename = DMSP_ann4$filename,
                         DMSP_cleavage = apply(DMSP_ann4[, grepl("ddd|ALMA", colnames(DMSP_ann4))], 1, 
                                             function(x) ifelse(sum(x>0) >= min_genes_1, 1, 0)),
                         DMSP_demethylation = apply(DMSP_ann4[, grepl("dmd", colnames(DMSP_ann4))], 1, 
                                                  function(x) ifelse(sum(x>0) >= min_genes_2, 1, 0)))
  
  #! select genes
  DMSP_ann6 <- DMSP_ann3[(DMSP_ann3$filename %in% DMSP_ann5$filename[DMSP_ann5$DMSP_cleavage > 0] &
                           grepl("ddd|ALMA", DMSP_ann3$gene_ann)) |
                          (DMSP_ann3$filename %in% DMSP_ann5$filename[DMSP_ann5$DMSP_demethylation > 0] 
                           & grepl("dmd", DMSP_ann3$gene_ann)), ]
  DMSP_ann6$trait <- ifelse(grepl("ddd|ALMA", DMSP_ann6$gene_ann), "DMSP_cleavage",
                           ifelse(grepl("dmd", DMSP_ann6$gene_ann), "DMSP_demethylation", NA))
  
  #! save output
  write.table(DMSP_ann6, file = blast_DMSP_tab, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}


parse_DMSP(DMSP_files = snakemake@input[["DMSP_files"]],
           blast_DMSP_tab = snakemake@output[["blast_DMSP_tab"]]
)