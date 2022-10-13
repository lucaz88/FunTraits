parse_phytohormones = function(phyhorm_files, suppl_tab, 
                               blast_phytohormones_tab) {
  #! load packages
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressMessages(suppressWarnings(library(readxl)))
  
  #! import data
  phyhorm_ann1 = lapply(phyhorm_files, function(i) 
    tryCatch(read.delim(i, h=F, quote="", fill=F, check.names = F), 
             error = function(x) cat("No phytohormone related genes in", i,"\n")))
  names(phyhorm_ann1) = gsub(".tsv", "", basename(phyhorm_files))
  phyhorm_ann1 = do.call(rbind, phyhorm_ann1)
  colnames(phyhorm_ann1) = c("qseqid","sseqid","pident","length","mismatch","gapopen",
                          "qstart","qend","sstart","send","evalue","bitscore","annotation") # standard BLAST output format 6
  colnames(phyhorm_ann1)[colnames(phyhorm_ann1) == "qseqid"] = "locus_tag"
  phyhorm_ann1 = cbind(filename=gsub("\\.[0-9]+$", "", row.names(phyhorm_ann1)), phyhorm_ann1)
  phyhorm_ann1$gene_ann = gsub("\\|.*", "", phyhorm_ann1$annotation)
  phyhorm_ann1$gene_ann = gsub("K01502,", "", phyhorm_ann1$gene_ann) # K01502 has same gene sequence as K01501, but K01502 is not in KEGG map00380
  
  #! filtering
  phyhorm_ann2 = phyhorm_ann1
  phyhorm_ann2 = phyhorm_ann2[phyhorm_ann2$evalue < 5, ] # Zhang et al., 2019
  phyhorm_ann2 = phyhorm_ann2[phyhorm_ann2$bitscore > 60, ] # Zhang et al., 2019
  
  #! parse blast output
  phyhorm_ann3 = do.call(rbind, lapply(split(phyhorm_ann2, f = phyhorm_ann2$locus_tag), function(i) {
    #! find best hit
    i2 = i[i$evalue == min(i$evalue) & i$bitscore == max(i$bitscore), , drop=F]
    #! solve promiscuity (i.e. label as unclear if ties have different gene names)
    if (nrow(i2) > 1) {
      i2 = i2[1, , drop=F]
      if (length(unique(i$gene_ann)) == 1) {
        i2[1, "gene_ann"] = unique(i$gene_ann)
      } else { i2[1, "gene_ann"] = "unclear_ann" }
    }
    return(i2)
  }))
  
  #! check for any path in genomes using the annotated KOs
  phyhorm_list = as.data.frame(read_xlsx(suppl_tab, col_names = T, sheet = 1, skip = 3)) 
  phyhorm_list = phyhorm_list[!is.na(phyhorm_list$`Interaction traits`), ] # remove extra lines in Excel table
  phyhorm_list = split(phyhorm_list$`KEGG Orthology`, f = phyhorm_list$`Interaction traits`)
  phyhorm_list = lapply(phyhorm_list, function(i) strsplit(gsub("\\n| ", "", i), "\\|"))
  phyhorm_list = lapply(seq_along(phyhorm_list), function(i) data.frame(phyhorm=names(phyhorm_list)[i], KO=do.call(cbind, phyhorm_list[[i]])))
  phyhorm_list = do.call(rbind, phyhorm_list)
  phyhorm_list$len = sapply(phyhorm_list$KO, function(i) length(strsplit(i, ",")[[1]]))
  
  phyhorm_ann4 = lapply(phyhorm_list$KO, function(i) {
    i2 = strsplit(i, ",")[[1]] 
    with(phyhorm_ann2[phyhorm_ann2$gene_ann %in% i2, ], table(filename, gene_ann))
  })
  phyhorm_ann4 = lapply(seq_along(phyhorm_ann4), function(i) {
    i2 = apply(phyhorm_ann4[[i]], 1, function(j) sum(j>0))
    if (length(i2) > 0) {
      data.frame(phyhorm_list[i, ], filename=names(i2), ann_len=as.integer(i2))
    } # else {   data.frame(phyhorm_list[i,], filename=NA, ann_len=NA) }
  })
  phyhorm_ann4 = do.call(rbind, phyhorm_ann4)
  
  
  ##-- check for completeness (use same rules as for KMs)
  len_breaks = c(3) # set to NULL if no gap is allowed (i.e. allowed_gaps = 0)
  allowed_gaps = c(0,1)
  
  #! estimate gaps allowed for each path
  if (!is.null(len_breaks)) {
    len_bin <- .bincode(phyhorm_ann4$ann_len, breaks = c(0, len_breaks, Inf), right = F, include.lowest = T)
  } else {
    len_bin <- .bincode(phyhorm_ann4$ann_len, breaks = c(0, Inf), right = F, include.lowest = T)
  }
  allowed_gaps <- allowed_gaps[len_bin]
  
  #! add 'allowed_gaps' as bonus
  phyhorm_ann5 = phyhorm_ann4
  phyhorm_ann5$ann_len = phyhorm_ann4$ann_len + allowed_gaps
  
  #! filter for completeness
  phyhorm_ann6 = phyhorm_ann5
  phyhorm_ann6$compl = phyhorm_ann6$ann_len/phyhorm_ann6$len
  phyhorm_ann6$compl[phyhorm_ann6$compl >=1] = 1
  phyhorm_ann6$compl[phyhorm_ann6$compl <1] = 0
  phyhorm_ann6 = phyhorm_ann6[phyhorm_ann6$compl > 0, ]
  ##--
  
  
  #! select genes
  if (nrow(phyhorm_ann6) > 0) {
    ann_mask = phyhorm_ann6 %>% 
      group_by(filename, phyhorm) %>%
      summarise(gene_ann = unique(unlist(strsplit(KO, ",")))) %>%
      group_by(filename, gene_ann) %>%
      summarise(trait = paste0(unique(phyhorm), collapse = ";"))
    phyhorm_ann7 = merge(phyhorm_ann3, ann_mask, by = c("filename", "gene_ann"), all.x=F)
    phyhorm_ann7 = phyhorm_ann7[, c(colnames(phyhorm_ann3), "trait")]
  } else {
    phyhorm_ann7 = cbind(phyhorm_ann3[0, ], trait=c())
  }
  
  #! save output
  write.table(phyhorm_ann7, file = blast_phytohormones_tab, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}


parse_phytohormones(phyhorm_files = snakemake@input[["phyhorm_files"]],
                    suppl_tab = snakemake@input[["suppl_tab"]],
                    blast_phytohormones_tab = snakemake@output[["blast_phytohormones_tab"]]
)