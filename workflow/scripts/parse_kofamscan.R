parse_kofamscan = function(ko_files, 
                           KEGG_KO_tab) {
  #! load packages
  suppressMessages(suppressWarnings(library(webchem)))
  suppressMessages(suppressWarnings(library(tidyverse)))
  
  #! import data
  KO_ann1 = lapply(ko_files, function(i) read.delim(i, h=T, quote="", fill=F, check.names = F))
  KO_ann2 = lapply(KO_ann1, function(i) {
    i2 = strsplit(i[2:nrow(i), ], "[ ]+")
    i3 = lapply(i2, function(j) c(j[1:6], paste(j[7:length(j)], collapse = " ")))
    i4 = do.call(rbind, i3)
    colnames(i4) = c("sign","locus_tag","gene_ann","thrshld","score","E.value","KO.definition")
    i5 = as.data.frame(apply(i4, 2, unlist))
    return(i5)
  })
  names(KO_ann2) = gsub("_ko.txt", "", basename(ko_files))
  KO_ann2 = do.call(rbind, KO_ann2)
  KO_ann2 = data.frame(filename=gsub("\\.[0-9]+$", "", row.names(KO_ann2)), KO_ann2)
  KO_ann2 = KO_ann2[, c("filename","locus_tag","thrshld","score","E.value","sign","KO.definition","gene_ann")]
  
  #! filtering
  KO_ann3 = KO_ann2[grepl("\\*", KO_ann2$sign), ] # keep only signif
  
  #! resolve promiscuity
  KO_ann4 = KO_ann3 %>% group_by(locus_tag) %>%
    # keep best match that is not an 'uncharacterized protein' if there is any
    filter( if (any(KO.definition != "uncharacterized protein")) { 
      KO.definition != "uncharacterized protein"
    } else KO.definition == KO.definition ) %>%
    filter(score == max(score)) %>%
    # keep match with highest threshold
    filter( if (length(thrshld) > 1) { 
      thrshld == max(thrshld)
    } else thrshld == thrshld )
  #---
  
  #! save output
  write.table(KO_ann4, file = KEGG_KO_tab, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}


parse_kofamscan(ko_files = snakemake@input[["ko_files"]], 
                KEGG_KO_tab = snakemake@output[["KEGG_KO_tab"]]
)