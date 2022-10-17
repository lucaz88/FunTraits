parse_dbCAN = function(CAZy_folder, 
                       CAZy_subs, 
                       CAZy_table) {
  #! load packages
  # library(rvest)
  
  #! import raw data
  CAZy_ann1 = lapply(CAZy_folder, function(i) 
    read.delim(file.path(i, "overview.txt"), h=T, sep = "\t"))
  names(CAZy_ann1) = sapply(strsplit(CAZy_folder, "/"), function(i) rev(i)[1])
  CAZy_ann1 = do.call(rbind, CAZy_ann1)
  CAZy_ann1 = cbind(filename=gsub("\\.[0-9]+$", "", row.names(CAZy_ann1)), CAZy_ann1)
  colnames(CAZy_ann1)[colnames(CAZy_ann1) == "Gene.ID"] = "locus_tag"
  
  #! filtering
  CAZy_ann2 = CAZy_ann1[CAZy_ann1$X.ofTools >= 2, ] # positive hit with 2 tools out of 3
  
  #! parsing data
  CAZy_ann2$gene_ann = apply(CAZy_ann2[, colnames(CAZy_ann2) %in% c("HMMER","eCAMI","DIAMOND")], 1, function(i) {
    clean_ann = unique(gsub("\\(.*", "", i[i != "-"]))
    paste0(clean_ann, collapse = ";")
  })
  # uniq_CAZy_fam = unique(CAZy_ann2$gene_ann)
  # uniq_CAZy_fam = strsplit(uniq_CAZy_fam, "\\+|;")
  # uniq_CAZy_fam = lapply(uniq_CAZy_fam, unique)
  # uniq_CAZy_fam_desc = lapply(uniq_CAZy_fam, function(i) {
  #   ii = sapply(i, function(j) {
  #     j2 = read_html(url(paste0('http://www.cazy.org/', j,'.html')), encoding = "UTF-8", verbose=T)
  #     j3 = as.data.frame(html_table(html_elements(ii, "table"))[[1]])
  #     j4 = j3[grepl("Activities in ", j3[, 1]), 2]
  #     return(j4)
  #   })
  # })
  cazy_substrates = read.delim(CAZy_subs, h=F)
  CAZy_ann2$trait = sapply(CAZy_ann2$gene_ann, function(i) {
    i2 = strsplit(i, "\\+|;")[[1]]
    i2_sub = paste0(unique(cazy_substrates$V1[cazy_substrates$V2 %in% i2]), collapse = "|")
    ifelse(i2_sub != "", i2_sub, NA)
  })
  
  #! save output
  write.table(CAZy_ann2, file = CAZy_table, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}


parse_dbCAN(CAZy_folder = snakemake@input[["CAZy_folder"]],
            CAZy_subs = snakemake@params[["CAZy_subs"]],
            CAZy_table = snakemake@output[["CAZy_table"]]
)