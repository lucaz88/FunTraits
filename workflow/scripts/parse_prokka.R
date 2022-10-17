parse_prokka = function(gff_folder, 
                        GFF_table) {
  #! load packages
  suppressMessages(suppressWarnings(library(rtracklayer)))
  
  #! set variables
  gff_list = sapply(gff_folder, function(i)
    list.files(i, pattern = ".gff$", full.names = T, recursive = T))
  
  #! parsing
  gff_tab = lapply(gff_list, readGFF)
  gff_tab = lapply(gff_tab, function(i) as.data.frame(i[i$type != "gene", c("locus_tag","type","start","end","strand","gene","product","db_xref","eC_number")]))
  # gff_tab = lapply(gff_tab, function(i) as.data.frame(i[i$type == "CDS", c("locus_tag","start","end","strand","gene","product","db_xref","eC_number")]))
  names(gff_tab) = sapply(strsplit(gff_list, "/"), function(i) rev(i)[2])
  gff_tab = do.call(rbind, gff_tab)
  # colnames(gff_tab)[colnames(gff_tab) == "ID"] = "geneID"
  gff_tab = data.frame(filename=gsub("\\.[0-9]+$", "", row.names(gff_tab)), gff_tab)

  write.table(gff_tab, file = GFF_table, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}


parse_prokka(gff_folder = snakemake@input[["gff_folder"]], 
             GFF_table = snakemake@output[["GFF_table"]]
)