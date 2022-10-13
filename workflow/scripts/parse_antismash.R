parse_antismash = function(html_folder, gff_folder, 
                           antiSMASH_tab) {
  #! load packages
  suppressMessages(suppressWarnings(library(rvest)))
  suppressMessages(suppressWarnings(library(rtracklayer)))
  
  #! set variables
  sm_list = list.files(html_folder, pattern = "index.html$", full.names = T, recursive = T)
  gff_list = list.files(gff_folder, pattern = ".gff$", full.names = T, recursive = T)
  gff_list = gff_list[match(sapply(strsplit(sm_list, "/"), function(i) rev(i)[2]), 
                            sapply(strsplit(gff_list, "/"), function(i) rev(i)[2]))]
  
  #! import data - code tested with AntiSMASH v6.0.1 and r-rvest v1.0.2
  sm_ann1 = lapply(sm_list, function(i) tryCatch(read_html(i), error = function(x) cat("No output for ", i,"!\n")) )
  sm_ann2 = lapply(sm_ann1, function(i) {
    ann_tabs = html_table(html_elements(i, "table"))
    # there is one table for each contig in genome fasta (uses prokka ids) 
    # plus final summary table and one table for each annotated sm pathway
    # we want to keep only the contig tables (!include "From" and "To" columns)
    if (length(ann_tabs) != 0) {
      which_contig = unlist(lapply(lapply(ann_tabs, colnames), function(j) 
        any(j %in% c("From","To"))))
      contig_ann_tabs = ann_tabs[which_contig]
      gff_locus = html_text(html_elements(i, "strong")) # retrieve the contig names
      gff_locus = gff_locus[which_contig]
      names(contig_ann_tabs) = gff_locus[1:length(contig_ann_tabs)]
      if (sum(which_contig) > 1) { # genomes with only one contig will not have a summary table
        contig_ann_tabs = contig_ann_tabs[-length(contig_ann_tabs)] # discard the summary tab
      }
      contig_ann_tabs = do.call(rbind, contig_ann_tabs)
      contig_ann_tabs = data.frame(contig_id=gsub("\\..*", "", row.names(contig_ann_tabs)), contig_ann_tabs)
      
      return(contig_ann_tabs)
    } 
  })
  names(sm_ann2) = sapply(strsplit(sm_list, "/"), function(i) rev(i)[2])
  sm_ann3 = sm_ann2[unlist(lapply(sm_ann2, length)) > 0]
  sm_ann4 = lapply(seq_along(sm_ann3), function(i) { # map annotations to CDS
    i_sm = sm_ann3[[i]]
    i_sm$From = as.integer(gsub(",", "", i_sm$From))
    i_sm$To = as.integer(gsub(",", "", i_sm$To))
    i_gff = as.data.frame(readGFF(gff_list[i]))
    i_gff = i_gff[i_gff$type == "CDS", ]
    i_gff$seqid = sapply(strsplit(as.character(i_gff$seqid), "\\|"), function(j) rev(j)[[1]])
    i_sm = lapply(seq_len(nrow(i_sm)), function(j) {
      j_sm = i_sm[j, , drop=T]
      i_gff_j = i_gff[i_gff$seqid == j_sm["contig_id"] &
                        (i_gff$start <= j_sm["From"] | i_gff$start < j_sm["To"]) &
                        (i_gff$end > j_sm["From"] | i_gff$end >= j_sm["To"]), c("seqid","locus_tag")]
      merge(i_gff_j, j_sm, by.x="seqid", by.y="contig_id", all=T)
    })
    i_sm = do.call(rbind, i_sm)
    return(i_sm)
  })
  names(sm_ann4) = names(sm_ann3)
  sm_ann4 = do.call(rbind, sm_ann4)
  sm_ann4 = cbind(filename=gsub("\\.[0-9]+$", "", row.names(sm_ann4)), sm_ann4)
  sm_ann5 = sm_ann4[!grepl("unknown|,", sm_ann4$Type), ]
  
  #! parse output
  sm_ann5 = sm_ann5[, c(colnames(sm_ann5)[colnames(sm_ann5) != "Type"], "Type")]
  colnames(sm_ann5)[colnames(sm_ann5) == "Type"] = "trait"
  
  #! save output
  write.table(sm_ann5, file = antiSMASH_tab, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}


parse_antismash(html_folder = snakemake@input[["html_folder"]],
                gff_folder = snakemake@input[["gff_folder"]],
                antiSMASH_tab = snakemake@output[["antiSMASH_tab"]]
)