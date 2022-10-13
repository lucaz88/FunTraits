# #___ for testing
# setwd("~/VSC5_data/genomics/B03_G00002_Cornus/")
# QC_report_dir <- "1_cleaned/99_QC_clean/Cornus_mas_v1.0"
# out_dir <- "1_cleaned/99_QC_clean/Cornus_mas_v1.0"
# QC_report_dir <- "1_cleaned/00_QC_raw/Cornus_mas_v1.0"
# out_dir <- "1_cleaned/00_QC_raw/Cornus_mas_v1.0"
# #___ for testing



#___Parse command line arguments - calling via Rscript
QC_report_dir <- commandArgs(TRUE)[1]



# plot_len_acc_fq <- function(QC_report_dir) {
  
  #! load packages & functions
  # suppressMessages(suppressWarnings(library(ShortRead)))
  
  
  #! import data
  QC_report <- list.files(QC_report_dir, pattern = ".zip", recursive = T, full.names = T)
  smpl_name <- gsub(".zip", "", basename(QC_report))
  qc_tab <- read.delim(unz(QC_report, paste0(smpl_name, "/fastqc_data.txt")), h=F)
  
  
  #! parse data
  idx_start <- which(grepl(">>", qc_tab$V1) & ! grepl(">>END_MODULE", qc_tab$V1))
  idx_end <- which(grepl(">>", qc_tab$V1) & grepl(">>END_MODULE", qc_tab$V1))
  qc_tabs <- lapply(seq_along(idx_start), function(i) {
    i_tab <- qc_tab[idx_start[i]:idx_end[i], ]
    colnames(i_tab) <- gsub("\\#", "", i_tab[2, ])
    i_tab <- i_tab[-c(1:2, nrow(i_tab)), ]
  })
  
  smpl_accuracy <- qc_tabs[[3]]
  smpl_accuracy$exp_err <- 10^-(as.integer(smpl_accuracy$Quality)/10)
  smpl_accuracy$accuracy <- 1-smpl_accuracy$exp_err
  smpl_len <- qc_tabs[[7]]
  smpl_len$Count <- as.numeric(smpl_len$Count) 
  smpl_len$Length_mean <- sapply(strsplit(smpl_len$Length, "-"), function(i) 
    round(mean(as.integer(i)), 0))
  
  
  #! plot
  pdf(file.path(QC_report_dir, "accuracy_length.pdf"), width = 8, height = 4)
  par(mfrow=c(1,2), 
      mar=c(7,4,3,0))
  tot_reads = paste(round(sum(as.integer(smpl_accuracy$Count))/1e6, 1), "Mb")
  hist(smpl_accuracy$accuracy, breaks = 10, freq = F, 
       xlab = "", ylab = paste0("#reads (", tot_reads, ")"), main = "Accuracy")
  barplot(smpl_len$Count, names.arg = smpl_len$Length_mean, 
          xlab = "", ylab = "", main = "Length - bp", 
          las = 2)
  dev.off()
# }



# plot_len_acc_fq(QC_report_dir = snakemake@input[["QC_report_dir"]]
# )