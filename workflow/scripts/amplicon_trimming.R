amplicon_trimming <- function(parsed_dir, is_paired, file_ext,
                          minLen, truncLen, maxEE, insert_len, 
                          ncore, ASV_dir, ASV_token) {
  
  #! load packages
  library(dada2)
  library(ggplot2); library(gridExtra)

  
  #! set variables
  dir.create(ASV_dir, recursive = T)
  fastqs <- list.files(path = parsed_dir, pattern = paste0(file_ext, "$"),
                      full.names = T, recursive = T)
  fastqs <- fastqs[!grepl("unknown", fastqs)]
  fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
  
  
  #! ASV inference
  if (is_paired) {     #!!! paired-end
    
    ### prepare input data
    fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
    fnRs <- gsub("R1", "R2", fnFs) # Just the reverse read files
    sample.names <- gsub(file_ext, "", basename(fnFs))
    
    
    #### Filter & Trimming
    
    ## visually check same samples to decide whether or not seq have to be trimmed for low quality
    rdm_smpl <- sample(1:length(sample.names), 4)
    #! trick to plot the quality profiles and exit 
    if (!file.exists(file.path(ASV_dir, "parsed_seq_QUAL.pdf"))) {
      QC1 <- plotQualityProfile(fnFs[rdm_smpl])
      QC2 <- plotQualityProfile(fnRs[rdm_smpl])
      ggsave(file.path(ASV_dir, "parsed_seq_QUAL.pdf"), 
             grid.arrange(QC1, QC2, newpage = T), width = 8, height = 12)
      system2(command = "touch", args = ASV_token)
      return(cat(paste0("\nQuality report created:\n", 
                  file.path(ASV_dir, "parsed_seq_QUAL.pdf\n\n"),
                  "Check it, adjust the trimming parameter in config file, ",
                  "remove the token (", ASV_token,") and re-run the workflow\n\n\n")))
    }
    
    ## create output folder
    filtFs <- file.path(dir_filt, paste0(sample.names, "_R1.fastq"))
    filtRs <- file.path(dir_filt, paste0(sample.names, "_R2.fastq"))
    
    ## Filtration
    filt_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                             minLen = minLen, truncLen = truncLen,
                             maxN = 0, maxEE = maxEE, truncQ = 2, rm.phix = T,
                             compress = T, multithread = ncore, verbose = T)
    QC1_filt <- plotQualityProfile(filtFs[rdm_smpl])
    QC2_filt <- plotQualityProfile(filtRs[rdm_smpl])
    ggsave(file.path(ASV_dir, "parsed&trimmed_seq_QUAL.pdf"), 
           grid.arrange(QC1_filt, QC2_filt, newpage = TRUE), width = 8, height = 12)
    #! remove empty samples (badly sequenced or blanks)
    filtFs <- filtFs[filt_out[, 2] > 0]
    filtRs <- filtRs[filt_out[, 2] > 0]
    
    
  } else if (!is_paired) { #!!! single-end
    
    ### prepare input data
    fnFs <- fastqs
    sample.names <- gsub(file_ext, "", basename(fnFs))
    
    
    #### Filter & Trimming
    ## visually check same samples to decide whether or not seq have to be trimmed for low quality
    rdm_smpl <- sample(1:length(sample.names), 4)
    #! trick to plot the quality profiles and exit 
    if (!file.exists(file.path(ASV_dir, "parsed_seq_QUAL.pdf"))) {
      QC1 <- plotQualityProfile(fnFs[rdm_smpl])
      ggsave(file.path(ASV_dir, "parsed_seq_QUAL.pdf"), 
             QC1, width = 8, height = 12)
      system2(command = "touch", args = ASV_token)
      return(cat(paste0("\nQuality report created:\n", 
                        file.path(ASV_dir, "parsed_seq_QUAL.pdf\n\n"),
                        "Check it, adjust the trimming parameter in config file, ",
                        "remove the token (", ASV_token,") and re-run the workflow\n\n\n")))
    }
    
    ## create output folder
    filtFs <- file.path(dir_filt, paste0(sample.names, ".fastq.gz"))
    
    ## Filtration
    filt_out <- filterAndTrim(fnFs, filtFs,
                             minLen = minLen, truncLen = truncLen,
                             maxN = 0, maxEE = maxEE, truncQ = 2, rm.phix = T,
                             compress = T, multithread = ncore, verbose = T)
    QC1_filt <- plotQualityProfile(filtFs[rdm_smpl])
    ggsave(file.path(ASV_dir, "parsed&trimmed_seq_QUAL.pdf"), 
           QC1_filt, width = 8, height = 6)
    #! remove empty samples (badly sequenced or blanks)
    filtFs <- filtFs[filt_out[, 2] > 0]
    

  } else {
    
    stop('Invalid argument provided for option "paired". Accepted values are "TRUE" or "FALSE"')
    
    
  }
}


amplicon_trimming(parsed_dir = snakemake@input[["parsed_dir"]],
              is_paired = snakemake@params[["is_paired"]],
              file_ext = snakemake@params[["file_ext"]],
              minLen = snakemake@params[["minLen"]],
              truncLen = snakemake@params[["truncLen"]],
              maxEE = snakemake@params[["maxEE"]],
              insert_len = snakemake@params[["insert_len"]],
              ncore = snakemake@params[["ncore"]],
              ASV_dir = snakemake@output[["ASV_dir"]],
              ASV_token = snakemake@output[["ASV_token"]]
)