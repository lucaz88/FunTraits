print_quality_profiles <- function(dir_parsed, is_paired, file_ext,
                                   dir_ASV, QC_token, ASV_arc) {
  cat("\n### PLOTTING QC\n")

    

  #! load packages
  suppressMessages(suppressWarnings(library(dada2)))
  suppressMessages(suppressWarnings(library(ggplot2)))
  suppressMessages(suppressWarnings(library(gridExtra)))
  
  
  
  #! set variables
  #!!! assumes "clip" is the last parsing module to be executed!!!
  fastqs <- list.files(path = file.path(dir_parsed, "clipped"), 
                       pattern = paste0(file_ext, "$"),
                       full.names = T, recursive = T)
  fastqs <- fastqs[!grepl("unknown", fastqs)]
  fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
  
  dir.create(dir_ASV, recursive = T, showWarnings = F)
  
  
  
  
  #! ASV inference
  if (is_paired) {     #!!! paired-end
    
    ##! prepare input data
    fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
    fnRs <- gsub("R1", "R2", fnFs) # Just the reverse read files
    sample.names <- gsub(file_ext, "", basename(fnFs))
    
    
    ##! Filter & Trimming
    ###! plot QC for inspection
    rdm_smpl <- sample(1:length(sample.names), 4)
    QC1 <- plotQualityProfile(fnFs[rdm_smpl])
    QC2 <- plotQualityProfile(fnRs[rdm_smpl])
    ggsave(file.path(dir_ASV, "parsed_seq_QUAL.pdf"), 
           arrangeGrob(QC1, QC2), width = 8, height = 12)
    
  } else if (!is_paired) { #!!! single-end
    
    ##! prepare input data
    fnFs <- fastqs
    sample.names <- gsub(file_ext, "", basename(fnFs))
    
    
    ##! Filter & Trimming
    ###! plot QC for inspection
    rdm_smpl <- sample(1:length(sample.names), 4)
    QC1 <- plotQualityProfile(fnFs[rdm_smpl])
    ggsave(file.path(dir_ASV, "parsed_seq_QUAL.pdf"), 
           QC1, width = 8, height = 12)
    
  } else {
    
    stop('\n\nInvalid argument provided for option "paired". Accepted values are "TRUE" or "FALSE"\n\n')
    
  }
  
  
  ##! create token and print instruction message
  system2(command = "touch", args = QC_token)
  system2(command = "touch", args = ASV_arc)
  
  return(cat(paste0("\nQuality report created:\n\n", 
                    "(1) Check ", file.path(dir_ASV, "parsed_seq_QUAL.pdf\n"),
                    "(2) Adjust the trimming parameter in config file\n",
                    "(3) Delete the token '", ASV_arc, "'\n",
                    "(4) Re-run the workflow\n\n\n")))
}


print_quality_profiles(dir_parsed = snakemake@input[["dir_parsed"]],
                       is_paired = snakemake@params[["is_paired"]],
                       file_ext = snakemake@params[["file_ext"]],
                       dir_ASV = snakemake@output[["dir_ASV"]],
                       QC_token = snakemake@output[["QC_token"]],
                       ASV_arc = snakemake@output[["ASV_arc"]]
)