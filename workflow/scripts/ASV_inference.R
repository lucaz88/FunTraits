ASV_inference <- function(dir_parsed, is_paired, file_ext,
                          minLen, truncLen, maxEE, BAND_SIZE, insert_len, ncore,
                          dir_ASV, ASV_fasta, ASV_tab) {

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
  
  dir_trim = file.path(dir_parsed, "trimmed")
  dir.create(dir_trim, recursive = T, showWarnings = F)
  dir.create(dir_ASV, recursive = T, showWarnings = F)
  
  
  
  #! ASV inference
  if (is_paired) {     #!!! paired-end
    
    ##! prepare input data
    fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
    fnRs <- gsub("R1", "R2", fnFs) # Just the reverse read files
    sample.names <- gsub(file_ext, "", basename(fnFs))
    
    
    ##! Filter & Trimming
    ###! plot QC
    rdm_smpl <- sample(1:length(sample.names), 4)
    QC1 <- plotQualityProfile(fnFs[rdm_smpl])
    QC2 <- plotQualityProfile(fnRs[rdm_smpl])
    ggsave(file.path(dir_ASV, "parsed_seq_QUAL.pdf"), 
           arrangeGrob(QC1, QC2), width = 8, height = 12)
    ###! create output folder
    filtFs <- file.path(dir_trim, paste0(sample.names, "_R1.fastq.gz"))
    filtRs <- file.path(dir_trim, paste0(sample.names, "_R2.fastq.gz"))
    ###! Filtration
    cat("\n### TRIMMING READS\n")
    filt_out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                              minLen = eval(parse(text = truncLen)), 
                              truncLen = eval(parse(text = truncLen)),
                              maxEE = eval(parse(text = maxEE)),
                              maxN = 0, truncQ = 2, rm.phix = T,
                              compress = T, multithread = ncore)
    print(head(filt_out))
    QC1_filt <- plotQualityProfile(filtFs[rdm_smpl])
    QC2_filt <- plotQualityProfile(filtRs[rdm_smpl])
    ggsave(file.path(dir_ASV, "parsed&trimmed_seq_QUAL.pdf"), 
           arrangeGrob(QC1_filt, QC2_filt), width = 8, height = 12)
    filtFs <- filtFs[filt_out[, 2] > 0] # remove empty samples (badly sequenced or blanks)
    filtRs <- filtRs[filt_out[, 2] > 0] # remove empty samples (badly sequenced or blanks)
    sample.names <- sample.names[filt_out[, 2] > 0]
    
    
    ##! Learn error rates
    cat("\n### MODELLING ERROR RATES\n")
    errF <- learnErrors(filtFs, multithread=ncore, randomize = T, nbases = 2e8, verbose=T) # randomize as we have samples from multiple sequencing
    errR <- learnErrors(filtRs, multithread=ncore, randomize = T, nbases = 2e8, verbose=T) # , MAX_CONSIST = 15
    # sanity check: if the black line (estimated rates) fits the observed rates (points), 
    # and if the black line (estimated rates) drop with increased quality as expected under the nominal definition of the Q-value (red line)
    errF_plot <- plotErrors(errF, nominalQ=T)
    errR_plot <- plotErrors(errR, nominalQ=T)
    ggsave(file.path(dir_ASV, "parsed&trimmed_seq_ERR.pdf"), 
           arrangeGrob(errF_plot, errR_plot), width = 8, height = 10)
    
    
    ##! Sample Inference
    cat("\n### ASV INFERENCE\n")
    dadaFs <- dada(filtFs, err=errF, multithread=ncore, BAND_SIZE=BAND_SIZE, pool=T, verbose=T) 
    dadaRs <- dada(filtRs, err=errR, multithread=ncore, BAND_SIZE=BAND_SIZE, pool=T, verbose=T)
    names(dadaFs) <- sample.names
    names(dadaRs) <- sample.names
    
    
    ##! Merge paired-ends
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap=12, verbose=T)
    head(mergers[[1]]) # Inspect the merger data.frame from the first sample
    
    
    ##! SEQ table
    seqtab <- makeSequenceTable(mergers)
    # seqtab2 <- seqtab
    seqtab2 <- seqtab[, nchar(colnames(seqtab)) %in% eval(parse(text = insert_len))] # remove sequence that appear too long or short
    
    
    ##! Remove chimeras
    cat("\n### REMOVING CHIMERAS\n")
    seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=ncore, verbose=T)
    print(sum(seqtab.nochim)/sum(seqtab)) # if the % of good seq is too low try check for primers/adapters still unremoved
    
    
    ##! Export
    cat("\n### EXPORTING RESULTS\n")
    ###! ASVs table
    write.table(data.frame(ASVs=paste0("ASV_", 1:length(colnames(seqtab.nochim))), 
                           seq=colnames(seqtab.nochim), t(seqtab.nochim)), 
                ASV_tab, col.names = T, row.names = F, quote = F, sep = "\t")
    ###! ref sequences for ASVs
    write.table(cbind(paste0(">ASV_", 1:length(colnames(seqtab.nochim))), colnames(seqtab.nochim)),
                ASV_fasta, col.names=F, row.names=F, sep="\n", quote=F) 
    ###! R projects
    save(sample.names, filtFs, filtRs, filt_out, errF, errR, dadaFs, dadaRs,
         mergers, seqtab, seqtab.nochim, 
         file = file.path(dir_ASV, "dada2_PE.RData"))

    
  } else if (!is_paired) { #!!! single-end
    
    ##! prepare input data
    fnFs <- fastqs
    sample.names <- gsub(file_ext, "", basename(fnFs))
    
    
    ##! Filter & Trimming
    ###! plot base quality for inspection
    rdm_smpl <- sample(1:length(sample.names), 4)
    QC1 <- plotQualityProfile(fnFs[rdm_smpl])
    ggsave(file.path(dir_ASV, "parsed_seq_QUAL.pdf"), 
           QC1, width = 8, height = 12)
    ###! create output folder
    filtFs <- file.path(dir_trim, paste0(sample.names, ".fastq.gz"))
    ###! Filtration
    cat("\n### TRIMMING READS\n")
    filt_out <- filterAndTrim(fnFs, filtFs,
                              minLen = eval(parse(text = truncLen)), 
                              truncLen = eval(parse(text = truncLen)),
                              maxEE = eval(parse(text = maxEE)),
                              maxN = 0, truncQ = 2, rm.phix = T,
                              compress = T, multithread = ncore)
    print(head(filt_out))
    QC1_filt <- plotQualityProfile(filtFs[rdm_smpl])
    ggsave(file.path(dir_ASV, "parsed&trimmed_seq_QUAL.pdf"), 
           QC1_filt, width = 8, height = 6)
    filtFs <- filtFs[filt_out[, 2] > 0] #! remove empty samples (badly sequenced or blanks)
    sample.names <- sample.names[filt_out[, 2] > 0]
    
    
    ##! Learn error rates
    cat("\n### MODELLING ERROR RATES\n")
    errF <- learnErrors(filtFs, multithread=ncore, randomize = T, nbases = 2e8, verbose=T) # randomize as we have samples from multiple sequencing
    # sanity check: if the black line (estimated rates) fits the observed rates (points), 
    # and if the black line (estimated rates) drop with increased quality as expected under the nominal definition of the Q-value (red line)
    errF_plot <- plotErrors(errF, nominalQ=T)
    ggsave(file.path(dir_ASV, "parsed&trimmed_seq_ERR.pdf"), 
           errF_plot, width = 8, height = 5)
    
    
    ##! Sample Inference
    cat("\n### ASV INFERENCE\n")
    dadaFs <- dada(filtFs, err=errF, multithread=ncore, BAND_SIZE=BAND_SIZE, pool=T, verbose=T)
    names(dadaFs) <- sample.names
    
    
    ##! SEQ table
    seqtab <- makeSequenceTable(dadaFs)
    # seqtab2 <- seqtab
    seqtab2 <- seqtab[, nchar(colnames(seqtab)) %in% eval(parse(text = insert_len))] # remove sequence that appear anomally long or short
    
    
    ##! Remove chimeras
    cat("\n### REMOVING CHIMERAS\n")
    seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=ncore, verbose=T)
    print(sum(seqtab.nochim)/sum(seqtab)) # if the % of good seq is too low try check for primers/adapters still unremoved
    
    
    ##! Export
    cat("\n### EXPORTING RESULTS\n")
    ###! ASVs table
    write.table(data.frame(ASVs=paste0("ASV_", 1:length(colnames(seqtab.nochim))), 
                           seq=colnames(seqtab.nochim), t(seqtab.nochim)), 
                ASV_tab, col.names = T, row.names = F, quote = F, sep = "\t")
    ###! ref sequences for ASVs
    write.table(cbind(paste0(">ASV_", 1:length(colnames(seqtab.nochim))), colnames(seqtab.nochim)),
                ASV_fasta, col.names=F, row.names=F, sep="\n", quote=F)
    ###! R projects
    save(sample.names, filtFs, filt_out, errF, dadaFs,
         seqtab, seqtab.nochim, 
         file = file.path(dir_ASV, "dada2_SE.RData"))

  } else {
    
    stop('\n\nInvalid argument provided for option "paired". Accepted values are "TRUE" or "FALSE"\n\n')
    
    
    
  }
}


ASV_inference(dir_parsed = snakemake@input[["dir_parsed"]],
              is_paired = snakemake@params[["is_paired"]],
              file_ext = snakemake@params[["file_ext"]],
              minLen = snakemake@params[["minLen"]],
              truncLen = snakemake@params[["truncLen"]],
              maxEE = snakemake@params[["maxEE"]],
              BAND_SIZE = snakemake@params[["BAND_SIZE"]],
              insert_len = snakemake@params[["insert_len"]],
              ncore = snakemake@params[["ncore"]],
              dir_ASV = snakemake@output[["dir_ASV"]],
              ASV_fasta = snakemake@output[["ASV_fasta"]],
              ASV_tab = snakemake@output[["ASV_tab"]]
)