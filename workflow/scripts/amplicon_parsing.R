amplicon_reorient <- function(list_in, is_paired, file_ext, forA, revA,
                              ncore, list_out) {

  #! load packages
  suppressMessages(suppressWarnings(library(ShortRead)))
  suppressMessages(suppressWarnings(library(future.apply)))
  
  
  
  #! re-orient reads
  if (is_paired) {     # paired-end
    # Typically, Illumina sequencing generates files where forward and reversed reads
    # are mixed 50%-50% in both R1 and R2 files
    
    plan(multicore, workers = ncore)
    future_sapply(seq_along(list_in), function(i) {
      
      
      ##! set variables
      iR1 <- list_in[i]
      iR2 <- gsub("R1", "R2", iR1)
      iR1_out <- list_out[i]
      iR2_out <- gsub("R1", "R2", iR1_out)
      dir_tmp <- gsub(file_ext, "", iR1_out)
      dir.create(dir_tmp, recursive = T, showWarnings = F)
      
      
      ##! split forward reads and rev-complements
      system2(command = "cutadapt", args = c(paste0("-g one=", forA), paste0("-g two=", revA),
                                             paste0("-G one=", forA), paste0("-G two=", revA),
                                             "--overlap 10", "--action=retain", "-e 0.15",
                                             paste0("-o ", dir_tmp, "/{name}.R1.fastq.bz2"),
                                             paste0("-p ", dir_tmp, "/{name}.R2.fastq.bz2"), iR1, iR2),
              stdout=F
      )
      
      
      # ##! re-orient making the rev-comp
      # system2(command = "reformat.sh", args = c(paste0("in=", dir_tmp, "/two.R1.fastq.bz2"),
      #                                           paste0("out=", dir_tmp, "/two.R1.RC.fastq.bz2"),
      #                                           "rcomp"),
      #         stdout=F, stderr=F
      # )
      # system2(command = "reformat.sh", args = c(paste0("in=", dir_tmp, "/two.R2.fastq.bz2"),
      #                                           paste0("out=", dir_tmp, "/two.R2.RC.fastq.bz2"),
      #                                           "rcomp"),
      #         stdout=F, stderr=F
      # )
      
      
      ##! merge back
      # system2(command = "cat", args = c(paste0(dir_tmp, "/one.R1.fastq.bz2"),
      #                                   paste0(dir_tmp, "/two.R1.RC.fastq.bz2"),
      #                                   paste(">", iR1_out)))
      # system2(command = "cat", args = c(paste0(dir_tmp, "/one.R2.fastq.bz2"),
      #                                   paste0(dir_tmp, "/two.R2.RC.fastq.bz2"),
      #                                   paste(">", iR2_out)))
      system2(command = "cat", args = c(paste0(dir_tmp, "/one.R1.fastq.bz2"),
                                        paste0(dir_tmp, "/two.R2.fastq.bz2"),
                                        paste(">", iR1_out)))
      system2(command = "cat", args = c(paste0(dir_tmp, "/one.R2.fastq.bz2"),
                                        paste0(dir_tmp, "/two.R1.fastq.bz2"),
                                        paste(">", iR2_out)))
      
      
      ##! clean-up
      system2(command = "rm", args = c("-rf", dir_tmp))
    })
    
  } else if (!is_paired) { # single-end
    plan(multicore, workers = ncore)
    future_sapply(seq_along(list_in), function(i) {
      
      
      ##! set variables
      iR1 <- list_in[i]
      iR1_out <- list_out[i]
      dir_tmp <- gsub(file_ext, "", iR1_out)
      dir.create(dir_tmp, recursive = T, showWarnings = F)
      
      
      ##! split forward reads and rev-complements
      system2(command = "cutadapt", args = c(paste0("-g one=", forA), paste0("-g two=", revA),
                                             "--overlap 10", "--action=retain", "-e 0.15",
                                             paste0("-o ", dir_tmp, "/{name}.fastq.bz2"), iR1),
              stdout=F
      )
      
      
      # ##! re-orient making the rev-comp - PacBio
      # system(paste("bzip2 -d", file.path(dir_tmp, "two.fastq.bz2")))
      # writeFastq(reverseComplement(readFastq(file.path(dir_tmp, "two.fastq"))),
      #            file.path(dir_tmp, "two_RC.fastq"), compress = F)
      # system(paste("bzip2 -z", file.path(dir_tmp, "two_RC.fastq")))
      
      ##! re-orient making the rev-comp - all other technologies
      system2(command = "reformat.sh", args = c(paste0("in=", dir_tmp, "/two.fastq.bz2"),
                                                paste0("out=", dir_tmp, "/two_RC.fastq.bz2"),
                                                "rcomp"),
              stdout=F, stderr=F
      )
      
      
      ##! merge back
      system2(command = "cat", args = c(paste0(dir_tmp, "/one.fastq.bz2"),
                                        paste0(dir_tmp, "/two_RC.fastq.bz2"),
                                        paste(">", iR1_out)))
      
      
      ##! clean-up
      system2(command = "rm", args = c("-rf", dir_tmp))
    })
    
  } else {
    stop('Invalid argument provided for option "paired". Accepted values are "TRUE" or "FALSE"')
    
  }
}


amplicon_clip <- function(list_in, is_paired, file_ext, forA, revA,
                          max_err, min_len, ncore, list_out) {
  
  #! load packages
  suppressMessages(suppressWarnings(library(ShortRead)))
  suppressMessages(suppressWarnings(library(future.apply)))
  
  
  
  #! remove adapters/primers
  if (is_paired) {     # paired-end
    plan(multicore, workers = ncore)
    future_sapply(seq_along(list_in), function(i) {
      
      ##! set variables
      iR1 <- list_in[i]
      iR2 <- gsub("R1", "R2", iR1)
      iR1_out <- list_out[i]
      iR2_out <- gsub("R1", "R2", iR1_out)
      iR_stdout <- gsub(file_ext, "_log.txt", iR1_out)
      
      
      # ##! remove linked-adapters in sequences of both + and - strand
      # #! use with framed amplicons where both adapter have to be present in sequences (e.g. V4 16S or 18S)
      # system2(command = "cutadapt", args = c(
      #   paste0("-g ", forA, "...", reverseComplement(DNAString(revA))),
      #   paste0("-g ", revA, "...", reverseComplement(DNAString(forA))),
      #   paste0("-G ", forA, "...", reverseComplement(DNAString(revA))),
      #   paste0("-G ", revA, "...", reverseComplement(DNAString(forA))),
      #   paste("-e", max_err), paste("--minimum-length", min_len),
      #   "--overlap 10", "--discard-untrimmed", "--pair-filter=any", "--max-n 0", 
      #   paste("-o", iR1_out), paste("-p", iR2_out), iR1, iR2),
      #   stdout = iR_stdout)
      
      
      ##! remove adapters in sequences of both + and - strand
      #! use with amplicons that might be longer than sequences, rev adapter might be missing (e.g. ITS)
      system2(command = "cutadapt", args = c(
        paste("-g", forA, "-a", reverseComplement(DNAString(revA))),
        paste("-g", revA, "-a", reverseComplement(DNAString(forA))),
        paste("-G", forA, "-A", reverseComplement(DNAString(revA))),
        paste("-G", revA, "-A", reverseComplement(DNAString(forA))),
        paste("-e", max_err), paste("--minimum-length", min_len),
        "--overlap 10", "--discard-untrimmed", "--pair-filter=any", "--max-n 0", 
        paste("-o", iR1_out), paste("-p", iR2_out), iR1, iR2),
        stdout = iR_stdout)
    })
    
  } else if (!is_paired) { # single-end
    plan(multicore, workers = ncore)
    future_sapply(seq_along(list_in), function(i) {
      
      ##! set variables
      iR1 <- list_in[i]
      iR1_out <- list_out[i]
      iR_stdout <- gsub(file_ext, "_log.txt", iR1_out)
      
      
      # ##! remove linked-adapters in sequences of both + and - strand
      # #! use with framed amplicons where both adapter have to be present in sequences (e.g. V4 16S or 18S)
      # system2(command = "cutadapt", args = c(
      #   paste0("-g ", forA, "...", reverseComplement(DNAString(revA))),
      #   paste0("-g ", revA, "...", reverseComplement(DNAString(forA))),
      #   paste("-e", max_err), paste("--minimum-length", min_len),
      #   "--overlap 10", "--discard-untrimmed", "--max-n 0", 
      #   paste("-o", iR1_out), iR1),
      #   stdout = iR_stdout)
      
      
      ##! remove adapters in sequences of both + and - strand
      #! use with amplicons that might be longer than sequences, rev adapter might be missing (e.g. ITS)
      system2(command = "cutadapt", args = c(
        paste("-g", forA, "-a", reverseComplement(DNAString(revA))),
        paste("-g", revA, "-a", reverseComplement(DNAString(forA))),
        paste("-e", max_err), paste("--minimum-length", min_len),
        "--overlap 10", "--discard-untrimmed", "--max-n 0", 
        paste("-o", iR1_out), iR1),
        stdout = iR_stdout)
    })
    
  } else {
    stop('Invalid argument provided for option "paired". Accepted values are "TRUE" or "FALSE"')
    
  }
}


amplicon_parsing <- function(dir_raw, parse_modules, is_paired, file_ext, 
                             forA, revA, max_err, min_len, ncore, dir_parsed, parsing_token) {
  
  #! set variables
  list_raw <- list.files(dir_raw, pattern = paste0(file_ext,"$"),
                         full.names = T, recursive = T)
  list_raw <- list_raw[!grepl("unknown", list_raw)] # remove reads that did NOT match to any known barcode
  
  dir_reorient <- file.path(dir_parsed, "reoriented")
  list_reorient <- file.path(dir_reorient, basename(list_raw))
  
  dir_clip <- file.path(dir_parsed, "clipped")
  list_clip <- file.path(dir_clip, basename(list_raw))
  
  
  
  #! run parsing
  if (identical(c("reorient", "clip"), parse_modules)) {
    
    cat("\n### RE-ORIENT READS\n")
    dir.create(dir_reorient, recursive = T, showWarnings = F)
    amplicon_reorient(list_in = list_raw, is_paired = is_paired, file_ext = file_ext, 
                      forA = forA, revA = revA, 
                     ncore = ncore, list_out = list_reorient)
    
    cat("\n### CLIP READS\n")
    dir.create(dir_clip, recursive = T, showWarnings = F)
    amplicon_clip(list_in = list_reorient, is_paired = is_paired, file_ext = file_ext,
                  forA = forA, revA = revA, max_err = max_err, min_len = min_len, 
                  ncore = ncore, list_out = list_clip)

    
  } else if (identical(c("clip"), parse_modules)) {
    
    cat("\n### CLIP READS\n")
    dir.create(dir_clip, recursive = T, showWarnings = F)
    amplicon_clip(list_in = list_raw, is_paired = is_paired, file_ext = file_ext,
                  forA = forA, revA = revA, max_err = max_err, min_len = min_len, 
                  ncore = ncore, list_out = list_clip)
    
    
  } else {
    stop(paste0('\n\nInvalid parsing module/s provided. Available modules include:\n',
         '(1) reorient, clip\n', '(2) clip\n\n'))
   
     
  }
  
  system2(command = "touch", args = parsing_token)
  
}


amplicon_parsing(dir_raw = snakemake@input[["dir_raw"]],
                 parse_modules = snakemake@params[["parse_modules"]],
                 is_paired = snakemake@params[["is_paired"]],
                 file_ext = snakemake@params[["file_ext"]],
                 forA = snakemake@params[["forA"]],
                 revA = snakemake@params[["revA"]],
                 max_err = snakemake@params[["max_err"]],
                 min_len = snakemake@params[["min_len"]],
                 ncore = snakemake@params[["ncore"]],
                 dir_parsed = snakemake@output[["dir_parsed"]],
                 parsing_token = snakemake@output[["parsing_token"]]
)







####### OLD CODE still to be implemented:


# ### (optional) PacBio make CCS & demultiplex ----
# 
# ## make CCS
# mkdir 0_ccs
# 
# # conda install -c bioconda pbcss==3.4 # pbcss=3.3 for 2017 chemistry
# find 0_raw -type f -name "*.subreads.bam" | while read file; do ccs -j 20 --minPredictedAccuracy 0.99 $file 0_ccs/$(basename $file | sed 's/\.bam/\.fastq/'); done # for css 3.3
# gzip 0_ccs/*.fastq
# 
# # conda install -c bioconda pbcss # for >2017 chemistry
# find 0_raw -type f -name "*.subreads.bam" | while read file; do ccs -j 20 --min-rq 0.99 0_ccs/$file $(basename $file | sed 's/\.bam/\.fastq.gz/'); done
# 
# 
# ## demultiplex
# mkdir 0_demux
# cat 0_ccs/*.fastq.gz > 0_ccs/ccs_acc099.fastq.gz
# 
# # conda install -c bioconda lima
# lima -j 20 --same --ccs --min-score 80 --split 0_ccs/ccs_acc099.fastq.gz barcodes.fa 1_demultiplexed/demux.fastq.gz
# 
# 
# 
# ### (optional) convert BAM to fastq files ----
# x = list.files(dir_raw, pattern = ".bam$", full.names = T, recursive = T)
# 
# sapply(x, function(i) {
#   i2 = gsub(".bam$", ".fq", i)
#   system(paste("samtools bam2fq -@ 7", i, ">", i2))
#   system(paste("bzip2 -z", i2))
# })
# 
# 
# 
# ### (optional) merge input ----
# 
# x = list.files("0_CSS/", pattern = ".fastq.gz$", full.names = T, recursive = T)
# 
# system(paste("cat", paste(x, collapse = " "), "> 0_CSS/css_acc999.fastq.bz2"))
# 
# 
# 
# ### (optional) merge fasta and qual files ----
# #!!! for mrDNA use the windows app they made
# x = list.files(dir_raw, pattern = ".fasta", full.names = T, recursive = T)
# 
# library(future.apply); plan(multicore, workers = ncore) # new (easier) parallelization
# future_sapply(x, function(i) {
#   i2 = gsub(".fasta", "", i)
#   system(paste0("reformat.sh in=", i2, ".fasta qfin=", i2, ".qual out=", i2, ".fq vint=t overwrite=t"))
#   # system(paste0("~/myscript/python/convert_fastaqual_fastq.py -f ", i2, ".fasta -q ", i2, ".qual -o 0_raw")) # tweak for mrDNA files
# })
# 
# 
# 
# ### (optional) demultiplexing - with barcodes ----
# barc_file = "barcodes.fa"
# dir.create(dir_demux, recursive = T)
# 
# #!!! single-end
# x = list.files(dir_raw, pattern = ".fastq.bz2$", full.names = T, recursive = T)
# sapply(x, function(i) {
#   system2(command = "cutadapt", args = c('-e 0.15', '--no-indels', paste0('--cores=', ncore), 
#                                          paste0('-g file:', barc_file), "--overlap 16", # adjust to barcode length
#                                          paste0('-o ', file.path(dir_demux, '{name}.fastq.gz')), i),
#           stdout = file.path(dir_demux, gsub(".fastq.gz$", ".log", basename(i))))
# })
# 
# #!!! paired-end ???still to code???
# x = list.files(dir_raw, pattern = "R1.fastq.bz$", full.names = T, recursive = T)
# sapply(x, function(i) {
#   i_R2 = gsub("R1", "R2", i)
#   system2(command = "cutadapt", args = c('-e 0.15', '--no-indels', paste0('--cores=', ncore), 
#                                          paste0('-g file:', barc_file), "--overlap 16", # adjust to barcode length
#                                          paste0('-o ', file.path(dir_demux, '{name}.1.fastq.gz')),
#                                          paste0('-p ', file.path(dir_demux, '{name}.2.fastq.gz')), 
#                                          i, i_R2),
#           stdout = file.path(dir_demux, gsub(".fastq.gz$", ".log", basename(i))))
# })
# 
# 
# 
# ### (optional) demultiplexing - with header names ----
# dir.create(dir_demux, recursive = T)
# 
# #!!! single-end
# x = list.files(dir_raw, pattern = ".fastq.bz2$", full.names = T, recursive = T)
# sapply(x, function(i) {
#   system2(command = "demuxbyname.sh", args = c(paste0("in=", i),
#                                                paste0("out=", file.path(dir_demux, '%.fastq.gz')),
#                                                paste0("outu=", file.path(dir_demux, 'unknown.fastq.gz')),
#                                                paste0("stats=", file.path(dir_demux, gsub(".fastq.gz$", ".log", basename(i)))),
#                                                "prefixmode=t", "delimiter=::", paste0("threads=", ncore)),
#           stdout = F)
# })
# 
# #!!! paired-end ???still to code???
# x = list.files(dir_raw, pattern = "R1.fastq.bz$", full.names = T, recursive = T)
# sapply(x, function(i) {
#   i_R2 = gsub("R1", "R2", i)
#   system2(command = "cutadapt", args = c('-e 0.15', '--no-indels', paste0('--cores=', ncore), 
#                                          paste0('-g file:', barc_file), "--overlap 16", # adjust to barcode length
#                                          paste0('-o ', file.path(dir_demux, '{name}.1.fastq.gz')),
#                                          paste0('-p ', file.path(dir_demux, '{name}.2.fastq.gz')), 
#                                          i, i_R2),
#           stdout = file.path(dir_demux, gsub(".fastq.gz$", ".log", basename(i))))
# })
# 
# 
# 
# #!!! check cutadapt output
# x_log1 = list.files(dir_clip, pattern = "_log.txt$", full.names = T, recursive = T)
# x_log2 = lapply(x_log1, function(i) {
#   i2 = readLines(i)
#   i3 = c(i2[grepl("Total read pairs processed:", i2)],
#          i2[grepl("Pairs that were too short:", i2)],
#          i2[grepl("Pairs with too many N:", i2)],
#          i2[grepl("Pairs written \\(passing filters\\):", i2)])
#   c(gsub("_log.txt", "", basename(i)), readr::parse_number(i3))
# })
# x_log3 = as.data.frame(do.call(rbind, x_log2))
# x_log3[, 2:5] = apply(x_log3[, 2:5], 2, as.numeric)
# x_log3$eff = x_log3[, 5]/x_log3[, 2]
# x_log3$no_primer = x_log3[, 2] - x_log3[, 5]
# x_log3 = x_log3[, c(1:2,7,3:6)]
# colnames(x_log3) = c("Sample", "Raw", "no_primer", "too_short", "with_N", "Trimmed", "efficiency")
# x_log = x_log3
# write.table(x_log, "clipping_stats.txt", 
#             col.names = T, row.names = F, quote = F, sep = "\t")