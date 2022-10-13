#___ for testing
setwd("~/VSC_data/lucaz/genomics/G00003_Insect_cell_lines/")
in_gbk <- "_input/virus/baculovirus_AcMNPV-WIV-Syn1.gb"
ref_dir <- "_input/genomes"
ref_ext <- "*.fa$"
recursive <- T
out_dir <- "test_blast"
#___ for testing


tblastn_from_gbk <- function(in_gbk, ref_dir, ref_ext, recursive=F, out_dir) {
  
  #! load packages
  suppressMessages(suppressWarnings(library(read.gb)))
  suppressMessages(suppressWarnings(library(read.gb)))
  
  
  #! import data
  gbk_file <- read.gb(in_gbk, DNA = T, Type = "full")
  ref_files <- list.files(path = ref_dir, pattern = ref_ext, full.names = T, recursive = recursive)
  dir.create(out_dir, recursive = T)
  
  
  #! parse data
  faa_file <- unlist(extract.gb(gbk_file, "CDS")[[1]], recursive = F)
  faa_file <- lapply(faa_file, function(z) z[z[, 1] %in% c("CDS","protein_id","product","translation"), 2])
  faa_file <- as.data.frame(do.call(rbind, faa_file))
  colnames(faa_file) <- c("CDS","protein_id","product","translation")

  
  #! make query file
  query_file <- file.path(out_dir, "query.faa")
  write.table(cbind(paste0(">", faa_file$product),
                    faa_file$translation),
              query_file, col.names=F, row.names=F, sep="\n", quote=F) 
  
  
  #! blast
  ref_db <- file.path(out_dir, "refDB")
  blast_out <- file.path(out_dir, "blastout.tsv")
  if (!file.exists(file.path(out_dir, "refDB.ndb"))) {
    system2(command = "cat", args = c(ref_files, 
                                      paste0("> ", ref_db)))
    system2(command = "makeblastdb",
            args = c(paste("-in", ref_db),
                     "-dbtype nucl"))
    
    system2(command = "rm", args = ref_db)
  }
  system2(command = "tblastn", 
          args = c(paste("-query", query_file),
                   paste("-db", ref_db),
                   paste("-out", blast_out),
                   "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'",
                   "-evalue 1e-5", paste("-num_threads", ncore), "-max_target_seqs 10"
          ))
  
  
  ASV_tax <- read.delim(blast_out, sep = "\t", h=F)
  
}