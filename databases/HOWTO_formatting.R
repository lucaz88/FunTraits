### PR2 databases

#! version <= 4.12
raw_pr2 = system("cut -f2,20,21,22,23,24,25,26,34 /media/lucaz/DATA/DBs_repository/PR2/pr2_version_4.12.0_merged.tsv", intern = T)
raw_pr2 = strsplit(raw_pr2, "\t")
raw_pr2 = do.call(rbind, raw_pr2)
colnames(raw_pr2) = raw_pr2[1, ]
raw_pr2 = as.data.frame(raw_pr2[-1, ])
write.table(cbind(paste0(">", raw_pr2$pr2_accession, " ",
                         apply(raw_pr2[, 2:8], 1, paste0, collapse=";")),
                  raw_pr2$sequence),
            file = "/media/lucaz/DATA/DBs_repository/PR2/pr2_version_4.12.0_merged.fasta",
            col.names=F, row.names=F, sep="\n", quote=F)

#! version > 4.12
raw_pr2 = as.data.frame(readxl::read_xlsx("/media/lucaz/DATA/DBs_repository/PR2/pr2_version_4.13.0_merged.xlsx", col_names = T, col_types = "text"))
write.table(cbind(paste0(">", raw_pr2$pr2_accession, " ",
                         apply(raw_pr2[, 2:9], 1, paste0, collapse=";")),
                  raw_pr2$sequence),
            file = "/media/lucaz/DATA/DBs_repository/PR2/pr2_version_4.13.0_merged.fasta",
            col.names=F, row.names=F, sep="\n", quote=F)
system2(command = "makeblastdb",
        args = c("-in /media/lucaz/DATA/DBs_repository/PR2/pr2_version_4.13.0_merged.fasta",
                 "-dbtype nucl"))