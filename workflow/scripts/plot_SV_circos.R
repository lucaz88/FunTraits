# #___ for testing
# setwd("~/VSC_data/lucaz/genomics/G00003_Insect_cell_lines/")
# ref_fai <- "_input/genomes/ref/tni_chr_v1.0.fai"
# ref_gff <- "_input/genomes/ref/tni_gene_v1.gff3"
# target_SV_svim <- "_results/SV_svim/FNL.1"
# in_ref_annotation <- "_input/genomes/ref/tni_annotation_v1.csv"
# plot_file <- "_plots/circos_SV.pdf"
# min_len <- 1000
# #___ for testing



plot_SV_circos <- function(ref_fai, ref_gff, target_SV_svim, # in_ref_annotation
                     plot_file,
                     min_len) {
  
  #! load packages
  suppressMessages(suppressWarnings(library(circlize)))
  suppressMessages(suppressWarnings(library(stringr)))
  suppressMessages(suppressWarnings(library(RColorBrewer)))
  suppressMessages(suppressWarnings(library(rtracklayer)))
  
  
  #! parse ref fai faile
  ref_fai_tab <- read.table(ref_fai, header = F)
  names(ref_fai_tab) <- c("NAME","LENGTH","OFFSET","LINEBASES","LINEWIDTH")
  ref_fai_tab <- ref_fai_tab[order(as.integer(gsub("chr", "", ref_fai_tab$NAME, ignore.case = T))), ]
  ref_fai_tab$NAME <- factor(ref_fai_tab$NAME, levels = unique(ref_fai_tab$NAME))
  
  
  #! parse ref gff file
  ref_gff <- as.data.frame(readGFF(ref_gff))
  ref_gff_exon <- ref_gff[ref_gff$type == "exon", c("seqid","start","end","ID")]
  # ref_ann <- read.delim(in_ref_annotation, h=F)
  # ref_gff_exon$ID <- gsub("-.*", "", ref_gff_exon$ID)
  # ref_gff_exon$annotation <- ref_ann$V2[match(ref_gff_exon$ID, ref_ann$V1)]
  
  
  #! parse svf file
  svf_file <- read.table(file.path(target_SV_svim, "variants.vcf"), header = F)
  svf_file <- svf_file[svf_file$V7 == "PASS", ]
  svf_tab <- data.frame(chr=svf_file$V1,
                       start=abs(as.numeric(svf_file$V2)),
                       end=abs(as.numeric(str_match(svf_file$V8, "END=\\s*(.*?)\\s*;")[, 2])),
                       value=abs(as.numeric(str_match(svf_file$V8, "SVLEN=\\s*(.*?)\\s*;")[, 2])),
                       type=str_match(svf_file$V8, "SVTYPE=\\s*(.*?)\\s*;")[, 2],
                       seq=apply(svf_file[, c("V4","V5")], 1, function(z) z[which.max(nchar(z))]))
  svf_tab <- svf_tab[svf_tab$type != "BND", ] # complex rearrangements with breakends
  svf_tab$value[svf_tab$type == "INV"] <- svf_tab$end[svf_tab$type == "INV"] - svf_tab$start[svf_tab$type == "INV"]
  svf_tab$type_col <- brewer.pal(n = length(unique(svf_tab$type)), name = "Set1")[match(svf_tab$type, unique(svf_tab$type))]
  svf_tab <- svf_tab[svf_tab$value >= min_len, ]
  svf_tab_list <- split(svf_tab, f = factor(svf_tab$type, levels = c("DEL","INS","DUP:INT","DUP:TANDEM","INV")))
  svf_tab_list <- svf_tab_list[lapply(svf_tab_list, nrow)>0] # remove missing type of SV
  
  
  #! Plotting
  
  ##! initialize plot
  pdf(file=plot_file,
      width=20,
      height=20,
      pointsize = 35
  )
  # dotrow_h = 0.1
  chr_gap <- 0
  circos.par("cell.padding" = c(0, 0),
             "track.margin" = c(0, 0),
             "gap.degree" = chr_gap,
             "start.degree" = 90,
             "track.height" = 0.05
  )
  circos.initialize(sectors = ref_fai_tab$NAME,
                    xlim = as.matrix(cbind(start = 0,
                                           end = ref_fai_tab$LENGTH))
  )
  
  
  ##! add labels
  circos.track(ylim = c(0,1), 
               bg.col = "grey", bg.border = NA, 
               panel.fun = function(x, y) {
                 circos.text(x = CELL_META$xcenter, y = CELL_META$ycenter, 
                             labels = CELL_META$sector.index, cex = 0.65)
               })
  
  
  ##! add axis
  brk <- seq(from=0, to=max(ref_fai_tab$LENGTH), by=25e5)
  circos.track(track.index = get.current.track.index(), 
               bg.border=F,
               panel.fun=function(x, y) {
                 circos.axis(h="top", 
                             major.at=brk, major.tick.length = .5,
                             labels=c(NA, paste0(brk/1e6, "")[-1]), # Mb
                             labels.cex=0.55, 
                             labels.pos.adjust = F, # col=col_text, labels.col=col_text, 
                             lwd=0.7, 
                             labels.facing="clockwise") # clockwise inside
               })
  
  
  ##! add gene density
  circos.genomicDensity(ref_gff_exon,
                        col = c("#4f5b5d"),
                        track.height = 0.1)
  
  
  ##! add SV tracks
  #! using dots
  # lapply(svf_tab_list, function(iSV_type) {
  #   circos.genomicTrack(data=iSV_type, 
  #                       track.height=0.1, 
  #                       bg.border=T,
  #                       panel.fun=function(region, value, ...) {
  #                         circos.genomicPoints(region, 
  #                                              value, 
  #                                              col = unique(iSV_type$type_col),
  #                                              pch = 15)
  #                         # circos.segments(x0=CELL_META$xlim[1], 
  #                         #                 x1=CELL_META$xlim[2], 
  #                         #                 y0=min(iSV_type$value), 
  #                         #                 y1=min(iSV_type$value), 
  #                         #                 lwd=2, lty="11", col="#4f5b5d")
  #                         # circos.segments(x0=CELL_META$xlim[1], 
  #                         #                 x1=CELL_META$xlim[2], 
  #                         #                 y0=max(iSV_type$value), 
  #                         #                 y1=max(iSV_type$value), 
  #                         #                 lwd=2, lty="11", col="#4f5b5d")
  #                       })
  #   circos.yaxis(at=c(min(iSV_type$value), max(iSV_type$value)), 
  #                labels.cex = .3,
  #                sector.index = "chr0"
  #                )
  # })
  
  # #!!! useing bands
  # lapply(svf_tab_list[1], function(iSV_type) {
  #   circos.genomicTrack(data=iSV_type,
  #                       ylim = c(0, 1),
  #                       track.height=0.1,
  #                       bg.border=NA,
  #                       bg.col="#eeeeee",
  #                       panel.fun=function(region, value, ...) {
  #                         circos.genomicRect(region, value,
  #                                            ybottom = 0, ytop = 1,
  #                                            col = unique(iSV_type$type_col),
  #                                            border = NA
  #                         )
  #                       })
  # })
  
  #! using density line
  lapply(svf_tab_list, function(iSV_type) {
    circos.genomicDensity(iSV_type, 
                          count_by = "percent", # "percent", "number"
                          col = unique(iSV_type$type_col), 
                          track.height = 0.1)
    circos.yaxis(at=c(min(svf_tab$value), max(svf_tab$value)), 
                 labels.cex = .3,
                 sector.index = "chr0"
    )
  })
  
  
  # ##! add rearrangements
  # rcols <- scales::alpha(ifelse(sign(nuc1$start - nuc1$end) != sign(nuc2$start - nuc2$end), "#f46d43", "#66c2a5"), alpha=0.4)
  # circos.genomicLink(nuc1, nuc2, col=rcols, border=NA)
  
  
  ##! legend
  legend(x = "bottomleft", horiz=F,
         legend=svf_tab$type[!duplicated(svf_tab$type)], 
         col = svf_tab$type_col[!duplicated(svf_tab$type)],
         pch=15, pt.cex = 2, bty="n", cex=.8) # inset = c(0, .05)
  legend(x = "center", horiz=F,
         legend=strsplit(target_SV_svim, "/")[[1]][3], 
         col = NA, cex=1.3)
  
  
  ##! close plot
  circos.clear()
  dev.off()
}


plot_SV_circos(ref_fai = snakemake@input[["ref_fai"]], 
          ref_gff = snakemake@input[["ref_gff"]],
          target_SV_svim = snakemake@input[["target_SV_svim"]], 
          # in_ref_annotation = snakemake@input[["ref_annotation"]], #! not yet implemented in the plotting
          plot_file = snakemake@output[["plot_file"]],
          min_len = snakemake@params[["min_len"]]
)
