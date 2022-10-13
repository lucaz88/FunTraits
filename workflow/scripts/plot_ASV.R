### downs analyses ----
???still to be coded!!!
dir.create(dir_res, recursive = T)

library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2"); library(ggrepel)
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")



### preparing data

## metadata
smpl_meta = read.delim("../Atacama_meta.tsv", header = T, check.names = F)
row.names(smpl_meta) = smpl_meta$label
smpl_meta$cols = palette.colors(n = length(unique(smpl_meta$Site)))[match(smpl_meta$Site, unique(smpl_meta$Site)[c(3:5,1,2,3)])]


## ASV table, sequences and Taxonomy
# ASVtab = read.delim(file.path(dir_tab, "ASV_tab_wTAX_nocontam.tsv"), h=T)
ASVtab = read.delim(file.path(dir_tab, "ASV_tab_wTAX.tsv"), h=T)

seqs = Biostrings::DNAStringSet(ASVtab$seq)
names(seqs) = ASVtab$ASVs


#! for SSU (6-level based taxonomy)
tax_tab = data.frame(t(as.data.frame(lapply(strsplit(ASVtab$Taxonomy, ";"), function(x) x[1:6]))),
                      row.names = ASVtab$ASVs)
colnames(tax_tab) = c("Kingdom","Phylum","Class","Order","Family","Genus")
# #! for LSU (variable levels)
# max_rank = max(stringr::str_count(ASVtab$Taxonomy, ";"))
# tax_tab = data.frame(t(as.data.frame(lapply(strsplit(ASVtab$Taxonomy, ";"), function(x) x[1:max_rank]))),
#                       row.names = ASVtab$ASVs)
# colnames(tax_tab) = paste0("Rank_", 1:max_rank)
# #! add II taxonomy
# tax_tab = cbind(tax_tab, ASVtab[, (which(colnames(ASVtab) == "Taxonomy")+1):ncol(ASVtab)])


ASVtab = ASVtab[, -2] # drop seq info
row.names(ASVtab) = ASVtab$ASVs # paste0(ASVtab$ASVs, "-", ASVtab$Taxonomy)
ASVtab = ASVtab[, -c(1, ncol(ASVtab))]
colnames(ASVtab) = row.names(smpl_meta)[match(gsub("AD.|\\..*", "", colnames(ASVtab)),
                                              smpl_meta$Tubenummer)]
ASVtab = ASVtab[, match(row.names(smpl_meta), colnames(ASVtab))]
dim(ASVtab); sum(ASVtab)


## sequences tree
# library(DECIPHER); library(phangorn)
# alignment = AlignSeqs(seqs, anchor=NA, processors = 7)
# phang.align = phyDat(as(alignment, "matrix"), type="DNA")
# dm = dist.ml(phang.align)
# treeNJ = NJ(dm) # Note, tip order != sequence order
# fit = pml(treeNJ, data=phang.align)
# fitGTR = update(fit, k=4, inv=0.2)
# fitGTR = optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                     rearrangement = "stochastic", control = pml.control(trace = 0))
# saveRDS(fitGTR, file.path(dir_tab, "fitGTR.RDS"))
fitGTR = readRDS(file.path(dir_tab, "fitGTR.RDS"))
my_tree = fitGTR$tree
#! root tree using the longest branch with a tip
pick_new_outgroup = function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT =
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup = treeDT[which.max(length)]$id
  return(new.outgroup) }
out.group = pick_new_outgroup(phy_tree(my_tree))
my_tree = ape::root(my_tree, outgroup=out.group, resolve.root=T)



### create phyloseq object
(ps = phyloseq(tax_table(as.matrix(tax_tab)), sample_data(smpl_meta),
                otu_table(ASVtab, taxa_are_rows = T)))
# (ps = phyloseq(tax_table(as.matrix(tax_tab)), sample_data(smpl_meta),
#                otu_table(ASVtab, taxa_are_rows = T), phy_tree(my_tree)))



### filtering and transformations

## taxonomic filt
(ps = subset_taxa(ps, !is.na(Phylum)))
(ps = subset_taxa(ps, Order != "Chloroplast"))
# (ps = subset_taxa(ps, !is.na(Rank_1)))
# (ps = subset_taxa(ps, Rank_1 != "uncl."))


## prevalence filt
(prevalenceThreshold = 0.1 * nsamples(ps))
keepTaxa = apply(otu_table(ps), 1, function(x) sum(x > 0)) 
(ps = prune_taxa(names(keepTaxa)[keepTaxa > prevalenceThreshold], ps))


## metadata filt
(ps = subset_samples(ps, is_blank == F & Site == "PP"))
# (ps = subset_samples(ps, is_blank == F))
# (ps = subset_samples(ps, is_blank == F & Site != "PP"))



### transformations
ps_trans = ps


## variance stabilization
# deseq_counts = DESeqDataSetFromMatrix(ps@otu_table@.Data, colData = smpl_meta, design = ~label) 
deseq_counts = phyloseq_to_deseq2(ps, design = ~label) 
deseq_counts_vst = varianceStabilizingTransformation(deseq_counts, fitType="mean") # mean parametric
vst_trans_count_tab = assay(deseq_counts_vst)
vst_trans_count_tab[vst_trans_count_tab < 0] = 0 # remove negative values as they were close to 0 in the raw data
ps_trans@otu_table@.Data = vst_trans_count_tab


## rel. abundance transformation
ps_trans@otu_table@.Data = decostand(ps_trans@otu_table@.Data, method = "total", MARGIN = 2)
write.table(data.frame(ASV=row.names(ps_trans@otu_table@.Data), ps_trans@otu_table@.Data, 
                       Taxonomy=apply(ps_trans@tax_table@.Data, 1, paste0, collapse = ";")),
            file.path(dir_tab, "ASV_table_wTAX_PP_std.tsv"), 
            col.names = T, row.names = F, quote = F, sep = "\t")


## clr transformation
ps_trans = microbiome::transform(ps_trans, transform = "clr")
write.table(data.frame(ASV=row.names(ps_trans@otu_table@.Data), ps_trans@otu_table@.Data, 
                       Taxonomy=apply(ps_trans@tax_table@.Data, 1, paste0, collapse = ";")),
            file.path(dir_tab, "ASV_table_wTAX_PP_clr.tsv"), 
            col.names = T, row.names = F, quote = F, sep = "\t")



# ### alpha-deiversity
# rarecurve(t(ASVtab), step=100, col=smpl_meta$cols, lwd=2, ylab="ASVs", label=F)
# # and adding a vertical line at the fewest seqs in any sample
# abline(v=(min(rowSums(t(ASVtab)))))
# 
# count_tab_phy = otu_table(ASVtab, taxa_are_rows=T)
# tax_tab_phy = tax_tab # tax_table(tax_tab)
# ASV_physeq = phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
# 
# p = plot_richness(ASV_physeq, x="Site", color="Depth", measures=c("Chao1", "Shannon")) + 
#   # scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) +
#   theme(legend.title = element_blank())
# # Chao1 is richness and Shannon is diversity



### Taxonomy plots

#! to highlight main bacterial taxa
tax_group = ps@tax_table@.Data[, 2]
tax_group[grepl("Proteobacteria", ps@tax_table@.Data[, 2])] = ps@tax_table@.Data[grepl("Proteobacteria", ps@tax_table@.Data[, 2]), 3]

# #! to highlight Fungi
# tax_group = ps@tax_table@.Data[, 7]
# tax_group[grepl("Dikarya", ps@tax_table@.Data[, 7])] = ps@tax_table@.Data[grepl("Dikarya", ps@tax_table@.Data[, 7]), 8]

# #! to  highlight Chytrid specificity
# tax_group = sapply(strsplit(ps@tax_table@.Data[, 19], ";"), "[[", 2)
# tax_group[tax_group == "NA"] = "unknown" 
# tax_group[as.numeric(ps@tax_table@.Data[, 20]) < 99] = "unknown" # filter on blast identity


major_taxa_counts_tab = rowsum(otu_table(ps), group = tax_group)
row.names(major_taxa_counts_tab)[is.na(row.names(major_taxa_counts_tab))] = "uncl."
major_taxa_proportions_tab = apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)
dim(major_taxa_proportions_tab)
# reduce taxa by keeping those that make up >5% in any individual sample; discarded will be lumped as 'Others"
temp_filt_major_taxa_proportions_tab = data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 5, ])
filtered_proportions = colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
major_taxa_proportions_tab = rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)
dim(major_taxa_proportions_tab)
plot_table = melt(as.matrix(major_taxa_proportions_tab))


## aesthetic adds
plot_table$value[plot_table$value == 0] = NA
plot_table$value[plot_table$X1 %in% c("Ankyra","Breviata anathema","Choanoflagellata","Costaria",
                                      "Metazoa","uncl.","Other")] = NA
plot_table$lifestyle = ps@sam_data$lifestyle[match(plot_table$X2, ps@sam_data$label)]
plot_table = plot_table[plot_table$lifestyle == "PA", ]
plot_table$lake_chain = ps@sam_data$lake_chain[match(plot_table$X2, ps@sam_data$label)]
plot_table$lake_chain = factor(plot_table$lake_chain, levels = unique(plot_table$lake_chain)[c(2,3,1,4)])
plot_table$lake_chain = mapvalues(plot_table$lake_chain, from = levels(plot_table$lake_chain),
                                  to = c("High-conn","Mid-conn","Low-conn","No-conn"))
plot_table$month = ps@sam_data$month[match(plot_table$X2, ps@sam_data$label)]
plot_table$month = factor(plot_table$month, levels = unique(plot_table$month)[c(2,1,3)])
plot_table$month = mapvalues(plot_table$month, from = levels(plot_table$month),
                             to = c("March","April","May"))
plot_table$lake = ps@sam_data$lake[match(plot_table$X2, ps@sam_data$label)]


## stacked bars
pdf(file=file.path(dir_res, "Bar_plot_fungi.pdf"), width=9, height=5) # font size 10
ggplot(plot_table, aes(x = lake, y = value, fill = X1)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(month~lake_chain, scales = "free_x", space = "free") +
  labs(x="", y="Relative abundance (%)", fill="Main taxa") +
  scale_x_discrete(limits = rev) +
  scale_fill_brewer(type = "qual") +
  # coord_flip() +
  theme_minimal() +
  theme(legend.position = "right", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  ) 
dev.off()


# ## bubble grid
# ggplot(plot_table, aes(x = Var2, y = Var1, colour = Var1, size = value)) +
#   geom_point() +
#   # geom_text(aes(label = value), colour = "white", size = 3) +
#   scale_x_discrete(position = "bottom", ) +
#   # scale_size_continuous(range = c(10, 30)) + # Adjust as required.
#   # scale_color_brewer(palette = "Set2") +
#   scale_size_area() +
#   labs(x = NULL, y = NULL) +
#   theme(legend.position = "right",
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         axis.ticks = element_blank()) + guides(color = F)



### Beta-div
d.meth = "wunifrac" # unifrac wunifrac bray euclidean
vst_dist = phyloseq::distance(ps_trans, method = d.meth) 

## Hierarchical clustering
vst_clust = hclust(vst_dist, method="ward.D2") # ward.D2 complete
vst_dend = as.dendrogram(vst_clust, hang=0.1)
labels_colors(vst_dend) = get_variable(ps_trans, "cols")[order.dendrogram(vst_dend)]
plot(vst_dend, ylab=paste0("VST - ", dmeth))


## Hierarchical clustering - with bootstrap
library(pvclust)
vst_clust = pvclust(data = otu_table(ps_trans), 
                    method.dist=function(x) {
                      x2 = prune_taxa(row.names(x), ps_trans)
                      phyloseq::distance(x2, method = d.meth)
                    }, 
                    method.hclust="ward.D2", nboot=100, 
                    r = seq(from = 0.4, to = 1, length.out = 6), parallel=F)
plot(vst_clust, ylab=d.meth, main="")
pvrect(vst_clust, alpha=0.9)


## PCoA
vst_pcoa = ordinate(ps_trans, method="MDS", distance=vst_dist)
eigen_vals = vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

p = plot_ordination(ps_trans, vst_pcoa, color="Site") + 
  geom_point(size=3) + labs(col="type") + 
  geom_text_repel(aes(label=get_variable(ps_trans, "label")), direction = "both", force = 2) +
  # geom_text(aes(label=rownames(smpl_meta), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=get_variable(ps_trans, "cols"), 
                     limits=get_variable(ps_trans, "Site")) + 
  theme(legend.position="none")
p$coordinates$clip = "off"
p


## dbRDA
ps_trans@sam_data = decostand(ps_trans@sam_data, method = "normalize", MARGIN = 2) # use [, c()] to normalize only specific columns
test_var = as.formula(paste0("~ ", paste(names(ps_trans@sam_data)[c(1,3,5,6)], collapse = "+"))) # select metadata variables that you want to include in the dbRDA 

vst_dbrda = ordinate(ps_trans, "CAP", distance=vst_dist, formula = test_var)
eigen_vals = vst_dbrda$CCA$eig # allows us to scale the axes according to their magnitude of separating apart the samples
arrowmat = vegan::scores(vst_dbrda, display = "bp")
arrowdf = data.frame(labels = rownames(arrowmat), arrowmat)

p = plot_ordination(ps_trans, vst_dbrda, color="Site") +
  geom_point(size=3) + labs(col="type") +
  geom_text_repel(aes(label=get_variable(ps_trans, "label")), direction = "both", force = 2) +
  # geom_text(aes(label=rownames(smpl_meta), hjust=0.3, vjust=-0.4)) +
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") +
  geom_segment(mapping = aes(xend = CAP1, yend = CAP2, x = 0, y = 0,
                             shape = NULL, color = NULL),
               size = .5, data = arrowdf, color = "gray",
               arrow = arrow(length = unit(0.02, "npc"))) +
  geom_text(mapping = aes(x = 1.3 * CAP1, y = 1.3 * CAP2, shape = NULL,
                          color = NULL, label = labels),
            size = 4, data = arrowdf,show.legend = FALSE) +
  scale_color_manual(values=get_variable(ps_trans, "cols"),
                     limits=get_variable(ps_trans, "Site")) +
  theme(legend.position="none")
p$coordinates$clip = "off"
p

#! model stats
anova(vst_dbrda)
anova(vst_dbrda, by="axis", perm.max=1000)
anova(vst_dbrda, by="terms", perm.max=1000)



# ### heatmap
# library(ComplexHeatmap)
# library(vegan)
# 
# seqtab.nochim.wtax = seqtab.nochim.wtax
# seqtab.nochim.wtax = SpiecEasi::clr(seqtab.nochim.wtax)
# 
# row.names(seqtab.nochim.wtax) = sapply(strsplit(row.names(seqtab.nochim.wtax), ";"), function(i) paste(i[c(1, 6:length(i))], collapse = ";"))
# row.names(seqtab.nochim.wtax) = gsub("Eukaryota;", "", row.names(seqtab.nochim.wtax))
# smpl_map = as.data.frame(readxl::read_excel("../AC18 Proben 2018 DNA IGB.xls", skip = 2))
# colnames(seqtab.nochim.wtax) = smpl_map$`AC 18 Proben`[match(gsub("AD.|\\..*", "", colnames(seqtab.nochim.wtax)), smpl_map$Tubenummer)]
# 
# # write.table(data.frame(ASV=row.names(seqtab.nochim.wtax), seqtab.nochim.wtax),
# #             file = "18S-ASVtab-decontam-15reads-trans.txt", 
# #             quote = F, row.names = F, col.names = T, sep = "\t")
# 
# pdf(file="18S-hm-decontam-15reads.pdf", width=18, height=20) # font size 10
# draw(Heatmap(seqtab.nochim.wtax,
#              cluster_rows = hclust(vegdist(seqtab.nochim.wtax, "eucl"), "complete"), # eucl for clr data (Quinn et al. 2018)
#              row_dend_width = unit(.1, "npc"),
#              cluster_columns = hclust(vegdist(t(seqtab.nochim.wtax), "eucl"), "complete"),
#              column_dend_height = unit(.1, "npc"), 
#              row_names_max_width = unit(.7, "npc")),
#      heatmap_legend_side = "left")
# dev.off()