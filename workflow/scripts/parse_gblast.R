parse_gblast = function(transp_folder, 
                        transp_set, 
                        BioV_transp_tab) {
  #! load packages
  suppressMessages(suppressWarnings(library(stringr)))
  suppressMessages(suppressWarnings(library(webchem)))
  
  #! import raw data
  transp_ann1 = lapply(transp_folder, function(i) 
    read.delim(file.path(i, "results.tsv"), h=T, sep = "\t"))
  names(transp_ann1) = sapply(strsplit(transp_folder, "/"), function(i) rev(i)[1])
  transp_ann1 = do.call(rbind, transp_ann1)
  transp_ann1 = cbind(filename=gsub("\\.[0-9]+$", "", row.names(transp_ann1)), transp_ann1)
  colnames(transp_ann1)[colnames(transp_ann1) == "X.Query_id"] = "locus_tag"
  
  #! filtering
  e_val = 1e-6 # specify e-value threshold
  TMP_over = 1 # min number of trans-membrane alpha-helical counts
  transp_ann2 = transp_ann1[transp_ann1$TM_Overlap_Score >= TMP_over &
                              transp_ann1$TM_Overlap_Score != "None" &
                              transp_ann1$e.value <= e_val, ]
  
  #! parsing data
  transp_ann3 = data.frame(transp_ann2[, c("filename","locus_tag","Match_length","e.value",
                                           "X._identity","Query_Length","Hit_Length","Query_Coverage",
                                           "Hit_Coverage","Query_n_TMS","Hit_n_TMS","TM_Overlap_Score")],
                           accession=sapply(strsplit(transp_ann2$Hit_desc, " "), function(x) paste(x[c(1:2)], collapse = " ")), 
                           Hit_desc=sapply(strsplit(transp_ann2$Hit_desc, " "), function(x) paste(x[-c(1:2)], collapse = " ")),
                           Predicted_Substrate=gsub("CHEBI:[0-9]+;", "", transp_ann2[, "Predicted_Substrate"]),
                           CHEBI=sapply(str_extract_all(transp_ann2[, "Predicted_Substrate"], "CHEBI:[0-9]+"), paste0, collapse=";"),
                           transp_ann2[, c("Family_Abrv","Hit_tcid")])
  colnames(transp_ann3)[colnames(transp_ann3) == "Hit_tcid"] = "gene_ann"
  transp_ann3$trait = paste0(transp_ann3$gene_ann, "|[", gsub(" ", "", transp_ann3$Predicted_Substrate), "]")
  
  
  ##-- select a specific transporter set
  
  if (transp_set == "ALL") {
    transp_ann4 = transp_ann3
    
  } else if (transp_set == "B1_B12_B3_B7_SIDERO") {
    chebi_sid = get_chebiid("siderophore")
    transp_ann4 = transp_ann3[sapply(transp_ann3$CHEBI, function(i) {
      any(strsplit(i, ";")[[1]] %in%
            c("CHEBI:15956", "CHEBI:41236", "CHEBI:13905", "CHEBI:22882", "CHEBI:22884", "CHEBI:3108" #B7
              , "CHEBI:30411" #B12
              , "CHEBI:17439", "CHEBI:60496", "CHEBI:48820", "CHEBI:3979", "CHEBI:14041", "CHEBI:23435" #CYANO-B12
              , "CHEBI:17154", "CHEBI:44258", "CHEBI:7556", "CHEBI:14645", "CHEBI:25521" #B3-nicotinamide
              , "CHEBI:15940", "CHEBI:44319", "CHEBI:7559", "CHEBI:25538" #B3-nicotinic acid
              , "CHEBI:9532" #B1-2P
              , "CHEBI:9533" #B1-1P
              , "CHEBI:18385", "CHEBI:46393", "CHEBI:9530", "CHEBI:15227", "CHEBI:26941" #B1-noP
              , chebi_sid$chebiid #siderophores
            )
      )}), ]
    
    # } else if (transp_set == "VITAMINS_SIDERO") {
    #   chebi_sid = get_chebiid("siderophore")
    #   chebi_vit = get_chebiid("vitamin", from = "all", max_res = 500) #??? not all vitamins are called vitamin!
    #   transp_ann4 = transp_ann3[sapply(transp_ann3$CHEBI, function(i) {
    #     any(strsplit(i, ";")[[1]] %in%
    #           c(chebi_sid$chebiid, chebi_vit$chebiid
    #           )
    #     )}), ]
    
  } else if (transp_set == "WITH_SUBSTRATE") {
    transp_ann4 = transp_ann3[transp_ann3$Predicted_Substrate != "None", ]
    
  } else if (transp_set == "NO_GENERIC") {
    bad_cpd = c("None","hydron","ion","cation","anion","proton","sodium(1+)",
                "potassium(1+)","calcium(2+)","lithium(1+)","chloride","electron",
                "fluoride","barium(2+)","strontium(2+)","cadmium(2+)","lead(2+)",
                "cobalt(2+)","mercury(2+)","nickel(2+)","zinc(2+)","dioxygen",
                "carbon dioxide","rubidium(1+)","iron(3+)","iron(2+)","manganese(2+)",
                "magnesium(2+)","metal cation","molecule","copper(1+)","zinc ion",
                "water","metabolite","dicarboxylic acid dianion","copper cation",
                "surfactant","vanadium oxoanion","monocarboxylic acid anion",
                "calcium ion","inorganic anion","base","inorganic cation",
                "monoatomic monocation","sodium atom","lithium atom","copper(2+)",
                "chromate(2-)","silver(1+)","selenite(2-)","detergent",
                "hydrogencarbonate","arsenite(3-)","arsenate(3-)",
                "arsenic molecular entity","benzalkonium chloride",
                "sodium dodecyl sulfate","CCCP","antimonite","aluminium(3+)",
                "sodium tungstate","molybdate","bromide","sodium iodide",
                "chlorate","bromate","periodate","thiocyanate","tetrafluoroborate(1-)",
                "nitric acid","silicic acid","sodium chloride","potassium chloride",
                "tungstate","silicate(4-)","silicon atom","caesium(1+)","mercury(0)")
    transp_ann4 = transp_ann3[!sapply(transp_ann3$Predicted_Substrate, function(z) all(strsplit(z, ", ")[[1]] %in% bad_cpd)), ]
    
  } else { stop("You have chose an invalid transporter set filter.\n
                Options are: ALL - B1_B12_B3_B7_SIDERO - WITH_SUBSTRATE - NO_GENERIC") }
  ##---
  
  
  #! save output
  write.table(transp_ann4, file = BioV_transp_tab, 
              col.names = T, row.names = F, quote = F, sep = "\t")
}


parse_gblast(transp_folder = snakemake@input[["transp_folder"]],
             transp_set = snakemake@params[["transp_set"]],
             BioV_transp_tab = snakemake@output[["BioV_transp_tab"]]
)