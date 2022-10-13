KMdiagram_fetcher = function(ncore=1, path) {
  #    "ncore": degree of parallelization, can be > of real number of threads but don't exagerate otherwise the fetching script will rise an error
  #    "path": path where to save (as .rds) the fetched KM structure object
  
  
  ## load packages
  library(stringr)
  library(xml2)
  library(httr)
  library(KEGGREST)
  library(future.apply)
  
  
  ## get list of KM names
  KM_list <- keggList("module")
  KM_list <- str_extract(names(KM_list), "M[0-9]{5}")
  
  
  ## fetch diagrams
  plan("multicore", workers = ncore)
  fetched_KM <- future_lapply(KM_list, function(imod) { 
    txt <- url(paste0("http://www.genome.jp/kegg-bin/show_module?", imod))
    #??? issue with next command
    # txt1 <- with_config(user_agent(txt$request$options$useragent), 
    #                    read_html(txt, encoding = "UTF-8", verbose=T))
    txt1 <- read_html(txt, encoding = "UTF-8", verbose=T)
    txt2 <- gsub(".*<map id=\"module\" name=\"module\">|</map>.*", "", txt1)
    txt3 <- strsplit(txt2, "<area shape=\"rect\" ")[[1]][-1]
    txt4 <- t(as.data.frame(strsplit(txt3, '\\" ')))[, c(1,3,5), drop=F]
    txt5 <- data.frame(row.names = NULL, KM_lvl=gsub("_.*", "", gsub('id=\"', "", txt4[,1])), KO_id=gsub(".*_", "", gsub('id=\"', "", txt4[,1])),
                      KO=gsub('title="| .*', "", txt4[,2]), t(as.data.frame(strsplit(gsub('coords="|">.*', "", txt4[,3]), ","))))
    names(txt5)[4:7] <- c("X1","Y1","X2","Y2")
    #! remove optional terms '-K00000'
    txt_str <- keggGet(paste0("md:", imod))[[1]]$DEFINITION
    optKOs <- unlist(c(gsub("\\-", "", str_extract_all(txt_str, "\\-K[0-9]{5}")[[1]]),
                      str_extract_all(str_extract_all(txt_str, "\\-\\([K0-9,]*[^)]\\)")[[1]], "K[0-9]{5}")))
    txt_str_clean <- gsub(paste0("-", optKOs, collapse = "|"), "", txt_str)
    txt6 <- txt5[!txt5$KO %in% optKOs, ]
    #---
    
    KM_len <- length(unique(txt6$KM_lvl))
    
    #! if single-level KM check whether it just contain a single macromolecular complex (simplify completness/specificity analysis)
    if (KM_len == 1 & str_count(paste(txt_str, collapse = " "), "\\(|\\)|,| ") == 0) { just_macro = TRUE
    } else { just_macro = FALSE }
    #---
    
    KM_len_macro <- ""; if (just_macro == TRUE) { KM_len_macro = str_count(txt_str_clean, "K[0-9]{5}") }
    iKM <- list(just_macro=just_macro, length=KM_len, length_macro=KM_len_macro, str=txt6)
    
    return(list(optKOs, iKM))
  } )
  
  KM_str <- lapply(fetched_KM, function(i) i[[2]])
  names(KM_str) <- KM_list
  KM_str <- KM_str[order(as.integer(gsub("M", "", KM_list)))]
  opt_KOs <- lapply(fetched_KM, function(i) i[[1]]) # list of otional KOs in modules
  names(opt_KOs) <- KM_list
  opt_KOs <- opt_KOs[order(as.integer(gsub("M", "", KM_list)))]
  
  
  ## add KM description
  KM_htext <- readLines("https://www.genome.jp/kegg-bin/download_htext?htext=ko00002.keg&format=htext&filedir=")
  KM_htext_tab <- data.frame(matrix(NA, nrow = length(grep("^D", KM_htext)), ncol = 4))
  count <- 0
  for (i in 1:length(KM_htext)) {
    if (length(grep("^A", KM_htext[i])) > 0) {
      i_A <- sapply(strsplit(KM_htext[i], "<|>"), "[[", 3)
    }
    if (length(grep("^B ", KM_htext[i])) > 0) {
      i_B <- sapply(strsplit(KM_htext[i], "<|>"), "[[", 3)
    }
    if (length(grep("^C", KM_htext[i])) > 0) {
      i_C <- sapply(strsplit(KM_htext[i], "    "), "[[", 2)
    }
    if (length(grep("^D", KM_htext[i])) > 0) {
      i_D <- sapply(strsplit(KM_htext[i], "      "), "[[", 2)
      count <- count+1
      KM_htext_tab[count, ] <- c(i_A,i_B,i_C,i_D)
    }
  }
  KM_names <- names(KM_str)
  KM_str <- lapply(seq_along(KM_str), function(j) {
    j2 <- as.list(KM_htext_tab[grepl(names(KM_str)[j], KM_htext_tab[, 4]), ])
    j2 <- lapply(j2, function(z) gsub("[ ]+", "_", z))
    names(j2) <- paste0("KM_level", 1:4)
    j2$KM_desc <- gsub("_\\[.*", "", j2$KM_level4)
    j3 <- c(KM_str[[j]], j2)
    return(j3)
  })
  names(KM_str) <- KM_names
  
  
  ## save output
  dir.create(path, recursive = T, showWarnings = F)
  saveRDS(KM_str, file = file.path(path, "KM_str.rds"))
  # saveRDS(opt_KOs, file = file.path(path, "opt_KOs.rds"))
}

KMreco = function(indata, KM_str, len_breaks=NULL, allowed_gaps=c(0), ncore=1) {
  #!!!    "indata": file must be a presence/absence (i.e. 1/0) matrix/data.frame 
  #       with genomes as rows and annotated KOs as columns
  #!!!    "KM_str": contains KM diagrams, obtained with 'KMdiagram_fetcher' function
  #!!!    "len_breaks": vector of integer of open-right limits for binning of KM length; e.g. c(3, 10)
  #!!!    "allowed_gaps": vector of integer specfing nomber of gaps allowed in each bin; it has length "breaks"+1; e.g. c(0,1,2) 
  #!!!    "ncore": degree of parallelization
  
  
  #! libraries
  library(future.apply)
  library(tidyverse)
  
  
  #! set variables
  indata = as.data.frame(indata)
  options(stringsAsFactors = F)
  
  
  
  ##### Parse genomes to get KM completness
  plan("multicore", workers = ncore)
  KMgnm_compl = future_lapply(row.names(indata), function(gnm) { # loop thru genomes
    KM_compl = data.frame(matrix(ncol = 3, nrow = 0))
    names(KM_compl) = c("trait", "gene_count", "gene_ann")
    
    # get KO list
    i_gnm = indata[gnm, , drop=F]
    i_gnmKO = colnames(i_gnm)[i_gnm > 0]
    
    # loop through all KM
    for (KM in names(KM_str)) { 
      i_KM_str = KM_str[[KM]]
      
      #1 look which KOs are present in the genome
      i_KM_str$str$exist = i_KM_str$str$KO %in% i_gnmKO
      if (sum(i_KM_str$str$exist) == 0) { #!!! when no KO of the KM is present in the genome
        i_KM_compl = data.frame(trait = KM,
                                gene_count = 0,
                                gene_ann = NA)
        KM_compl = rbind.data.frame(KM_compl, i_KM_compl)
      } else { #!!! when at least some KOs of the KM are present in the genome
        
        #2 per each KM level get completness
        if (i_KM_str$just_macro == T) { # when single-level macromolecular complex
          i_KM_existKO = i_KM_str$str$exist
          
        } else { # when is not a single-level macromolecular complex
          i_KM_existKO = c()
          for (i_lvl in unique(i_KM_str$str$KM_lvl)) { # loop through levels of that KM
            i_lvl_KM_str = i_KM_str$str[i_KM_str$str$KM_lvl %in% i_lvl, , drop=F]
            i_lvl_KM_str[, c(4:7)] = apply(i_lvl_KM_str[, c(4:7)], 2, function(x) as.numeric(x)) # fix
            
            if (nrow(i_lvl_KM_str) == 1) { #!!! if level has single-KO
              i_KM_existKO = c(i_KM_existKO, i_lvl_KM_str$exist)
              
            } else { #!!! if level has multiple-KOs
              #2a     find first row of KOs
              sol_starts = which(i_lvl_KM_str$Y1 == min(i_lvl_KM_str$Y1))
              #2b    looks for presence of gaps between gene1 X2 - gene2 X1
              sol_x2_lim = c()
              if (length(sol_starts) > 1) { # multiple hypotetic solution
                for (l in 1:(length(sol_starts)-1)) {
                  if (i_lvl_KM_str[sol_starts[l], "X2"] != i_lvl_KM_str[sol_starts[l+1], "X1"] & # when there is a gap
                      i_lvl_KM_str[sol_starts[l+1], "X1"] - i_lvl_KM_str[sol_starts[l], "X2"] == 4) { #!! fix fake gap that can be with nested redundant level of different length (e.g. M00003)
                    sol_x2_lim = c(sol_x2_lim, i_lvl_KM_str[sol_starts[l], "X2"])
                    if (l == length(sol_starts)-1) { sol_x2_lim = c(sol_x2_lim, i_lvl_KM_str[sol_starts[l+1], "X2"]) } # add last level solution
                  } 
                }
              }
              if (length(sol_x2_lim) == 0) { sol_x2_lim = c(sol_x2_lim, i_lvl_KM_str[sol_starts[length(sol_starts)], "X2"]) } # when there is only one solution or none of hypotetic split was good 
              #2c    split level in all possible solutions - X-wise
              i_lvl_KM_str2 = list()
              for (l in 1:length(sol_x2_lim)) {
                i_lvl_KM_str2[[l]] = i_lvl_KM_str[i_lvl_KM_str$X2 <= sol_x2_lim[l], ]
                i_lvl_KM_str = i_lvl_KM_str[!i_lvl_KM_str$X2 <= sol_x2_lim[l], ] # remove parsed element from level
              }
              #2d    split each solution Y-wise
              for (l in 1:length(i_lvl_KM_str2)) {
                i_lvl_KM_str2[[l]] = split(i_lvl_KM_str2[[l]], f = as.factor(i_lvl_KM_str2[[l]]$Y1))
              }
              #2e    pull out all possible KO-combinations from each solution
              sol_paths_exist = list()
              for (l in 1:length(i_lvl_KM_str2)) { # loop thru solutions - for completness
                sol_paths_exist[[l]] = expand.grid(lapply(i_lvl_KM_str2[[l]], function(z) z$exist)) # use z$KO for debug
              }
              #2f    discard each solution with missing KOs
              sol_paths_exist2 = list()
              for (l in 1:length(sol_paths_exist)) {
                tmp_exist = sol_paths_exist[[l]][apply(sol_paths_exist[[l]], 1, function(y) sum(y) == length(y)), ]
                if (length(tmp_exist) > 0) { # keep only sol without 0
                  sol_paths_exist2 = c(sol_paths_exist2, list(tmp_exist))
                }
              }
              #2g    look whether there is a solution for this KM-level
              if (length(unlist(sol_paths_exist2)) > 0) {
                sol_exist = TRUE
              } else { sol_exist = FALSE}
              
              i_KM_existKO = c(i_KM_existKO, sol_exist)
            }
          }
        }
        
        #3 estimate KM completeness
        i_KM_compl = data.frame(trait = KM,
                                gene_count = sum(i_KM_existKO),
                                gene_ann = paste0(i_KM_str$str$KO[i_KM_str$str$exist], collapse = ";"))
        KM_compl = rbind.data.frame(KM_compl, i_KM_compl)
      }
    }
    
    KM_compl = cbind(filename=gnm, KM_compl)
    
    return(KM_compl)
  })
  KMgnm_compl = do.call(rbind, KMgnm_compl)
  KMgnm_compl$gene_count = as.numeric(KMgnm_compl$gene_count)
  
  
  
  ### apply completeness rules
  
  #! get KM lengths
  KM_len = unlist(lapply(KM_str, "[[", "length"))
  
  
  #! parse macromolecular complex KMs with same completness rules as other KMs
  macro_complx = unlist(lapply(KM_str, "[[", "just_macro"))
  KM_len[macro_complx] = as.integer(unlist(lapply(KM_str, "[[", "length_macro")))[macro_complx]
  
  
  #! estimate gaps allowed for each KM
  if (!is.null(len_breaks)) {
    KM_len_bin = .bincode(KM_len, breaks = c(0, len_breaks, Inf), right = F, include.lowest = T)
  } else {
    KM_len_bin = .bincode(KM_len, breaks = c(0, Inf), right = F, include.lowest = T)
  }
  KM_allowed_gaps = allowed_gaps[KM_len_bin]
  
  
  #! sort vectors
  KMgnm_compl$allowed_gaps = KM_allowed_gaps[match(KMgnm_compl$trait, names(KM_len))]
  KMgnm_compl$KM_len = KM_len[match(KMgnm_compl$trait, names(KM_len))]
  
  
  #! check for completeness
  KMgnm_compl$complete = (KMgnm_compl$gene_count + KMgnm_compl$allowed_gaps) >= KMgnm_compl$KM_len
  KMgnm_compl$complete_perc = (KMgnm_compl$gene_count + KMgnm_compl$allowed_gaps)/KMgnm_compl$KM_len
  KMgnm_compl$complete_perc[KMgnm_compl$complete_perc > 1] = 1
  
  
  #! pasre output
  KMgnm_compl_filt = KMgnm_compl[KMgnm_compl$complete_perc > (KMgnm_compl$allowed_gaps/KMgnm_compl$KM_len), ]
  KMgnm_compl_filt = KMgnm_compl_filt %>% separate_rows(gene_ann, sep = ";")
  KMgnm_compl_filt = KMgnm_compl_filt[, c(colnames(KMgnm_compl_filt)[!colnames(KMgnm_compl_filt) %in% c("gene_ann","trait")], "gene_ann","trait")]
  return(KMgnm_compl_filt) 
}