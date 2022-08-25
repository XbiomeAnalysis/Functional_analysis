# this script is used to get associated gene family for specific pathways 
# http://172.16.100.150:8080/browse/GCIP-238
# input file includes 
# humann2-2.8.2/humann2/data/pathways/metacyc_reactions_level4ec_only.uniref.bz2 - this is reaction to gene family file 
# humann2-2.8.2/humann2/data/pathways/metacyc_pathways - this is pathway to reaction file 

pathway2gene <- function(target_pathway){
  
  library(data.table)
  library(dplyr)
  library(reshape2)
  
  #target_pathway <- "BRANCHED-CHAIN-AA-SYN-PWY"

  # pathway to rxn 
  tmp1 <- read.csv("/home/xuxiaomin/project/standardized_analytics_workflow_R_function/demo_data/Funtional/metacyc_data/metacyc_pathways", header = F) %>% data.table(.) %>% .[grepl(target_pathway, V1)] 
  
  if(nrow(tmp1)==0){
    return(NULL)
  }
  
  tmp1$V1 <- gsub("\t"," ", tmp1$V1)
  tmp1 <- as.data.frame(stringi::stri_split_fixed(tmp1, " ", simplify = TRUE)) %>% melt(.,id.var = "V1", value.name = "reaction")
  tmp1$variable <- NULL
  
  # rxn to gene family 
  tmp2 <- read.csv("/home/xuxiaomin/project/standardized_analytics_workflow_R_function/demo_data/Funtional/metacyc_data/metacyc_reactions_level4ec_only.uniref", header = F) %>% data.table(.)
  
  tmp4 <- data.table()
  for(i in tmp1$reaction){
    tmp3 <- tmp2 %>% .[grepl(i, V1)] 
    tmp4 <- rbindlist(list(tmp3,tmp4))
  }
  
  tmp4$V1 <- gsub("\t"," ", tmp4$V1)
  
  test <- tmp4[,c("reaction","ec","genefamily"):=data.table(str_split_fixed(V1," ",3))]
  
  splits <- max(lengths(strsplit(tmp4$genefamily, " ")))
  test2 <- tmp4[, c("reaction","genefamily")] %>% .[,paste0("V", 1:splits) := tstrsplit(genefamily, " ", fixed=T)] %>% melt(.,id.var = "reaction", value.name = "genefamily")
  test2$variable <- NULL
  test2 <- test2[!grepl("UniRef50", genefamily)] %>% .[!is.na(genefamily) | !genefamily == "0"] %>% unique(.)
  
  return(test2 %>% as.data.frame())
  
  #fwrite(test2, file = paste0(output_pwd, paste0(target_pathway,"_associated_reaction_genefamily.csv")), row.names = F)

}

