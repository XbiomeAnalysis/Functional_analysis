run_permanova_betadisp <- function(physeq, permanova =T, betadisp = T, vars){ 
  
  # distance calculation and dispersion comparison 
  
  otu <- microbiome::abundances(physeq) 
  
  meta <- microbiome::meta(physeq) 
  
  dis <- vegdist(t(otu)) 
  
  group <- meta[[vars]] 
  
  if(betadisp){ 
    
    mod <- betadisper(dis, group) 
    
    dis_p <- permutest(mod, pairwise = TRUE, permutations = 999) 
    
    betadisp_res <- data.frame("variable" = vars, "p_value"= dis_p$tab$`Pr(>F)`[1], "analysis" = "beta_dispersion_permutation999") 
    
  }else{betadisp_res <- data.frame()} 
  
  if(permanova){ 
    
    # adonis test 
    
    colnames(meta)[colnames(meta)==vars] <- "vars_col" 
    
    permanova <- adonis(t(otu) ~ vars_col, data = meta, permutations=999, method = "bray") 
    
    # results 
    
    permanova_res <- data.frame("variable" = vars, "p_value"= permanova$aov.tab$`Pr(>F)`[1], "R2"=permanova$aov.tab$R2[1], "analysis" = "permanova_permutation999") 
    
  }else{ 
    
    permanova_res <- data.frame() 
    
  } 
  
  return(list(betadisp_res=betadisp_res, permanova_res=permanova_res)) 
  
} 