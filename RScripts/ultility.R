## 1. Filter pathway data
aggregate_low_abundance <- function(input_data,threshold){
  
  cur_data_row = nrow(input_data)
  
  input_data[cur_data_row + 1, ] = 0
  rownames(input_data)[cur_data_row + 1] = "LowAbundanceSum"

  for (i in 1:cur_data_row) {
    for (j in 1:ncol(input_data)) {
      cur = input_data[i, j]
      if (cur < threshold) {
        rowname = rownames(input_data)[i]
        target_row = which(rownames(input_data) == "LowAbundanceSum")
        input_data[target_row, j] = input_data[target_row, j] + cur
        input_data[i, j] = 0
      }
    }
  }
  input_data = input_data %>% rownames_to_column("ID") %>%
    dplyr::filter(rowSums(.[, -1]) > 0) %>% column_to_rownames("ID")
  
  return(input_data)
}



## 2. Filter low-prevalence data
filter_prevalence <- function(otu_table, threshold=0.02,taxa_are_rows = TRUE){
  if(threshold>1 | threshold<0){
    stop("threshold [0-1]")
  }
  else{
    if(!taxa_are_rows){
      otu_table <- t(otu_table)
      if(ncol(otu_table) * threshold < 2){
        filtered_otu_table<-otu_table[(rowSums(!otu_table==0)>=2),]
      }
      else{
        filtered_otu_table<- otu_table[(rowSums(!otu_table==0)>=ceiling(ncol(otu_table)*threshold)),]
      }
    }
    else{
      if(ncol(otu_table) * threshold < 2){
        filtered_otu_table<-otu_table[(rowSums(!otu_table==0)>=2),]
      }
      else{
        filtered_otu_table<- otu_table[(rowSums(!otu_table==0)>=ceiling(ncol(otu_table)*threshold)),]
      }
    }
  }
  return (filtered_otu_table)
}


## 3. dispersion & permanova test
run_permanova_betadisp <- function(physeq, permanova =T, betadisp = T, vars){
  # distance calculation and dispersion comparison
  otu <- microbiome::abundances(physeq)
  meta <- microbiome::meta(physeq)
  dis <- vegan::vegdist(t(otu))
  group <- meta[[vars]]
  
  if(betadisp){
    mod <- vegan::betadisper(dis, group)
    dis_p <- vegan::permutest(mod, pairwise = TRUE, permutations = 999)
    betadisp_res <- data.frame("variable" = vars, "p_value"= dis_p$tab$`Pr(>F)`[1], "analysis" = "beta_dispersion_permutation999")
  }
  else{
    betadisp_res <- data.frame()
  }
  if(permanova){
    # adonis test
    colnames(meta)[colnames(meta)==vars] <- "vars_col"
    permanova <- vegan::adonis(t(otu) ~ vars_col, data = meta, permutations=999, method = "bray")
    # results
    permanova_res <- data.frame("variable" = vars, "p_value"= permanova$aov.tab$`Pr(>F)`[1], "R2"=permanova$aov.tab$R2[1], "analysis" = "permanova_permutation999")
  }
  else{
    permanova_res <- data.frame()
  }
  return(list(betadisp_res=betadisp_res, permanova_res=permanova_res))
}