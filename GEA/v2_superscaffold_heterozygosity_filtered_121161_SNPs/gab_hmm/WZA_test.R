setwd("/Users/gnocc/Desktop/Projects/Banff/GEA/v2_superscaffold_heterozygosity_filtered_121161_SNPs/")
norquay <- readRDS("likelihood_ratio_results_4.Norquay.rds")


 assign.pvalues <- function(array){
    #array <- sample(sw, 1000)
    pvalues <- array(0, length(array))

    ordered.indexes <- order(array, decreasing = TRUE)

    j <- length(array)
    for( i in ordered.indexes ){
      pvalues[i] <- j/length(array)
      j <- j-1
    }

    return(pvalues)
  }
  
  
 norquay$emp <- assign.pvalues(as.numeric(as.character(norquay$p_value)))