cluster_hmm <- function (to_cluster){

 cluster_num <- array (NA,length(to_cluster))
 cluster_num [1] <- 1

 for (pp in 2:length (to_cluster)){

  if (to_cluster[pp] != to_cluster[(pp-1)]){
   cluster_num[pp] <- cluster_num[(pp-1)] + 1
  } else {
   cluster_num[pp] <- cluster_num[(pp-1)]
  }
 }
 cluster_num
}


library(tidyverse)
library(GlmSimulatoR)
library(ggplot2)
library(stats)
library(depmixS4)
library(fitdistrplus)
library(spdep)
library(ape)
library(actuar)
library(data.table)
library(gridExtra)
library(ggplot2)
library(dplyr)


### READ INPUT DATA --  GEA RESULTS
my_results <- readRDS("likelihood_ratio_results_8.Old_Baldy.rds")
my_results$p_value <- as.numeric(as.character(my_results$p_value))
my_results$SNP <- gsub("_.*", "", my_results$SNP)

# Format Pg- to Pg.
lg <- read.table("final_order.txt")

# Perform left join with scaffolds map to add LG group to SNPs
merged_data <- left_join(my_results, lg, by = c("SNP" = "V4"))
colnames(merged_data) <- c("SNP", "p_value","LG","cM","cDNA")
pval <- merged_data$p_value
pval <- as.data.frame(pval)
colnames(pval) <- "V1"
ld_chroms <- merged_data$LG
ld_pos <- as.numeric(row.names(merged_data))
which_chroms <- unique(ld_chroms)

### Window size setting
window_size <- 1

# Number of states tested
max_states <- 15

# Results set up
gene_results <- array (NA, c((length(window_size)*length(which_chroms)*1),51))
the_crits <- array (NA,c((max_states),2))

pval_log_List <- rep(list(vector( "numeric" ,69)),length(which_chroms))
pos_List <- rep(list(vector( "numeric" ,69)),length(which_chroms))

window_pval_log_List <- rep(list(vector( "numeric" ,69)),length(which_chroms))
window_pos_List <- rep(list(vector( "numeric" ,69)),length(which_chroms))

pval_log_List2 <- rep(list(vector( "numeric" ,69)),length(which_chroms))
pos_List2 <- rep(list(vector( "numeric" ,69)),length(which_chroms))


# Removing possible NAs
for (j in 1:length (which_chroms)){

 pval_log_List[[j]] <- -1* log10( pval[ld_chroms == which_chroms[j],1])  ##convert p-values to -log10
 pos_List[[j]] <- ld_pos[ld_chroms == which_chroms[j]]

 pval_log_List2[[j]]<-pval_log_List[[j]][which(!is.na(pval_log_List[[j]]))]
 pos_List2[[j]]<-pos_List[[j]][which(!is.na(pval_log_List[[j]]))]

}


# Windows parameters set up
chromo_param <- array (NA,c((length (which_chroms)),1))
window_pos_List2<-rep(list(vector( "numeric" ,69)),length(which_chroms))


for (j in 1:length (which_chroms)){

 chromo_param[j,1]<-floor(length(pos_List2[[j]])/window_size)
}


# Save windows boundaries
window_pos_differences_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))
window_pos_low_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))
window_pos_high_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))



for (j in 1:length (which_chroms)){

 window_pos_differences_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,1]-1),byrow = TRUE)
 window_pos_low_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,1]-1),byrow = TRUE)
 window_pos_high_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,1]-1),byrow = TRUE)

 for(k in 1:((chromo_param[j,1]-1))){
  difference_pos<-(pos_List2[[j]][(k-1)*window_size+window_size])-(pos_List2[[j]][(k-1)*window_size+1])
  window_pos_differences_window3_List[[j]][1,k]<-difference_pos
  window_pos_low<-(pos_List2[[j]][(k-1)*window_size+1])
  window_pos_low_window3_List[[j]][1,k]<-window_pos_low
  window_pos_high<-(pos_List2[[j]][(k-1)*window_size+window_size])
  window_pos_high_window3_List[[j]][1,k]<-window_pos_high
 }

}




window_maxlogp_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))
window_minlogp_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))
window_meanlogp_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))
window_meanpos_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))

## Save windows p-values info
for (j in 1:length (which_chroms)){

 window_maxlogp_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,1]-1),byrow = TRUE)
 window_minlogp_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,1]-1),byrow = TRUE)
 window_meanlogp_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,1]-1),byrow = TRUE)
 window_meanpos_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,1]-1),byrow = TRUE)

 the_pos <- pos_List2[[j]]

 this_chrom_logp<-pval_log_List2[[j]][1:(window_size * (chromo_param[j,1]-1))]
 le_window<-window_size
 numb_le_window<-(chromo_param[j,1]-1)
 le_window_factor <- sort (array(1:numb_le_window,(numb_le_window * le_window)))
 minlogp<-tapply (this_chrom_logp[1:(le_window *numb_le_window)],le_window_factor,min,na.rm=TRUE)
 maxlogp<-tapply (this_chrom_logp[1:(le_window *numb_le_window)],le_window_factor,max,na.rm=TRUE)
 meanlogp<-tapply (this_chrom_logp[1:(le_window *numb_le_window)],le_window_factor,mean,na.rm=TRUE)
 the_meanpos<-tapply (the_pos[1:(le_window *numb_le_window)],le_window_factor,mean)
 window_minlogp_window3_List[[j]]<-minlogp
 window_maxlogp_window3_List[[j]]<-maxlogp
 window_meanlogp_window3_List[[j]]<-meanlogp
 window_meanpos_window3_List[[j]]<-the_meanpos


}


#reverse windows
window_meanpos_window3_List_rev<-window_meanpos_window3_List
window_meanlogp_window3_List_rev<-window_meanlogp_window3_List

for (j in 1:length (which_chroms)){	 #work one chromosome at a time

 window_meanpos_List_rev_Temp<-rev(window_meanpos_window3_List_rev[[j]])
 window_meanpval_log_rev_Temp<-rev(window_meanlogp_window3_List_rev[[j]])

 window_meanpos_window3_List_rev[[j]]<-window_meanpos_List_rev_Temp
 window_meanlogp_window3_List_rev[[j]]<-window_meanpval_log_rev_Temp

}



# HMM starts here - 1 seed

for(sd in 1:10){


 Mix_mod_List<-vector("list", (length(max_states)))
 Fit_Mix_mod_List<-vector("list", (length(max_states)))

 for (pp in 1:max_states){
  set.seed((1000+sd))

  Mix_mod_List_sublist<-vector("list", (length(which_chroms)))
  Fit_Mix_mod_List_sublist<-vector("list", (length(which_chroms)))

  for(i in 1:length(which_chroms)){
   Mix_mod_List_sublist[[i]]<-depmix(window_meanlogp_window3_List[[i]] ~ 1,nstates = pp, family = gaussian(link = "inverse"), ntimes=length((window_meanlogp_window3_List[[i]])))
   Fit_Mix_mod_List_sublist[[i]]<-fit(Mix_mod_List_sublist[[i]],em=em.control(maxit=2500000))
  }
  Mix_mod_List[[pp]]<-Mix_mod_List_sublist
  Fit_Mix_mod_List[[pp]]<-Fit_Mix_mod_List_sublist
 }


 the_crits_list<-vector("list", (length(max_states)))

 for(i in 1:length(which_chroms)){
  the_crits <- array (NA,c((max_states),2))
  for (pp in 1:(max_states)){
   the_crits[pp,1] <- AIC((Fit_Mix_mod_List[[pp]][[i]]))
   the_crits[pp,2] <- BIC((Fit_Mix_mod_List[[pp]][[i]]))
  }
  the_crits_list[[i]]<-the_crits
 }


 optimal_statenums<-rep(0, (length(which_chroms)))

 for(i in 1:length(which_chroms)){

  mod_diff <- the_crits_list[[i]][1:(max_states - 1),2] - the_crits_list[[i]][2:(max_states),2]
  diff_neg <- mod_diff < 0

  if (sum (diff_neg) == 0){
   best_mod <- max_states
  } else {
   best_mod <- min (which (diff_neg)) #find the first value that is negative
  }
  optimal_statenums[i]<-best_mod
 }

 optimal_statenums

 the_crits_list[[1]]


 testing_statenums2<-c(2,max_states)


 optimal_Mix_mod_List<-vector("list", (length(which_chroms)))
 optimal_Fit_Mix_mod_List<-vector("list", (length(which_chroms)))

 for(j in 1:length(which_chroms)){
  optimal_Mix_mod_List[[j]]<-Mix_mod_List[[(optimal_statenums[j])]][[j]]
  optimal_Fit_Mix_mod_List[[j]]<-Fit_Mix_mod_List[[(optimal_statenums[j])]][[j]]
 }



 est.states_variablestatenum<-vector("list", (testing_statenums2[2]-testing_statenums2[1]+1))

 for (k in testing_statenums2[1]:testing_statenums2[2]){

  est.states_variablestatenum_sublist<-vector("list", (length(which_chroms)))

  for(i in 1:length(which_chroms)){
   est.states_variablestatenum_sublist[[i]]<-(posterior(Fit_Mix_mod_List[[k]][[i]]))$state
  }
  est.states_variablestatenum[[(k-testing_statenums2[1]+1)]]<-est.states_variablestatenum_sublist

 }


 Fitted_Item_reorderstates_variablestatenum<-vector("list", (testing_statenums2[2]-testing_statenums2[1]+1))
 Fitted_Item_reorderstates_index_variablestatenum<-vector("list", (testing_statenums2[2]-testing_statenums2[1]+1))
 Fitted_Item_reordertrans_variablestatenum<-vector("list", (testing_statenums2[2]-testing_statenums2[1]+1))

 for (k in testing_statenums2[1]:testing_statenums2[2]){

  Fitted_Item_reorderstates<-vector("list", (length (which_chroms)))
  Fitted_Item_reorderstates_index<-vector("list", (length (which_chroms)))
  Fitted_Item_reordertrans<-vector("list", (length (which_chroms)))


  for(j in 1:length (which_chroms)){

   h_states<-getpars(Fit_Mix_mod_List[[k]][[j]])[seq((k+k*k+1),(k+k*k+2*k),2)]
   names(h_states)<-seq(1,k,1)
   h_states2<-h_states[order(h_states,decreasing=TRUE)]
   h_states2_index<-order(h_states,decreasing=TRUE)

   Fitted_Item_reorderstates[[j]]<-h_states2
   Fitted_Item_reorderstates_index[[j]]<-h_states2_index

   h_states_trans<-t(array((getpars(Fit_Mix_mod_List[[k]][[j]])[seq((k+1),(k+k*k),1)]),c(k,k)))
   index_statenames<-seq(1,k,1)
   rownames(h_states_trans) <- index_statenames
   colnames(h_states_trans) <- index_statenames
   h_states_trans_reorder<-h_states_trans[h_states2_index,h_states2_index]

   Fitted_Item_reordertrans[[j]]<-h_states_trans_reorder

  }

  Fitted_Item_reorderstates_variablestatenum[[(k-testing_statenums2[1]+1)]]<-Fitted_Item_reorderstates
  Fitted_Item_reorderstates_index_variablestatenum[[(k-testing_statenums2[1]+1)]]<-Fitted_Item_reorderstates_index
  Fitted_Item_reordertrans_variablestatenum[[(k-testing_statenums2[1]+1)]]<-Fitted_Item_reordertrans

 }



 est.states_variablestatenum_reorder<-vector("list", (testing_statenums2[2]-testing_statenums2[1]+1))

 for(k in testing_statenums2[1]:testing_statenums2[2]){

  est.states_variablestatenum_reorder_sub<-vector("list", (length (which_chroms)))

  for(j in 1:length(which_chroms)){
   est.states_variablestatenum_reorder_sub[[j]]<-est.states_variablestatenum[[(k-testing_statenums2[1]+1)]][[j]]
   for(i in 1:length ((Fitted_Item_reorderstates_index_variablestatenum[[(k-testing_statenums2[1]+1)]][[j]]))){
    est.states_variablestatenum_reorder_sub[[j]][(est.states_variablestatenum_reorder_sub[[j]])==(Fitted_Item_reorderstates_index_variablestatenum[[(k-testing_statenums2[1]+1)]][[j]][i])] <- (100+i)
   }
   est.states_variablestatenum_reorder_sub[[j]]<-est.states_variablestatenum_reorder_sub[[j]]-100
  }
  est.states_variablestatenum_reorder[[(k-testing_statenums2[1]+1)]]<-est.states_variablestatenum_reorder_sub

 }



 cluster_est.states_variablestatenum_reorder<-vector("list", (length (est.states_variablestatenum_reorder)))

 for(k in 1:length(est.states_variablestatenum_reorder)){
  cluster_est.states_variablestatenum_reorder_sub<-vector("list", (length (which_chroms)))

  for(j in 1:length (which_chroms)){
   cluster_est.states_variablestatenum_reorder_sub[[j]] <- cluster_hmm(est.states_variablestatenum_reorder[[k]][[j]])
  }
  cluster_est.states_variablestatenum_reorder[[k]]<-cluster_est.states_variablestatenum_reorder_sub
 }





 # window stats tables

 Allwindow_stat_List<-vector("list", (length (which_chroms)))


 for(j in 1:length (which_chroms)){
  Allwindow_stat_List[[j]]<- array (0,c((length(window_meanlogp_window3_List[[j]])),(8+25*(length(Fitted_Item_reorderstates_variablestatenum)))))

  windowpos<-seq(1,(length(window_meanlogp_window3_List[[j]])),1)
  Allwindow_stat_List[[j]][,1]<-windowpos

  for(i in 1:(length(window_meanlogp_window3_List[[j]]))){

   Allwindow_stat_List[[j]][i,2]<-pos_List2[[j]][(Allwindow_stat_List[[j]][i,1]-1)*window_size+1]
   Allwindow_stat_List[[j]][i,3]<-pos_List2[[j]][(Allwindow_stat_List[[j]][i,1])*window_size]
   Allwindow_stat_List[[j]][i,4]<-(pos_List2[[j]][(Allwindow_stat_List[[j]][i,1]-1)*window_size+1]+pos_List2[[j]][(Allwindow_stat_List[[j]][i,1])*window_size])/2


   for(k in 1:length(Fitted_Item_reorderstates_variablestatenum)){
    Allwindow_stat_List[[j]][i,7+4*(k-1)+1]<-est.states_variablestatenum_reorder[[k]][[j]][i]
    Allwindow_stat_List[[j]][i,7+4*(k-1)+2]<-cluster_est.states_variablestatenum_reorder[[k]][[j]][i]


   }
  }

  window_mean<-window_meanlogp_window3_List[[j]]
  Allwindow_stat_List[[j]][,5]<-window_mean

  window_max<-window_maxlogp_window3_List[[j]]
  Allwindow_stat_List[[j]][,6]<-window_max

  window_min<-window_minlogp_window3_List[[j]]
  Allwindow_stat_List[[j]][,7]<-window_min


 }



 # Get min of peak state and max of all non-peak states for each chromosomes
 minpeakANDmaxnonpeak_thres<-array (0,c((length (which_chroms)),4*length(Fitted_Item_reorderstates_variablestatenum)))

 for(j in 1:length(which_chroms)){

  for(k in 1:length(Fitted_Item_reorderstates_variablestatenum)){
   Temp_states<-sort(unique(Allwindow_stat_List[[j]][,(7+4*(k-1)+1)]),decreasing=FALSE)
   Temp_minmean_peakwindow<-min(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,(7+4*(k-1)+1)] %in% c(1)),5])
   Temp_maxmean_nonpeakwindow<-max(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,(7+4*(k-1)+1)] %in% (Temp_states[2:length(Temp_states)])),5])
   minpeakANDmaxnonpeak_thres[j,(4*(k-1)+1)]<-Temp_minmean_peakwindow
   minpeakANDmaxnonpeak_thres[j,(4*(k-1)+2)]<-Temp_maxmean_nonpeakwindow

  }
 }




 #Let's fill in I_window, I_window_rev, I_window3, I_window3_rev variables

 for(j in 1:length (which_chroms)){
  for(k in 1:length(Fitted_Item_reorderstates_variablestatenum)){
   # I_window
   for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
    Allwindow_stat_List[[j]][i,(7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<-0
    if(Allwindow_stat_List[[j]][i,(7+4*(k-1)+1)]==1){
     Allwindow_stat_List[[j]][i,(7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<-1
    }
   }



  }


 }

 # Save results


 saveRDS(Allwindow_stat_List, paste(c("8.Old_Baldy_ReOrder_window_",window_size,"_window_stats_seed_",(1000+sd),".rds"),collapse=""))

 saveRDS(testing_statenums2, paste(c("8.Old_Baldy_ReOrder_window_",window_size,"_state_number_seed_",(1000+sd),".rds"),collapse=""))

 saveRDS(the_crits_list, paste(c("8.Old_Baldy_ReOrder_window_",window_size,"_the_crits_list__seed_",(1000+sd),".rds"),collapse=""))

}
