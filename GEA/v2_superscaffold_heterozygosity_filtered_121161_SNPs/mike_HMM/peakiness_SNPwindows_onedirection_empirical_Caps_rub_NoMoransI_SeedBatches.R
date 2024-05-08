rm(list = ls())


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


.libPaths( "/data/michael_HMM/Rlibs" )

date()

library(tidyverse)
library(GlmSimulatoR)
library(ggplot2)
library(stats)
library(depmixS4)
library(fitdistrplus)
library(spdep)
library(ape)
library(actuar)




Enviro_Vars<-read.table ("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/Caps_rub_list.tsv")

#for(qqq in 1:(1)){
for(qqq in 1:(nrow(Enviro_Vars))){
  
  
  Enviro_Vars_of_Interest<-Enviro_Vars[qqq,1]
  
  
  setwd("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual")
  pval <- read.table ((paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_pval.tsv"),collapse="")))
  names1 <- read.table ((paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_gene_names2.tsv"),collapse="")))
  
  
  ld_chroms <- sapply (strsplit(as.character (names1[,1]),split = "__"),"[[",1)
  ld_pos <- as.numeric (names1[,2])
  
  window_sizes <- c(10,25,50,50,500)
  which_chroms <- unique (ld_chroms)
  
  gene_results <- array (NA, c((length(window_sizes)*length(which_chroms)*1),51))
  
  max_states <- 11
  the_crits <- array (NA,c((max_states),2))
  
  
  pval_log_List<-rep(list(vector( "numeric" ,69)),length(which_chroms))
  pos_List<-rep(list(vector( "numeric" ,69)),length(which_chroms))
  
  window_pval_log_List<-rep(list(vector( "numeric" ,69)),length(which_chroms))
  window_pos_List<-rep(list(vector( "numeric" ,69)),length(which_chroms))
  
  pval_log_List2<-rep(list(vector( "numeric" ,69)),length(which_chroms))
  pos_List2<-rep(list(vector( "numeric" ,69)),length(which_chroms))
  
  
  for (j in 1:length (which_chroms)){
    
    pval_log_List[[j]] <- -1* log10( pval[ld_chroms == which_chroms[j],1])  ##convert p-values to -log10
    pos_List[[j]] <- ld_pos[ld_chroms == which_chroms[j]]
    
    pval_log_List2[[j]]<-pval_log_List[[j]][which(!is.na(pval_log_List[[j]]))]
    pos_List2[[j]]<-pos_List[[j]][which(!is.na(pval_log_List[[j]]))]
    
  }
  
  
  chromo_param <- array (NA,c((length (which_chroms)),7))
  
  window_pos_List2<-rep(list(vector( "numeric" ,69)),length(which_chroms))
  
  
  for (j in 1:length (which_chroms)){
    chromo_param[j,1]<-max(pos_List2[[j]])
    chromo_param[j,2]<-length(pos_List2[[j]])
    chromo_param[j,3]<-floor(length(pos_List2[[j]])/window_sizes[1])
    chromo_param[j,4]<-floor(length(pos_List2[[j]])/window_sizes[2])
    chromo_param[j,5]<-floor(length(pos_List2[[j]])/window_sizes[3])
    chromo_param[j,6]<-floor(length(pos_List2[[j]])/window_sizes[4])
    chromo_param[j,7]<-floor(length(pos_List2[[j]])/window_sizes[5])
  }
  
  
  window_pos_differences_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))
  window_pos_low_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))
  window_pos_high_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))
  
  
  for (j in 1:length (which_chroms)){
    
    window_pos_differences_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,5]-1),byrow = TRUE)
    window_pos_low_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,5]-1),byrow = TRUE)
    window_pos_high_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,5]-1),byrow = TRUE)
    
    for(k in 1:((chromo_param[j,5]-1))){
      difference_pos<-(pos_List2[[j]][(k-1)*window_sizes[3]+window_sizes[3]])-(pos_List2[[j]][(k-1)*window_sizes[3]+1])
      window_pos_differences_window3_List[[j]][1,k]<-difference_pos
      window_pos_low<-(pos_List2[[j]][(k-1)*window_sizes[3]+1])
      window_pos_low_window3_List[[j]][1,k]<-window_pos_low
      window_pos_high<-(pos_List2[[j]][(k-1)*window_sizes[3]+window_sizes[3]])
      window_pos_high_window3_List[[j]][1,k]<-window_pos_high
    }
    
  }
  
  
  
  window_maxlogp_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))
  window_minlogp_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))
  window_meanlogp_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))
  window_meanpos_window3_List<-rep(list(matrix(0, nrow = 1, ncol = 3,byrow = TRUE)),length(which_chroms))
  
  
  
  for (j in 1:length (which_chroms)){
    
    window_maxlogp_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,5]-1),byrow = TRUE)
    window_minlogp_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,5]-1),byrow = TRUE)
    window_meanlogp_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,5]-1),byrow = TRUE)
    window_meanpos_window3_List[[j]]<-matrix(0, nrow = 1, ncol = (chromo_param[j,5]-1),byrow = TRUE)
    
    the_pos <- pos_List2[[j]]
    
    this_chrom_logp<-pval_log_List2[[j]][1:(window_sizes[3] * (chromo_param[j,5]-1))]
    le_window<-window_sizes[3]
    numb_le_window<-(chromo_param[j,5]-1)
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
  
  date()
  
 for(sd in 11:20){ 
  #---------- HMM Time------------------
  
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
  
  # now do reverse fitting too ---------------------------------------------------------------------
  
  
 # date()
  
 # Mix_mod_List_rev<-vector("list", (length(max_states)))
 # Fit_Mix_mod_List_rev<-vector("list", (length(max_states)))
  
 # for (pp in 1:max_states){
 #   set.seed(123)
    
 #   Mix_mod_List_rev_sublist<-vector("list", (length(which_chroms)))
 #   Fit_Mix_mod_List_rev_sublist<-vector("list", (length(which_chroms)))
    
 #   for(i in 1:length(which_chroms)){
 #     Mix_mod_List_rev_sublist[[i]]<-depmix(window_meanlogp_window3_List_rev[[i]] ~ 1,nstates = pp, family = gaussian(link = "inverse"), ntimes=length((window_meanlogp_window3_List_rev[[i]])))
 #     Fit_Mix_mod_List_rev_sublist[[i]]<-fit(Mix_mod_List_rev_sublist[[i]],em=em.control(maxit=2500000))
 #   }
 #   Mix_mod_List_rev[[pp]]<-Mix_mod_List_rev_sublist
 #   Fit_Mix_mod_List_rev[[pp]]<-Fit_Mix_mod_List_rev_sublist
 # }
  
  
  
 # the_crits_list_rev<-vector("list", (length(max_states)))
  
 # for(i in 1:length(which_chroms)){
 #   the_crits <- array (NA,c((max_states),2))
 #   for (pp in 1:(max_states)){
 #     the_crits[pp,1] <- AIC((Fit_Mix_mod_List_rev[[pp]][[i]]))
 #     the_crits[pp,2] <- BIC((Fit_Mix_mod_List_rev[[pp]][[i]]))
 #   }
 #   the_crits_list_rev[[i]]<-the_crits
 # }
  
  
 # optimal_statenums_rev<-rep(0, (length(which_chroms)))
  
#  for(i in 1:length(which_chroms)){
    
 #   mod_diff <- the_crits_list_rev[[i]][1:(max_states - 1),2] - the_crits_list_rev[[i]][2:(max_states),2]
 #   diff_neg <- mod_diff < 0
    
 #   if (sum (diff_neg) == 0){
 #     best_mod <- max_states	
 #   } else {
 #     best_mod <- min (which (diff_neg)) #find the first value that is negative
 #   }	
 #   optimal_statenums_rev[i]<-best_mod
 # }
  
  
  #fitting of a bunch of relevant state numbers forward-----------------------------------
  
  
  
 # testing_statenums<-sort(unique(c(optimal_statenums,optimal_statenums_rev)),decreasing=FALSE)
  testing_statenums2<-c(2,max_states)
  #testing_statenums2
  
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
  
  
 # the_crits_list_rev
  
  
  #fitting of a bunch of relevant state numbers reverse-----------------------------------
  
  
 # optimal_Mix_mod_List_rev<-vector("list", (length(which_chroms)))
 # optimal_Fit_Mix_mod_List_rev<-vector("list", (length(which_chroms)))
  
#  for(j in 1:length(which_chroms)){
#    optimal_Mix_mod_List_rev[[j]]<-Mix_mod_List_rev[[(optimal_statenums[j])]][[j]]
#    optimal_Fit_Mix_mod_List_rev[[j]]<-Fit_Mix_mod_List_rev[[(optimal_statenums[j])]][[j]]
#  }
  
  
  
#  est.states_rev_variablestatenum<-vector("list", (testing_statenums2[2]-testing_statenums2[1]+1))
  
#  for (k in testing_statenums2[1]:testing_statenums2[2]){
    
#    est.states_rev_variablestatenum_sublist<-vector("list", (length(which_chroms)))
    
#    for(i in 1:length(which_chroms)){
#      est.states_rev_variablestatenum_sublist[[i]]<-(posterior(Fit_Mix_mod_List_rev[[k]][[i]]))$state
#    }
#    est.states_rev_variablestatenum[[(k-testing_statenums2[1]+1)]]<-est.states_rev_variablestatenum_sublist
    
#  }
  
  
#  Fitted_rev_Item_reorderstates_variablestatenum<-vector("list", (testing_statenums2[2]-testing_statenums2[1]+1))
#  Fitted_rev_Item_reorderstates_index_variablestatenum<-vector("list", (testing_statenums2[2]-testing_statenums2[1]+1))
#  Fitted_rev_Item_reordertrans_variablestatenum<-vector("list", (testing_statenums2[2]-testing_statenums2[1]+1))
  
#  for (k in testing_statenums2[1]:testing_statenums2[2]){
    
#    Fitted_rev_Item_reorderstates<-vector("list", (length (which_chroms)))
#    Fitted_rev_Item_reorderstates_index<-vector("list", (length (which_chroms)))
#    Fitted_rev_Item_reordertrans<-vector("list", (length (which_chroms)))
    
    
#    for(j in 1:length (which_chroms)){
      
#      h_states<-getpars(Fit_Mix_mod_List_rev[[k]][[j]])[seq((k+k*k+1),(k+k*k+2*k),2)]
#      names(h_states)<-seq(1,k,1)
#      h_states2<-h_states[order(h_states,decreasing=TRUE)]
#      h_states2_index<-order(h_states,decreasing=TRUE)
#      
#      Fitted_rev_Item_reorderstates[[j]]<-h_states2
#      Fitted_rev_Item_reorderstates_index[[j]]<-h_states2_index
      
#      h_states_trans<-t(array((getpars(Fit_Mix_mod_List_rev[[k]][[j]])[seq((k+1),(k+k*k),1)]),c(k,k)))
#      index_statenames<-seq(1,k,1)
#      rownames(h_states_trans) <- index_statenames
#      colnames(h_states_trans) <- index_statenames
#      h_states_trans_reorder<-h_states_trans[h_states2_index,h_states2_index]
      
#      Fitted_rev_Item_reordertrans[[j]]<-h_states_trans_reorder
      
#    }
    
#    Fitted_rev_Item_reorderstates_variablestatenum[[(k-testing_statenums2[1]+1)]]<-Fitted_rev_Item_reorderstates
#    Fitted_rev_Item_reorderstates_index_variablestatenum[[(k-testing_statenums2[1]+1)]]<-Fitted_rev_Item_reorderstates_index
#    Fitted_rev_Item_reordertrans_variablestatenum[[(k-testing_statenums2[1]+1)]]<-Fitted_rev_Item_reordertrans
    
#  }
  
  
  
 # est.states_rev_variablestatenum_reorder<-vector("list", (testing_statenums2[2]-testing_statenums2[1]+1))
  
 # for(k in testing_statenums2[1]:testing_statenums2[2]){
    
 #   est.states_rev_variablestatenum_reorder_sub<-vector("list", (length (which_chroms)))
    
 #   for(j in 1:length(which_chroms)){
 #     est.states_rev_variablestatenum_reorder_sub[[j]]<-est.states_rev_variablestatenum[[(k-testing_statenums2[1]+1)]][[j]]
 #     for(i in 1:length ((Fitted_rev_Item_reorderstates_index_variablestatenum[[(k-testing_statenums2[1]+1)]][[j]]))){
 #       est.states_rev_variablestatenum_reorder_sub[[j]][(est.states_rev_variablestatenum_reorder_sub[[j]])==(Fitted_rev_Item_reorderstates_index_variablestatenum[[(k-testing_statenums2[1]+1)]][[j]][i])] <- (100+i)
 #     }
 #     est.states_rev_variablestatenum_reorder_sub[[j]]<-est.states_rev_variablestatenum_reorder_sub[[j]]-100
 #   }
 #   est.states_rev_variablestatenum_reorder[[(k-testing_statenums2[1]+1)]]<-est.states_rev_variablestatenum_reorder_sub
    
 # }
  
  
  
 # est.states_rev_variablestatenum_reorder_unreversed<-vector("list", (testing_statenums2[2]-testing_statenums2[1]+1))
  
 # for(k in 1:length(est.states_rev_variablestatenum_reorder_unreversed)){
 #   est.states_rev_variablestatenum_reorder_unreversed[[k]]<-est.states_rev_variablestatenum_reorder[[k]]
 #   for(j in 1:length (which_chroms)){
 #     est.states_rev_variablestatenum_reorder_unreversed[[k]][[j]]<-rev(est.states_rev_variablestatenum_reorder[[k]][[j]])
 #   }
 # }
  
  
  
#  cluster_est.states_rev_variablestatenum_reorder_unreversed<-vector("list", (length (est.states_rev_variablestatenum_reorder_unreversed)))
  
#  for(k in 1:length(est.states_rev_variablestatenum_reorder_unreversed)){
#    cluster_est.states_rev_variablestatenum_reorder_unreversed_sub<-vector("list", (length (which_chroms)))
    
#    for(j in 1:length (which_chroms)){
#      cluster_est.states_rev_variablestatenum_reorder_unreversed_sub[[j]] <- cluster_hmm(est.states_rev_variablestatenum_reorder_unreversed[[k]][[j]])
#    }
#    cluster_est.states_rev_variablestatenum_reorder_unreversed[[k]]<-cluster_est.states_rev_variablestatenum_reorder_unreversed_sub
#  }
  
  
  
  #window stats tables-----------------------------------------------------------------
  
  Allwindow_stat_List<-vector("list", (length (which_chroms)))
  
  
  for(j in 1:length (which_chroms)){
    Allwindow_stat_List[[j]]<- array (0,c((length(window_meanlogp_window3_List[[j]])),(8+25*(length(Fitted_Item_reorderstates_variablestatenum)))))
    
    windowpos<-seq(1,(length(window_meanlogp_window3_List[[j]])),1) 
    Allwindow_stat_List[[j]][,1]<-windowpos
    
    for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
      
      Allwindow_stat_List[[j]][i,2]<-pos_List2[[j]][(Allwindow_stat_List[[j]][i,1]-1)*window_sizes[3]+1]
      Allwindow_stat_List[[j]][i,3]<-pos_List2[[j]][(Allwindow_stat_List[[j]][i,1])*window_sizes[3]]
      Allwindow_stat_List[[j]][i,4]<-(pos_List2[[j]][(Allwindow_stat_List[[j]][i,1]-1)*window_sizes[3]+1]+pos_List2[[j]][(Allwindow_stat_List[[j]][i,1])*window_sizes[3]])/2
      
      
      for(k in 1:length(Fitted_Item_reorderstates_variablestatenum)){
        Allwindow_stat_List[[j]][i,7+4*(k-1)+1]<-est.states_variablestatenum_reorder[[k]][[j]][i]
        Allwindow_stat_List[[j]][i,7+4*(k-1)+2]<-cluster_est.states_variablestatenum_reorder[[k]][[j]][i]
 #       Allwindow_stat_List[[j]][i,7+4*(k-1)+3]<-est.states_rev_variablestatenum_reorder_unreversed[[k]][[j]][i]
 #       Allwindow_stat_List[[j]][i,7+4*(k-1)+4]<-cluster_est.states_rev_variablestatenum_reorder_unreversed[[k]][[j]][i]
        
      }
    }
    
    window_mean<-window_meanlogp_window3_List[[j]]
    Allwindow_stat_List[[j]][,5]<-window_mean
    
    window_max<-window_maxlogp_window3_List[[j]]
    Allwindow_stat_List[[j]][,6]<-window_max
    
    window_min<-window_minlogp_window3_List[[j]]
    Allwindow_stat_List[[j]][,7]<-window_min
    
    
  }
  
  
  
  #get min of peak state and max of all non-peak states for each chromosomes
  minpeakANDmaxnonpeak_thres<-array (0,c((length (which_chroms)),4*length(Fitted_Item_reorderstates_variablestatenum)))
  
  for(j in 1:length(which_chroms)){
    
    for(k in 1:length(Fitted_Item_reorderstates_variablestatenum)){
      Temp_states<-sort(unique(Allwindow_stat_List[[j]][,(7+4*(k-1)+1)]),decreasing=FALSE)
    #  Temp_states_rev<-sort(unique(Allwindow_stat_List[[j]][,(7+4*(k-1)+3)]),decreasing=FALSE)
      
      Temp_minmean_peakwindow<-min(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,(7+4*(k-1)+1)] %in% c(1)),5])
      Temp_maxmean_nonpeakwindow<-max(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,(7+4*(k-1)+1)] %in% (Temp_states[2:length(Temp_states)])),5])
      
  #    Temp_minmean_peakwindow_rev<-min(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,(7+4*(k-1)+3)] %in% c(1)),5])
  #    Temp_maxmean_nonpeakwindow_rev<-max(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,(7+4*(k-1)+3)] %in% (Temp_states_rev[2:length(Temp_states_rev)])),5])
      
      minpeakANDmaxnonpeak_thres[j,(4*(k-1)+1)]<-Temp_minmean_peakwindow
      minpeakANDmaxnonpeak_thres[j,(4*(k-1)+2)]<-Temp_maxmean_nonpeakwindow
   #   minpeakANDmaxnonpeak_thres[j,(4*(k-1)+3)]<-Temp_minmean_peakwindow_rev
   #   minpeakANDmaxnonpeak_thres[j,(4*(k-1)+4)]<-Temp_maxmean_nonpeakwindow_rev
    }
  }
  
  
  
  
  #lets fill in I_window, I_window_rev, I_window3, I_window3_rev variables
  
  for(j in 1:length (which_chroms)){
    for(k in 1:length(Fitted_Item_reorderstates_variablestatenum)){
      #I_window
      for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
        Allwindow_stat_List[[j]][i,(7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<-0
        if(Allwindow_stat_List[[j]][i,(7+4*(k-1)+1)]==1){
          Allwindow_stat_List[[j]][i,(7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<-1
        }
      }
      
#      I_window_rev
   # for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
   #     Allwindow_stat_List[[j]][i,(7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]<-0
   #     if(Allwindow_stat_List[[j]][i,(7+4*(k-1)+3)]==1){
   #       Allwindow_stat_List[[j]][i,(7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]<-1
   #     }
   #   }
      
      #I_window and I_window_rev intersection
    #  for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
    #    Allwindow_stat_List[[j]][i,(7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]<-0
    #   if(Allwindow_stat_List[[j]][i,(7+4*(k-1)+1)]==1 & Allwindow_stat_List[[j]][i,(7+4*(k-1)+3)]==1){
    #      Allwindow_stat_List[[j]][i,(7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]<-1
    #    }
    #  }
      
      
      #I_window3
  #    for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  #      if(Allwindow_stat_List[[j]][i,(7+4*(k-1)+1)]==1 & Allwindow_stat_List[[j]][i,5]>minpeakANDmaxnonpeak_thres[j,(4*(k-1)+2)]){
  #        Allwindow_stat_List[[j]][i,(7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<-1
  #      }
  #    }
      
      #I_window3_rev
  #    for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  #      if(Allwindow_stat_List[[j]][i,(7+4*(k-1)+3)]==1 & Allwindow_stat_List[[j]][i,5]>minpeakANDmaxnonpeak_thres[j,(4*(k-1)+4)]){
  #        Allwindow_stat_List[[j]][i,(7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]<-1
  #      }
  #    }
      
      #I_window3 and I_window3_rev intersection
  #    for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  #      if(Allwindow_stat_List[[j]][i,(7+4*(k-1)+1)]==1 & Allwindow_stat_List[[j]][i,5]>minpeakANDmaxnonpeak_thres[j,(4*(k-1)+2)] & Allwindow_stat_List[[j]][i,(7+4*(k-1)+3)]==1 & Allwindow_stat_List[[j]][i,5]>minpeakANDmaxnonpeak_thres[j,(4*(k-1)+4)]){
  #        Allwindow_stat_List[[j]][i,(7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]<-1
  #      }
  #    }
      
    } 
    
    
    #log10 of mean -log10 p value
   # for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
   #   Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1)]<-log10(Allwindow_stat_List[[j]][i,5]+10^(-4))
   # }
    
    
    for(k in 1:length(Fitted_Item_reorderstates_variablestatenum)){
      #phybrid
      
      #log10 of mean -log10 p value adjustment to I_window
     # for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
     #   Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(k-1)+1)]<-(Allwindow_stat_List[[j]][i,(7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]*Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1)])
     # }
      
      #log10 of mean -log10 p value adjustment to I_window_rev
     # for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
     #   Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(k-1)+2)]<-(Allwindow_stat_List[[j]][i,(7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]*Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1)])
     # }
      
      #log10 of mean -log10 p value adjustment to I_window and I_window_rev intersection
   #   for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  #      Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(k-1)+3)]<-(Allwindow_stat_List[[j]][i,(7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]*Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1)])
  #    }
      
      #log10 of mean -log10 p value adjustment to I_window3
  #    for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  #      Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<-(Allwindow_stat_List[[j]][i,(7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]*Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1)])
  #    }
      
      #log10 of mean -log10 p value adjustment to I_window3_rev
  #    for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  #      Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]<-(Allwindow_stat_List[[j]][i,(7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]*Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1)])
  #    }
      
      #log10 of mean -log10 p value adjustment to I_window3 and I_window3_rev intersection
  #    for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  #      Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]<-(Allwindow_stat_List[[j]][i,(7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]*Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1)])
  #    }
      
      
      #log10 of mean -log10 p value adjustment to I_window (no negative)
 #     for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
#        Allwindow_stat_List[[j]][i,(8+16*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<-Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(k-1)+1)]
 #       if(Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(k-1)+1)]<0){
  #        Allwindow_stat_List[[j]][i,(8+16*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<-0
  #      }
  #    }
      
      #log10 of mean -log10 p value adjustment to I_window_rev (no negative)
  #    for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  #      Allwindow_stat_List[[j]][i,(8+16*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]<-Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(k-1)+2)]
  #      if(Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(k-1)+2)]<0){
  #        Allwindow_stat_List[[j]][i,(8+16*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]<-0
  #      }
  #    }
      
      #log10 of mean -log10 p value adjustment to I_window and I_window_rev intersection (no negative)
 #     for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
#        Allwindow_stat_List[[j]][i,(8+16*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]<-Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(k-1)+3)]
#        if(Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(k-1)+3)]<0){
#          Allwindow_stat_List[[j]][i,(8+16*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]<-0
#        }
#      }
      
      
      #log10 of mean -log10 p value adjustment to I_window3 (no negative)
  #    for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  #      Allwindow_stat_List[[j]][i,(8+19*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<-Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]
  #      if(Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<0){
  #        Allwindow_stat_List[[j]][i,(8+19*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<-0
  #      }
  #    }
      
      
      #log10 of mean -log10 p value adjustment to I_window3_rev (no negative)
  #    for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  ##      Allwindow_stat_List[[j]][i,(8+19*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]<-Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]
  #      if(Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]<0){
  #        Allwindow_stat_List[[j]][i,(8+19*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]<-0
  #      }
  #    }
      
      
      #log10 of mean -log10 p value adjustment to I_window3 and I_window3_rev intersection (no negative)
 #     for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
#        Allwindow_stat_List[[j]][i,(8+19*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]<-Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]
 #       if(Allwindow_stat_List[[j]][i,(7+10*(length(Fitted_Item_reorderstates_variablestatenum))+1+3*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]<0){
  #        Allwindow_stat_List[[j]][i,(8+19*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]<-0
  #      }
  #    }
      
      
      
      #I_window5
  #    for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  #      Allwindow_stat_List[[j]][i,(8+22*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<-0
  #      if(Allwindow_stat_List[[j]][i,(7+4*(k-1)+1)]<3){
  #        Allwindow_stat_List[[j]][i,(8+22*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+1)]<-1
  #      }
  #    }
      
      #I_window5_rev
   #   for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  #      Allwindow_stat_List[[j]][i,(8+22*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]<-0
  #      if(Allwindow_stat_List[[j]][i,(7+4*(k-1)+3)]<3){
  #        Allwindow_stat_List[[j]][i,(8+22*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+2)]<-1
  #      }
  #    }
      
      #I_window5 and I_window5_rev intersection
  #    for(i in 1:(length(window_meanlogp_window3_List[[j]]))){
  #      Allwindow_stat_List[[j]][i,(8+22*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]<-0
  #      if(Allwindow_stat_List[[j]][i,(7+4*(k-1)+1)]<3 & Allwindow_stat_List[[j]][i,(7+4*(k-1)+3)]<3){
  #        Allwindow_stat_List[[j]][i,(8+22*(length(Fitted_Item_reorderstates_variablestatenum))+3*(k-1)+3)]<-1
  #      }
  #    }
      
    }
    
  }
  
  
  #Moran's I window tables-----------------------------------
  
  
  I_windows_results<-vector("list", (length (which_chroms)))
 # I_windows_results_3<-vector("list", (length (which_chroms)))
#  I_windows_results_phybrid<-vector("list", (length (which_chroms)))
#  I_windows_results_phybrid3<-vector("list", (length (which_chroms)))
#  I_windows_results_phybridnonneg<-vector("list", (length (which_chroms)))
#  I_windows_results_phybridnonneg3<-vector("list", (length (which_chroms)))
#  I_windows_results_5<-vector("list", (length (which_chroms)))
  
  
 # for(j in 1:length (which_chroms)){
    
#    dist_bp_mat_windows<-array (0,c((dim(Allwindow_stat_List[[j]])[1]),(dim(Allwindow_stat_List[[j]])[1])))
    
#    for(k in 1:((dim(Allwindow_stat_List[[j]])[1])-1)){
#      for(i in (k+1):(dim(Allwindow_stat_List[[j]])[1])){
#        dist_bp_mat_windows[k,i]<-(Allwindow_stat_List[[j]][i,4]-Allwindow_stat_List[[j]][k,4])
#      }
#    }
    #bottom left of diagonal
#    for(k in 2:((dim(Allwindow_stat_List[[j]])[1]))){
#      for(i in 1:(k-1)){
#        dist_bp_mat_windows[k,i]<-(Allwindow_stat_List[[j]][k,4]-Allwindow_stat_List[[j]][i,4])
#      }
#    }
    
#    dist_bp_mat_windows<-as.matrix(dist_bp_mat_windows)
#    dist_bp_mat_windows_weights<-1/dist_bp_mat_windows
#    diag(dist_bp_mat_windows_weights)<-0
    
    
    
#    I_windows_results_sub<-array (0,c((length(Fitted_Item_reorderstates_variablestatenum)),(7)))
    
#    for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
#      I_windows_results_sub[n,1]<-(Moran.I((Allwindow_stat_List[[j]][,7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+1]), dist_bp_mat_windows_weights))$observed
#      I_windows_results_sub[n,2]<-(Moran.I((Allwindow_stat_List[[j]][,7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+2]), dist_bp_mat_windows_weights))$observed
#      I_windows_results_sub[n,3]<-(Moran.I((Allwindow_stat_List[[j]][,7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3]), dist_bp_mat_windows_weights))$observed
#    }
    
#    I_windows_results[[j]]<-I_windows_results_sub
    
    
#    I_windows_results_3sub<-array (0,c((length(Fitted_Item_reorderstates_variablestatenum)),(7)))
    
#    for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
#      if(mean((Allwindow_stat_List[[j]][,7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+1]))>0){
#        I_windows_results_3sub[n,1]<-(Moran.I((Allwindow_stat_List[[j]][,7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+1]), dist_bp_mat_windows_weights))$observed
#      }
      
#      if(mean((Allwindow_stat_List[[j]][,7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+2]))>0){
#        I_windows_results_3sub[n,2]<-(Moran.I((Allwindow_stat_List[[j]][,7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+2]), dist_bp_mat_windows_weights))$observed
#      }
      
#      if(mean((Allwindow_stat_List[[j]][,7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3]))>0){
#        I_windows_results_3sub[n,3]<-(Moran.I((Allwindow_stat_List[[j]][,7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3]), dist_bp_mat_windows_weights))$observed
#      }
#    }
    
#    I_windows_results_3[[j]]<-I_windows_results_3sub
    
    
    
#    I_windows_results_phybridsub<-array (0,c((length(Fitted_Item_reorderstates_variablestatenum)),(3)))
    
#    for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
#      I_windows_results_phybridsub[n,1]<-(Moran.I((Allwindow_stat_List[[j]][,8+10*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3]), dist_bp_mat_windows_weights))$observed
#    }
#    I_windows_results_phybrid[[j]]<-I_windows_results_phybridsub
    
    
    
#    I_windows_results_phybrid3sub<-array (0,c((length(Fitted_Item_reorderstates_variablestatenum)),(3)))
    
#    for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
#      if(mean((Allwindow_stat_List[[j]][,8+13*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3]))>0){
#        I_windows_results_phybrid3sub[n,1]<-(Moran.I((Allwindow_stat_List[[j]][,8+13*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3]), dist_bp_mat_windows_weights))$observed
#      }
      
#    }
#    I_windows_results_phybrid3[[j]]<-I_windows_results_phybrid3sub
    
    
    
 #   I_windows_results_phybridnonnegsub<-array (0,c((length(Fitted_Item_reorderstates_variablestatenum)),(3)))
    
#    for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
#      I_windows_results_phybridnonnegsub[n,1]<-(Moran.I((Allwindow_stat_List[[j]][,8+16*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3]), dist_bp_mat_windows_weights))$observed
#    }
#    I_windows_results_phybridnonneg[[j]]<-I_windows_results_phybridnonnegsub
    
    
    
#    I_windows_results_phybridnonneg3sub<-array (0,c((length(Fitted_Item_reorderstates_variablestatenum)),(3)))
    
#    for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
      
#      if(mean((Allwindow_stat_List[[j]][,8+19*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3]))>0){
#        I_windows_results_phybridnonneg3sub[n,1]<-(Moran.I((Allwindow_stat_List[[j]][,8+19*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3]), dist_bp_mat_windows_weights))$observed
#      }
#    }
#    I_windows_results_phybridnonneg3[[j]]<-I_windows_results_phybridnonneg3sub
    
    
    
#    I_windows_results_5sub<-array (0,c((length(Fitted_Item_reorderstates_variablestatenum)),(5)))
    
#    for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
#      I_windows_results_5sub[n,1]<-(Moran.I((Allwindow_stat_List[[j]][,8+22*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3]), dist_bp_mat_windows_weights))$observed
#    }
    
#    I_windows_results_5[[j]]<-I_windows_results_5sub
    
    
    
#  }
  
 # dist_bp_mat_windows<-array (0,c((2),(2)))
  
  
  
  
  for(j in 1:length (which_chroms)){
    
    
   # for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
   #   I_windows_results[[j]][n,4]<- mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] %in% c(1)),5])
   #   I_windows_results[[j]][n,5]<- (mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] %in% c(1)),5])/mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,7+4*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] < 1),5]))
   #   I_windows_results[[j]][n,6]<-I_windows_results[[j]][n,4]*I_windows_results[[j]][n,3]
   #   I_windows_results[[j]][n,7]<-I_windows_results[[j]][n,5]*I_windows_results[[j]][n,3]
   # }
    
    
 #   for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
#      if(mean((Allwindow_stat_List[[j]][,7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3]))>0){
#        I_windows_results_3[[j]][n,4]<- mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] %in% c(1)),5])
 #       I_windows_results_3[[j]][n,5]<- (mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] %in% c(1)),5])/mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,7+7*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] < 1),5]))
#        I_windows_results_3[[j]][n,6]<-I_windows_results_3[[j]][n,4]*I_windows_results_3[[j]][n,3]
#        I_windows_results_3[[j]][n,7]<-I_windows_results_3[[j]][n,5]*I_windows_results_3[[j]][n,3]
#      }
#    }
    
    
 #   for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
#      I_windows_results_phybrid[[j]][n,2]<-I_windows_results_phybrid[[j]][n,1]*I_windows_results[[j]][n,4]
#      I_windows_results_phybrid[[j]][n,3]<-I_windows_results_phybrid[[j]][n,1]*I_windows_results[[j]][n,5]
#    }
    
    
#    for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
#      I_windows_results_phybrid3[[j]][n,2]<-I_windows_results_phybrid3[[j]][n,1]*I_windows_results_3[[j]][n,4]
#      I_windows_results_phybrid3[[j]][n,3]<-I_windows_results_phybrid3[[j]][n,1]*I_windows_results_3[[j]][n,5]
#    }
    
    
 #   for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
      
#      I_windows_results_phybridnonneg[[j]][n,2]<-I_windows_results_phybridnonneg[[j]][n,1]*mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,8+10*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] > 0),5])
#      I_windows_results_phybridnonneg[[j]][n,3]<-I_windows_results_phybridnonneg[[j]][n,1]*(mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,8+10*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] > 0),5])/mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,8+10*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] <= 0),5]))
#    }
    
    
 #   for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
      
  #    if(mean((Allwindow_stat_List[[j]][,8+13*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3]))>0){
  #      I_windows_results_phybridnonneg3[[j]][n,2]<-I_windows_results_phybridnonneg3[[j]][n,1]*mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,8+13*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] > 0),5])
  #      I_windows_results_phybridnonneg3[[j]][n,3]<-I_windows_results_phybridnonneg3[[j]][n,1]*(mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,8+13*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] > 0),5])/mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,8+13*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] <= 0),5]))
  #    }
      
      
  # }
    
    
    
    
 #   for(n in 1:length(Fitted_Item_reorderstates_variablestatenum)){
#      I_windows_results_5[[j]][n,2]<- mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,8+22*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] %in% c(1)),5])
#      I_windows_results_5[[j]][n,3]<- (mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,8+22*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] %in% c(1)),5])/mean(Allwindow_stat_List[[j]][which(Allwindow_stat_List[[j]][,8+22*(length(Fitted_Item_reorderstates_variablestatenum))+3*(n-1)+3] < 1),5]))
#      I_windows_results_5[[j]][n,4]<-I_windows_results_5[[j]][n,2]*I_windows_results_5[[j]][n,1]
#      I_windows_results_5[[j]][n,5]<-I_windows_results_5[[j]][n,3]*I_windows_results_5[[j]][n,1]
#    }
    
  }
  
  
  
#  optimal_statenums_store<-optimal_statenums-testing_statenums2[1]+1
  
  
  
  I_window_optimal_StatTable<-array (0,c((length(which_chroms)),(31)))
  
  for(j in 1:length(which_chroms)){
  #  I_window_optimal_StatTable[j,1]<-I_windows_results[[j]][optimal_statenums_store[j],1]
  #  I_window_optimal_StatTable[j,2]<-I_windows_results[[j]][optimal_statenums_store[j],2]
   # I_window_optimal_StatTable[j,3]<-I_windows_results[[j]][optimal_statenums_store[j],3]
   # I_window_optimal_StatTable[j,4]<-I_windows_results[[j]][optimal_statenums_store[j],4]
   # I_window_optimal_StatTable[j,5]<-I_windows_results[[j]][optimal_statenums_store[j],5]
   # I_window_optimal_StatTable[j,6]<-I_windows_results[[j]][optimal_statenums_store[j],6]
   # I_window_optimal_StatTable[j,7]<-I_windows_results[[j]][optimal_statenums_store[j],7]
    
   # I_window_optimal_StatTable[j,8]<-I_windows_results_3[[j]][optimal_statenums_store[j],1]
  #  I_window_optimal_StatTable[j,9]<-I_windows_results_3[[j]][optimal_statenums_store[j],2]
  #  I_window_optimal_StatTable[j,10]<-I_windows_results_3[[j]][optimal_statenums_store[j],3]
  #  I_window_optimal_StatTable[j,11]<-I_windows_results_3[[j]][optimal_statenums_store[j],4]
  ##  I_window_optimal_StatTable[j,12]<-I_windows_results_3[[j]][optimal_statenums_store[j],5]
  #  I_window_optimal_StatTable[j,13]<-I_windows_results_3[[j]][optimal_statenums_store[j],6]
  #  I_window_optimal_StatTable[j,14]<-I_windows_results_3[[j]][optimal_statenums_store[j],7]
    
  #  I_window_optimal_StatTable[j,15]<-I_windows_results_phybrid[[j]][optimal_statenums_store[j],1]
  #  I_window_optimal_StatTable[j,16]<-I_windows_results_phybrid[[j]][optimal_statenums_store[j],2]
  #  I_window_optimal_StatTable[j,17]<-I_windows_results_phybrid[[j]][optimal_statenums_store[j],3]
    
  #  I_window_optimal_StatTable[j,18]<-I_windows_results_phybrid3[[j]][optimal_statenums_store[j],1]
  #  I_window_optimal_StatTable[j,19]<-I_windows_results_phybrid3[[j]][optimal_statenums_store[j],2]
  #  I_window_optimal_StatTable[j,20]<-I_windows_results_phybrid3[[j]][optimal_statenums_store[j],3]
    
  #  I_window_optimal_StatTable[j,21]<-I_windows_results_phybridnonneg[[j]][optimal_statenums_store[j],1]
  #  I_window_optimal_StatTable[j,22]<-I_windows_results_phybridnonneg[[j]][optimal_statenums_store[j],2]
  #  I_window_optimal_StatTable[j,23]<-I_windows_results_phybridnonneg[[j]][optimal_statenums_store[j],3]
    
  #  I_window_optimal_StatTable[j,24]<-I_windows_results_phybridnonneg3[[j]][optimal_statenums_store[j],1]
  #  I_window_optimal_StatTable[j,25]<-I_windows_results_phybridnonneg3[[j]][optimal_statenums_store[j],2]
  #  I_window_optimal_StatTable[j,26]<-I_windows_results_phybridnonneg3[[j]][optimal_statenums_store[j],3]
    
  #  I_window_optimal_StatTable[j,27]<-I_windows_results_5[[j]][optimal_statenums_store[j],1]
  #  I_window_optimal_StatTable[j,28]<-I_windows_results_5[[j]][optimal_statenums_store[j],2]
  #  I_window_optimal_StatTable[j,29]<-I_windows_results_5[[j]][optimal_statenums_store[j],3]
  #  I_window_optimal_StatTable[j,30]<-I_windows_results_5[[j]][optimal_statenums_store[j],4]
  #  I_window_optimal_StatTable[j,31]<-I_windows_results_5[[j]][optimal_statenums_store[j],5]
  }
  
  
 # I_window_optimal_StatTable_reorder<-I_window_optimal_StatTable[c(1),]
  
  
  
  #data for export ------------------------------------
  
  
  saveRDS(Allwindow_stat_List, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_Allwindow_stat_List_Streamline_seed_",(1000+sd),".rds"),collapse=""))
 # saveRDS(I_windows_results, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_I_windows_results_Streamline.rds"),collapse=""))
 # saveRDS(I_windows_results_3, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_I_windows_results_3.rds"),collapse=""))
#  saveRDS(I_windows_results_phybrid, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_I_windows_results_phybrid.rds"),collapse=""))
#  saveRDS(I_windows_results_phybrid3, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_I_windows_results_phybrid3.rds"),collapse=""))
  #saveRDS(I_windows_results_phybridnonneg, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_I_windows_results_phybridnonneg.rds"),collapse=""))
  #saveRDS(I_windows_results_phybridnonneg3, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_I_windows_results_phybridnonneg3.rds"),collapse=""))
  #saveRDS(I_windows_results_5, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_I_windows_results_5.rds"),collapse=""))
  
#  saveRDS(optimal_statenums, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_optimal_statenums_Streamline.rds"),collapse=""))
#  saveRDS(optimal_statenums_rev, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_optimal_statenums_rev_Streamline.rds"),collapse=""))
  saveRDS(testing_statenums2, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_testing_statenums2_Streamline_seed_",(1000+sd),".rds"),collapse=""))
  
  saveRDS(the_crits_list, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_the_crits_list_Streamline_seed_",(1000+sd),".rds"),collapse=""))
#  saveRDS(the_crits_list_rev, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_the_crits_list_rev_Streamline.rds"),collapse=""))

  
 # write.csv(I_window_optimal_StatTable, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_I_window_optimal_StatTable_Streamline.csv"),collapse=""), row.names=FALSE)
 # write.csv(I_window_optimal_StatTable_reorder, paste(c("/data/michael_HMM/HMM_SNPfiles/New_folder/220927_Capsella_rubella_Weigel_Individual/",Enviro_Vars_of_Interest,"_ReOrder_window_",window_sizes[3],"_I_window_optimal_StatTable_reorder_Streamline.csv"),collapse=""), row.names=FALSE)

}
  
}




