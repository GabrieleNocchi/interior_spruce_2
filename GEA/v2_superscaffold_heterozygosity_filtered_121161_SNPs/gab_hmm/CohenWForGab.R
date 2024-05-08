library(tidyverse)
library(GlmSimulatoR)
library(ggplot2)
library(stats)
library(depmixS4)
library(fitdistrplus)
library(spdep)
library(ape)
library(actuar)
library(gplots)
library(viridisLite)
library(gridExtra)


ReplicateSeeds<-2
bumper_size<-0
ws<-1
mxs<-15
files <- list.files(pattern = "*_window_stats_seed_1001\\.rds")



for (file in files) {
file_name <- gsub("_ReOrder_window_1_window_stats_seed_1001.rds", "", file)

  
  All_Window_MassRep<-list()
  
  All_Window_Temp<-readRDS(file)
 
  ch_num<-length(All_Window_Temp)
  
  
  for(kk in 1:ch_num){
  
  All_Window_MassRep[[kk]]<-list()
  
  for(i in 1:(mxs-1)){
    All_Window_MassRep[[kk]][[i]]<-as.data.frame(array(NA,c((dim(All_Window_Temp[[kk]])[1]),(5+2*ReplicateSeeds))))
    All_Window_MassRep[[kk]][[i]][,1]<-All_Window_Temp[[kk]][,2]
    All_Window_MassRep[[kk]][[i]][,2]<-All_Window_Temp[[kk]][,3]
    All_Window_MassRep[[kk]][[i]][,3]<-All_Window_Temp[[kk]][,4]
    All_Window_MassRep[[kk]][[i]][,4]<-All_Window_Temp[[kk]][,5]
    All_Window_MassRep[[kk]][[i]][,5]<-All_Window_Temp[[kk]][,6]
    
    #HMM fit state and cluster number
    All_Window_MassRep[[kk]][[i]][,5+2*(1-1)+1]<-All_Window_Temp[[kk]][,7+4*(i-1)+1]
    All_Window_MassRep[[kk]][[i]][,5+2*(1-1)+2]<-All_Window_Temp[[kk]][,7+4*(i-1)+2]
  }
  
  
  if(ReplicateSeeds>1){
  for(jj in 2:ReplicateSeeds){
    All_Window_Temp<-readRDS(paste(file_name,"_ReOrder_window_",ws, "_window_stats_seed_",(1000+jj),".rds",sep=""))
    
    for(i in 1:(mxs-1)){
      All_Window_MassRep[[kk]][[i]][,5+2*(jj-1)+1]<-All_Window_Temp[[kk]][,7+4*(i-1)+1]
      All_Window_MassRep[[kk]][[i]][,5+2*(jj-1)+2]<-All_Window_Temp[[kk]][,7+4*(i-1)+2]
    }
  }
  }
  
  }


  
  
  
  All_Window_MassRep_HMMStates<-list()
  
  for(kk in 1:length(All_Window_MassRep)){
    All_Window_MassRep_HMMStates[[kk]]<-list()
    
  for(i in 1:(mxs-1)){
    All_Window_MassRep_HMMStates[[kk]][[i]]<-All_Window_MassRep[[kk]][[i]][,seq((5+1),(2*(ReplicateSeeds-1)+5+1),2)]
  }
  }
  
  


  All_Window_MassRep_HMMStates_PeakFilter<-All_Window_MassRep_HMMStates
  
  for(kk in 1:length(All_Window_MassRep_HMMStates)){
    All_Window_MassRep_HMMStates_PeakFilter[[kk]]<-list()
  for(i in 1:(mxs-1)){
    All_Window_MassRep_HMMStates_PeakFilter[[kk]][[i]]<-All_Window_MassRep_HMMStates[[kk]][[i]] %>% mutate(across(where(is.numeric), function(x) ifelse(x > 1, 0, x)))
  }
  }

  
 
date()
  




peakclusters_List<-list()
peakclusters_StartEndList<-list()
peakclusters_size<-list()


for(kk in 1:length(All_Window_MassRep)){
  peakclusters_List[[kk]]<-list()
  peakclusters_List_Temp<-rep(list(),length(All_Window_MassRep))
  peakclusters_size[[kk]]<-list()
  
for(jj in 1:(mxs-1)){
  
  peakclusters_List_Temp[[jj]]<-rep(list(),ReplicateSeeds)
  
  for(ii in 1:ReplicateSeeds){
    peakclusters_List_TempTemp<-unique(All_Window_MassRep[[kk]][[jj]][which(All_Window_MassRep[[kk]][[jj]][,(5+2*(ii-1)+1)] %in% c(1)),(5+2*(ii-1)+2)])
    peakclusters_List_Temp[[jj]][[ii]]<-peakclusters_List_TempTemp
  }
  
  peakclusters_List[[kk]][[jj]]<- peakclusters_List_Temp[[jj]]
  
}
  
  
  for(jj in 1:(mxs-1)){
    peakclusters_size[[kk]][[jj]]<-list()
    
    for(ii in 1:ReplicateSeeds){
      
      peakclusters_size[[kk]][[jj]][[ii]]<-rep(NA,length(peakclusters_List[[kk]][[jj]][[ii]]))
      
      for(mm in 1:length(peakclusters_List[[kk]][[jj]][[ii]])){  
        peakclusters_size[[kk]][[jj]][[ii]][mm]<-max(which(All_Window_MassRep[[kk]][[jj]][,(5+2*(ii-1)+2)] %in% peakclusters_List[[kk]][[jj]][[ii]][mm]))-min(which(All_Window_MassRep[[kk]][[jj]][,(5+2*(ii-1)+2)] %in% peakclusters_List[[kk]][[jj]][[ii]][mm]))+1
      }
      
    }
  } 

  
}





#find optimal geomp and accompanying chi^2 for each combo of sim, state, replicate



opt_selfprob<-list()
opt_chistat<-list()
opt_chipvalue<-list()

opt_cW<-list()

peak_tallies_list<-list()
expec_peak_tallies_list<-list()


for(kk in 1:length(peakclusters_size)){
  
  opt_selfprob[[kk]]<-array(NA,c((mxs-1),ReplicateSeeds))
opt_chistat[[kk]]<-array(NA,c((mxs-1),ReplicateSeeds))
opt_chipvalue[[kk]]<-array(NA,c((mxs-1),ReplicateSeeds))
opt_cW[[kk]]<-array(NA,c((mxs-1),ReplicateSeeds))
  
peak_tallies_list[[kk]]<-list()
expec_peak_tallies_list[[kk]]<-list()

for(jj in 1:(mxs-1)){
  
  peak_tallies_list[[kk]][[jj]]<-list()
  expec_peak_tallies_list[[kk]][[jj]]<-list()
  
  for(ii in 1:(ReplicateSeeds)){
    
    gprob_vect<-seq(1,999)/1000
    chi_chi_baby_varygeomp<-rep(10.^10,999)
    chi_chi_pepe_varygeomp<-rep(10.^10,999)
    
    gprob<-(1-1/mean(peakclusters_size[[kk]][[jj]][[ii]]))
    lower_testbound<-max(ceiling((gprob*1000)-400)+1,1)
    upper_testbound<-min(ceiling((gprob*1000)+450)+1,999)
    
    
    peak_tally_names<-as.numeric(dimnames(table(peakclusters_size[[kk]][[jj]][[ii]]))[[1]])
    peak_tallies_pre<-table(peakclusters_size[[kk]][[jj]][[ii]])
    
    peak_tallies<-rep(0,(max(peak_tally_names)+1))
    
    for(i in 1:length(peak_tallies_pre)){
      peak_tallies[peak_tally_names[i]]<-peak_tallies_pre[i]
    }
    
    expec_templist<-list()
    
    for(i in lower_testbound:(upper_testbound)){
      expec<-dgeom(0:(max(peakclusters_size[[kk]][[jj]][[ii]])-1),(1-gprob_vect[i]))*length(peakclusters_size[[kk]][[jj]][[ii]])
      expec <- c(expec, max(0,(sum(peak_tallies)-sum(expec))))
      
      #sum(expec)
      
      y <- peak_tallies
      
      chi_chi_baby_varygeomp[i]<-as.vector((chisq.test(y,p=expec/sum(expec)))$statistic)
      chi_chi_pepe_varygeomp[i]<-as.vector((chisq.test(y,p=expec/sum(expec)))$p.value)
      expec_templist[[i]]<-expec
    }
    
    
    
    opt_selfprob[[kk]][jj,ii]<-which(chi_chi_baby_varygeomp %in% min(chi_chi_baby_varygeomp,na.rm=TRUE))/1000
    opt_chistat[[kk]][jj,ii]<-chi_chi_baby_varygeomp[which(chi_chi_baby_varygeomp %in% min(chi_chi_baby_varygeomp,na.rm=TRUE))]
    opt_chipvalue[[kk]][jj,ii]<-chi_chi_pepe_varygeomp[which(chi_chi_baby_varygeomp %in% min(chi_chi_baby_varygeomp,na.rm=TRUE))]
    
    opt_cW[[kk]][jj,ii]<-sqrt(opt_chistat[[kk]][jj,ii]/sum(expec))
    
    peak_tallies_list[[kk]][[jj]][[ii]]<-peak_tallies
    expec_peak_tallies_list[[kk]][[jj]][[ii]]<-expec_templist[[which(chi_chi_baby_varygeomp %in% min(chi_chi_baby_varygeomp,na.rm=TRUE))]]
    
    
  }
}


}








#peakcluster_array_List<-list()

#for(kk in 1:length(peakclusters_size)){
#  
#  peakcluster_array_List[[kk]]<-list()
#  
#  for(tt in 1:(mxs-1)){
#    
#    peakcluster_array_List[[kk]][[tt]]<-list()
#    
#    for(nn in 1:ReplicateSeeds){
#      
#      peakcluster_index<-unique(All_Window_MassRep[[kk]][[tt]][which(All_Window_MassRep[[kk]][[tt]][,(5+2*(nn-1)+1)] %in% c(1)),(5+2*(nn-1)+2)])
#      peakcluster_array<-array(0,c(length(peakcluster_index),5))
#      
#      for(qq in 1:length(peakcluster_index)){
#        peakcluster_array[qq,1]<-All_Window_MassRep[[kk]][[tt]][min(which(All_Window_MassRep[[kk]][[tt]][,(5+2*(nn-1)+2)] %in% peakcluster_index[qq])),2]
#        peakcluster_array[qq,2]<-All_Window_MassRep[[kk]][[tt]][max(which(All_Window_MassRep[[kk]][[tt]][,(5+2*(nn-1)+2)] %in% peakcluster_index[qq])),3]
#        
#        peakcluster_array[qq,3]<-max(which(All_Window_MassRep[[kk]][[tt]][,(5+2*(nn-1)+2)] %in% peakcluster_index[qq]))-min(which(All_Window_MassRep[[kk]][[tt]][,(5+2*(nn-1)+2)] %in% peakcluster_index[qq]))+1
#       peakcluster_array[qq,4]<-min(which(All_Window_MassRep[[kk]][[tt]][,(5+2*(nn-1)+2)] %in% peakcluster_index[qq]))
#        peakcluster_array[qq,5]<-max(which(All_Window_MassRep[[kk]][[tt]][,(5+2*(nn-1)+2)] %in% peakcluster_index[qq]))
#        
#      }
#      peakcluster_array_List[[kk]][[tt]][[nn]]<-peakcluster_array
      
#    }
#  }
#}

#date()





#saveRDS(opt_selfprob,paste("./cohen/4.Norquay_ReOrder_window_",ws,"_bumper_0_opt_selfprob.rds",sep=""))
#saveRDS(opt_chistat,paste("./cohen/4.Norquay_ReOrder_window_",ws, "_bumper_0_opt_chistat.rds",sep=""))
saveRDS(opt_chipvalue,paste("./cohen/",file_name,"_ReOrder_window_",ws, "_bumper_0_opt_chipvalue.rds",sep=""))
saveRDS(opt_cW,paste("./cohen/",file_name,"_ReOrder_window_",ws,"_bumper_0_opt_cW.rds",sep=""))
}
#saveRDS(peakclusters_List,paste("./cohen/4.Norquay_ReOrder_window_",ws,"_bumper_0_opt_cW.rds",sep=""))
#saveRDS(peakclusters_size,paste("./cohen/4.Norquay_ReOrder_window_",ws,"_bumper_0_opt_cW.rds",sep=""))


#saveRDS(All_Window_MassRep,paste("./cohen/4.Norquay_ReOrder_window_",ws, "_window_stats_MassRep.rds",sep=""))
#saveRDS(All_Window_MassRep_HMMStates,paste("./cohen/4.Norquay_ReOrder_window_",ws, "_window_stats_MassRep_HMMStates.rds",sep=""))
#saveRDS(All_Window_MassRep_HMMStates_PeakFilter,paste("./cohen/4.Norquay_ReOrder_window_",ws, "_window_stats_MassRep_HMMStates_PeakFilter.rds",sep=""))


# saveRDS(peakcluster_array_List,paste("./cohen/4.Norquay_ReOrder_window_",ws,"_window_stats_peakcluster_array_List.rds",sep=""))

#saveRDS(expec_peak_tallies_list,paste("./cohen/4.Norquay_ReOrder_window_",ws,"_window_stats_expec_peak_tallies_list.rds",sep=""))
#saveRDS(peak_tallies_list,paste("./cohen/4.Norquay_ReOrder_window_",ws,"_window_stats_peak_tallies_list.rds",sep=""))





#chromosome
andrew<-3
#state num plus one
tristan<-8
#seed
tate<-1

transformedpeaktallies<-log10(peak_tallies_list[[andrew]][[tristan]][[tate]])
transformedpeaktallies[transformedpeaktallies < -10] <- NA


plot(All_Window_MassRep[[andrew]][[tristan]][,3],All_Window_MassRep[[andrew]][[tristan]][,4],col=All_Window_MassRep[[andrew]][[tristan]][,6])


plot(seq(1,length(expec_peak_tallies_list[[andrew]][[tristan]][[tate]]),1),log10(expec_peak_tallies_list[[andrew]][[tristan]][[tate]]),
     ylim=c(min(log10(expec_peak_tallies_list[[andrew]][[tristan]][[tate]])),max(log10(expec_peak_tallies_list[[andrew]][[tristan]][[tate]]),transformedpeaktallies,na.rm=T)),
     ylab="log10(count)",xlab="Peak Cluster size",
     sub=paste("Chromosome ",andrew,". State numb ",(1+tristan),". Seed = ",(1000+tate),sep=""),
     main=paste("CohenW = ",opt_cW[[andrew]][tristan,tate],"\n pvalue = ",opt_chipvalue[[andrew]][tristan,tate],sep=""))
points(seq(1,length(expec_peak_tallies_list[[andrew]][[tristan]][[tate]]),1),transformedpeaktallies,col=2)


opt_chipvalue[[10]]
#peakcluster_array_List[[10]][[8]][[1]]


  