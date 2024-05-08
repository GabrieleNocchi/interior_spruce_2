options(stringsAsFactors = F)
setwd("~/Dropbox/desktop/adaptree/snp_array")

frqs <- read.csv("Spr.Assn2.popfreq.new", header=T, sep="\t", na.strings=c("NA", "NaN", "-nan", ""),stringsAsFactors = FALSE)
#frqs <- read.csv("Spr.Assn2.popfrq", header=T, sep="\t", na.strings=c("NA", "NaN", "-nan", ""),stringsAsFactors = FALSE)

slpops <- read.csv("5_SeedlotObtained.csv", sep=",", header=T,stringsAsFactors = FALSE)
envs <- read.csv("spruce.pop.unscaled.envi", sep="\t", header=T,stringsAsFactors = FALSE)


envs$ID <- as.numeric (sapply (strsplit(envs$pop,split = "_"),"[[",2))
envs <- envs[,c("ID", colnames(envs)[c(24,2:23)])]



# final assembly ...
dat <- merge(frqs, slpops[ , c("ID", "SET")], by.x="pop", by.y="ID")


un1 <- sort (unique (dat$SET))
un2 <- c (un1[19:42],un1[46:52])
un3 <- c (un1[4:17],un1[19:42],un1[46:52])


dat2 <- dat[dat$SET %in% un3,]


dat2 <- merge(dat2, envs, by.x="pop", by.y="ID")

dat2 <- dat2[is.na(dat2$Frq) == F,]

#sprpops <- slpops [slpops$SPECIES == "SX",]

library ("lawstat")
library ("lme4")
library ("car")


the_pairs <- paste ("P",1:50,sep = "")
the_trans <- paste ("T",1:50,sep = "")

the_snps <- unique (dat2$SNP)

the_res <- array (NA, c (length (the_snps),9))

for (i in 1:length (the_snps)){
	
	
	sub1 <- dat2[which(dat2$SNP == the_snps[i]),]

	sub2 <- sub1[order(sub1$SET,sub1$ELEVATION),]

	#THE PAIRS
	sub3 <- sub2[sub2$SET %in% the_pairs,]
	
	low <- sub3[seq(1,(nrow(sub3)-1),by = 2),]
	high <- sub3[seq(2,nrow(sub3),by = 2),]

	mat1 <- array (NA, c(2,2,nrow(low)))
	
	mat1[1,1,] <- low$Frq * low$N
	mat1[1,2,] <- low$N - (low$Frq * low$N)

	mat1[2,1,] <- high$Frq * high$N
	mat1[2,2,] <- high$N - (high$Frq * high$N)
	
	cmh_res <- cmh.test(mat1)
	
	the_res [i,1] <- sub1[1,2]
	the_res [i,2] <- the_snps[i]
	the_res [i,4] <- cmh_res[[1]][3] #pvalue
	the_res [i,3] <- nrow(sub3)
	
	# par (mfcol = c (1,2))
	#these_trans <- unique (sub3$SET)
	# for (dd in 1:length (these_trans)){
		# sub4 <- sub3[sub3$SET == these_trans[dd],]
		# the_col <- array (NA,nrow(sub4))
		# the_col[sub4$species == "hybrid"] <- "purple"
		# the_col[sub4$species == "enge"] <- "red"
		# the_col[sub4$species == "glauca"] <- "blue"

		
		# if (dd == 1){
			# plot (sub4$ELEVATION,sub4$Frq, type = "l", pch = dd, xlim = range (sub3$ELEVATION),ylim = range (sub3$Frq),col = "black")
			# points (sub4$ELEVATION,sub4$Frq, type = "p", pch = dd, xlim = range (sub3$ELEVATION),ylim = range (sub3$Frq),col = the_col)

		# } else {
			# points (sub4$ELEVATION,sub4$Frq, type = "l", pch = dd,col = "black")	
			# points (sub4$ELEVATION,sub4$Frq, type = "p", pch = dd,col = the_col)	

		# }
	# }
		

	
	
	##the transects	
	sub3 <- sub2[sub2$SET %in% the_trans,]
	
	
	#mixed effects model:
	if (var(sub3$Frq) != 0){
		lm1 <- lmer (sub3$Frq ~ sub3$ELEVATION + (1 | sub3$SET))
		lm2 <- Anova(lm1)
	
		the_res[i,5] <- lm2[[3]][1]
	}
	#these_trans <- unique (sub3$SET)
	# for (dd in 1:length (these_trans)){
		# sub4 <- sub3[sub3$SET == these_trans[dd],]
		# the_col <- array (NA,nrow(sub4))
		# the_col[sub4$species == "hybrid"] <- "purple"
		# the_col[sub4$species == "enge"] <- "red"
		# the_col[sub4$species == "glauca"] <- "blue"

		
		# if (dd == 1){
			# plot (sub4$ELEVATION,sub4$Frq, type = "l", pch = dd, xlim = range (sub3$ELEVATION),ylim = range (sub3$Frq),col = "black")
			# points (sub4$ELEVATION,sub4$Frq, type = "p", pch = dd, xlim = range (sub3$ELEVATION),ylim = range (sub3$Frq),col = the_col)

		# } else {
			# points (sub4$ELEVATION,sub4$Frq, type = "l", pch = dd,col = "black")	
			# points (sub4$ELEVATION,sub4$Frq, type = "p", pch = dd,col = the_col)	

		# }
	# }
		
		
		
	##check ANOVA on hybrid ancestry:
	
	lmanc <- summary (lm (sub2$Frq ~ sub2$species))	
	lmanc2 <- anova (lm (sub2$Frq ~ sub2$species))	

	the_res[i,6] <- lmanc2[[5]][1]
	the_res[i,7] <- lmanc$r.squared
	
	
	lmMAT <- summary (lm (sub2$Frq ~ sub2$MAT))	
	lmMAT2 <- anova (lm (sub2$Frq ~ sub2$MAT))	

	the_res[i,8] <- lmMAT2[[5]][1]
	the_res[i,9] <- lmMAT$r.squared
				
				
		
}	


# # dd <- as.numeric (the_res[,4])
# ee <- as.numeric (the_res[,5])
# ind1 <- which (dd < 0.001 & ee < 0.001)

colnames(the_res) <- c ("contig","snp","num_pairs_pops","cmh_pvalue","transect_pvalue","hybrid_pvalue","hybrid_r2","MAT_pvalue","MAT_r2")

write.table (the_res,"results_SNP_array_pairs_transects_spruce_cmh.txt",row.names = F, col.names = T, quote = F)
