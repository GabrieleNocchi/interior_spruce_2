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


#sprpops <- slpops [slpops$SPECIES == "SX",]

library ("lawstat")
library ("lme4")
library ("car")

library (rgl)

##sqrt(2 * ((r1 * r2) / (r1 + r2))^2)
rad4 <- function(r1,r2){
	sqrt(2 * ((r1 * r2) / (r1 + r2))^2)	
}

the_pairs <- paste ("P",1:50,sep = "")
the_trans <- paste ("T",1:50,sep = "")

the_snps <- unique (dat2$SNP)

#results made by "analyze_spruce_freqs_bySNP.R"
the_res <- read.table ("results_SNP_array_pairs_transects_spruce_cmh.txt",T)

colnames(the_res) <- c("contig","SNP","N","pairs","transects","hybrid_p","hybrid_R","mat_p","mat_R")
#"pairs" is the p-value from a cmh test of the pairs of high/low populations
#transects is the p-value of a mixed effect model regressing each population's frequency on elevation and with transect as a random effect
#hybrid_p and hybrid_R are p and R^2 values for freq ~ species is a factor with three levels (white, engelmann, glauca)
#mat_p and mat_R are p and R^2 values for freq ~ MAT (mean annual temperature)



the_res$qnorm1 <- qnorm(the_res[,4])
the_res$qnorm2 <- qnorm(the_res[,5])
temp11 <- (the_res$qnorm1 + the_res$qnorm2) / sqrt(2)
the_res$elevation <- pnorm(temp11)

the_res$combined <- rad4((-1*log10(the_res[,4])),(-1*log10(the_res[,5])))

toplot <- the_res[which(the_res$elevation < 0.01 & the_res[,4] != 0),]
toplot <- the_res[which(the_res[,4] != 0),]

plot3d (log10(toplot$elevation),log10(toplot[,6]),log10(toplot[,8]),xlab = "elevation",ylab = "hybrid", zlab = "MAT")

plot3d (log10(toplot$elevation),log10(toplot$pairs),log10(toplot$transects),xlab = "elevation",ylab = "pairs", zlab = "transects")

plot3d (envs$ELEVATION, envs$LAT,envs$MAT)

plot3d (log10(toplot[,4]),log10(toplot[,6]),log10(toplot[,8]),xlab = "elevation",ylab = "hybrid", zlab = "MAT")

##pairs, transects, hybrid, elevation
par(mfcol = c (2,3))
plot (log10(toplot$hybrid_p),log10(toplot$elevation),xlab = "hybrid p",ylab = "combined p pairs + transects")
plot (log10(toplot$hybrid_p),log10(toplot$pairs),xlab = "hybrid p",ylab = "pairs p")
plot (log10(toplot$hybrid_p),log10(toplot$transects),xlab = "hybrid p",ylab = "transects p")
plot (log10(toplot$pairs),log10(toplot$transects),xlab = "pairs",ylab = "transects p")
plot (log10(toplot$hybrid_p),log10(toplot$mat_p),xlab = "hybrid p",ylab = "MAT p")
plot (log10(toplot$mat_p),log10(toplot$elevation),xlab = "MAT p",ylab = "combined p pairs + transects")


dd <- as.numeric (the_res[,4])
ee <- as.numeric (the_res[,5])
ind1 <- which (dd < 0.001 & ee < 0.001)


pdf ("plotter_spruce_snps.pdf")

for (i in 1:length (ind1)){
	
	
	sub1 <- dat2[which(dat2$SNP == the_snps[ind1[i]]),]

	sub2 <- sub1[order(sub1$SET,sub1$ELEVATION),]

	#THE PAIRS
	sub3 <- sub2[sub2$SET %in% the_pairs,]
	
	low <- sub3[seq(1,(nrow(sub3)-1),by = 2),]
	high <- sub3[seq(2,nrow(sub3),by = 2),]

	# mat1 <- array (NA, c(2,2,nrow(low)))
	
	# mat1[1,1,] <- low$Frq * low$N
	# mat1[1,2,] <- low$N - (low$Frq * low$N)

	# mat1[2,1,] <- high$Frq * high$N
	# mat1[2,2,] <- high$N - (high$Frq * high$N)
	
	# cmh_res <- cmh.test(mat1)
	
	# the_res [i,1] <- sub1[1,2]
	# the_res [i,2] <- the_snps[i]
	# the_res [i,4] <- cmh_res[[1]][3] #pvalue
	# the_res [i,3] <- nrow(sub3)
	
	par (mfcol = c (2,1))
	these_trans <- unique (sub3$SET)
	for (dd in 1:length (these_trans)){
		sub4 <- sub3[sub3$SET == these_trans[dd],]
		the_col <- array (NA,nrow(sub4))
		the_col[sub4$species == "hybrid"] <- "purple"
		the_col[sub4$species == "enge"] <- "red"
		the_col[sub4$species == "glauca"] <- "blue"

		
		if (dd == 1){
			plot (sub4$ELEVATION,sub4$Frq, type = "l", pch = dd, xlim = range (sub3$ELEVATION),ylim = range (sub3$Frq),col = "black")
			points (sub4$ELEVATION,sub4$Frq, type = "p", pch = dd, xlim = range (sub3$ELEVATION),ylim = range (sub3$Frq),col = the_col)

		} else {
			points (sub4$ELEVATION,sub4$Frq, type = "l", pch = dd,col = "black")	
			points (sub4$ELEVATION,sub4$Frq, type = "p", pch = dd,col = the_col)	

		}
	}
		

	
	
	##the transects	
	sub3 <- sub2[sub2$SET %in% the_trans,]
	
	
	#mixed effects model:
	if (var(sub3$Frq) != 0){
		lm1 <- lmer (sub3$Frq ~ sub3$ELEVATION + (1 | sub3$SET))
		lm2 <- Anova(lm1)
	
		the_res[i,5] <- lm2[[3]][1]
	}
	these_trans <- unique (sub3$SET)
	for (dd in 1:length (these_trans)){
		sub4 <- sub3[sub3$SET == these_trans[dd],]
		the_col <- array (NA,nrow(sub4))
		the_col[sub4$species == "hybrid"] <- "purple"
		the_col[sub4$species == "enge"] <- "red"
		the_col[sub4$species == "glauca"] <- "blue"

		
		if (dd == 1){
			plot (sub4$ELEVATION,sub4$Frq, type = "l", pch = dd, xlim = range (sub3$ELEVATION),ylim = range (sub3$Frq),col = "black")
			points (sub4$ELEVATION,sub4$Frq, type = "p", pch = dd, xlim = range (sub3$ELEVATION),ylim = range (sub3$Frq),col = the_col)

		} else {
			points (sub4$ELEVATION,sub4$Frq, type = "l", pch = dd,col = "black")	
			points (sub4$ELEVATION,sub4$Frq, type = "p", pch = dd,col = the_col)	

		}
	}
		
		
		
}	

dev.off()

