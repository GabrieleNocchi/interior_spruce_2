setwd("~/Dropbox/desktop/adaptree/snp_array")

blast1 <- read.table ("~/Dropbox/desktop/adaptree/snp_array/axiom_spruce_SNP_array_to_WS77111-v3.txt")

blast1 <- blast1[order(blast1[,1],blast1[,12], decreasing = T),]

map_one <- read.table ("~/Dropbox/desktop/adaptree/snp_array/gnavigator_full-GM-stats_summary.tsv", sep = "\t")

map_two <- read.table ("~/Dropbox/desktop/adaptree/snp_array/gnavigator_full-cDNA-stats_summary.tsv", sep = "\t")

map_three <- read.table ("~/Dropbox/desktop/adaptree/snp_array/Pglauca_geneticMap_v2.GCATfilter.tsv", sep = "\t")

the_queries <- unique (blast1[,1])

to_keep <- NULL

#keep top hit for each query (after ordering above)
for (i in 1:length (the_queries)){
	
	ind1 <- which (blast1[,1] == the_queries[i])
	
	to_keep <- c(to_keep,ind1[1])
	
}

blast2 <- blast1[to_keep,]

write.table (blast2,"axiom_spruce_SNP_array_to_WS77111-v3_top_hits.txt")

merg1 <- merge (map_two,map_three,by.x = "V1", by.y = "V3", all.x = T)

#combine info on blast from SNP array to white spruce v3 genome
merg2 <- merge (blast2,merg1,by.x = "V2",by.y = "V3",all.x = T)

#pick out hits from merg2 with contig__pos vs. SNP ID from elsewhere
merg2b <- merg2[grep ("__",merg2[,2]),]
merg2c <- merg2[grep ("__",merg2[,2],invert = T),]


split1 <- strsplit (merg2b[,2],split = "__")
chrom <- sapply (split1,"[[",1)
pos <- sapply (split1,"[[",2)

map_2b <- cbind(chrom,pos,merg2b$V4.y,merg2b$V2.y)


split1 <- strsplit (merg2c[,2],split = "__")
chrom <- sapply (split1,"[[",1)
pos <- array (NA, length (chrom))

map_2c <- cbind(chrom,pos,merg2c$V4.y,merg2c$V2.y)


out_map <- rbind (map_2b,map_2c)

write.table (out_map,"temp555.txt")
out_map <- read.table ("temp555.txt")

out_map <- out_map[order(out_map[,3],as.numeric (out_map[,4])),] #includes lots of repeated stuff

colnames(out_map)[3:4] <- c ("LG","cM")

out_map <- out_map[out_map[,3] != "unknownLG",]

out_map$ID <- paste (out_map$chrom,out_map$pos,sep = "__")


uniqID <- unique (out_map$ID)

new_out <- NULL

for (i in 1:length (uniqID)){

	ind1 <- which (out_map$ID == uniqID[i])
	
	if (length (ind1) > 1){
		if (length (unique (out_map[ind1,3])) == 1){ #if the SNP is only mapping to one linkage map
		
			this_out <- out_map[ind1[1],]
			this_out[4] <- median (out_map[ind1,4]) #take the median cM across all matches
			
		}
	} else {
		
		this_out <- out_map[ind1,]	
		new_out <- rbind (new_out, this_out)
	}
}

write.table (new_out,"mapping_to_linkage_map.txt",row.names = F, col.names = F, quote = F, sep = "\t")
#read in results from CMH:

library ("lawstat")
library ("lme4")
library ("car")

library (rgl)
sqrt(2 * ((r1 * r2) / (r1 + r2))^2)
rad4 <- function(r1,r2){
	sqrt(2 * ((r1 * r2) / (r1 + r2))^2)	
}

#results made by "analyze_spruce_freqs_bySNP.R"
the_res <- read.table ("results_SNP_array_pairs_transects_spruce_cmh.txt",T)

colnames(the_res) <- c("contig","SNP","N","pairs","transects","hybrid_p","hybrid_R","mat_p","mat_R")




the_res$qnorm1 <- qnorm(the_res[,4])
the_res$qnorm2 <- qnorm(the_res[,5])
temp11 <- (the_res$qnorm1 + the_res$qnorm2) / sqrt(2)
the_res$elevation <- pnorm(temp11)
the_res$combined <- rad4((-1*log10(the_res[,4])),(-1*log10(the_res[,5])))

the_res[,2] <- gsub (".","__",the_res[,2],fixed = T)

res1 <- merge (new_out,the_res,by.x = "ID",by.y = "SNP",all.x = T)

write.table (res1,"temp556.txt")
res1 <- read.table ("temp556.txt")

res2 <- res1[order (res1$LG,as.numeric (res1$cM)),]

res2$total_pos <- NA
res2$colour1 <- 1

maxpos <- tapply (res2$cM,res2$LG,max)
maxpos2 <- maxpos[1:11] 
names(maxpos2)[1:11] <- names(maxpos)[2:12]

maxpos2_names <- names(maxpos2)

for (i in 1:nrow(res2)){
	if (res2[i,4] != "LG01"){
		ind4 <- which (maxpos2_names == res2[i,4])
		sumlength <- sum(maxpos2[1:ind4])
		res2$total_pos[i] <- res2$cM[i] + sumlength
		res2$colour1[i] <- (ind4+1)
	}
}

plot (res2$total_pos,-1*log10(res2$mat_p), col = res2$colour1)



plot (res2$total_pos,-1*log10(res2$hybrid_p), col = res2$colour1)


par (mfcol = c (2,1))
plot (res2$total_pos,-1*log10(res2$elevation), col = res2$colour1)


plot (res2$total_pos,(res2$combined), col = res2$colour1)

plot (res2$total_pos,-1*log10(res2$mat_p), col = res2$colour1)


