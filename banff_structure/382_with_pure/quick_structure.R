library(pophelper)
a <- readQ("382_spruce_filtered_MAF_pure_fast_output.2.meanQ")
n <- read.table("names.txt")
row.names(a[[1]]) <- n$V1
# plotQMultiline(a, exportpath = getwd(), sortind = "Cluster2", useindlab = T,  indlabsize = 3, clustercol = c("orange","grey"))
#plotQ(a, exportpath = getwd(), sortind = "Cluster1", useindlab = T,  indlabsize = 3, clustercol = c("grey","orange"))


b <- read.table("joined_metadata_2.txt", sep = "\t")
b <- cbind(b$V3,b$V4)
b <- as.data.frame(b)
# plotQMultiline(a, exportpath = getwd(), sortind = "Cluster2", useindlab = T,  indlabsize = 3, clustercol = c("orange","grey"), grplab = b, ordergrp = T)




b <- read.table("joined_metadata.txt", sep = "\t")
b <- b$V4
b <- as.data.frame(b)
# plotQMultiline(a, exportpath = getwd(), sortind = "Cluster2", useindlab = T,  indlabsize = 3, clustercol = c("orange","grey"), grplab = b, ordergrp = T)





# Going Crazy!

b <- read.table("joined_metadata_2.txt", sep = "\t")
b <- cbind(b$V3,b$V4)
b <- as.data.frame(b)



matching_rows_b <- which(b$V1 == "1.Emerald")

# Step 2: Extract the corresponding rows in 'a'
matching_rows_a_1 <- a[[1]][matching_rows_b, , drop = FALSE]
matching_rows_a_1 <- list(matching_rows_a_1)
names(matching_rows_a_1) <- "Emerald"
matching_b_1 <- b[b$V1=="1.Emerald",]
plotQMultiline(matching_rows_a_1, exportpath = getwd(), sortind = "Cluster2", showindlab = FALSE, clustercol = c("orange","grey"), imgtype='pdf', grplabsize = 3, grplab = matching_b_1, ordergrp = T, height = 2, width = 10)



matching_rows_b <- which(b$V1 == "2.Cathedral")

# Step 2: Extract the corresponding rows in 'a'
matching_rows_a_1 <- a[[1]][matching_rows_b, , drop = FALSE]
matching_rows_a_1 <- list(matching_rows_a_1)
names(matching_rows_a_1) <- "Cathedral"
matching_b_1 <- b[b$V1=="2.Cathedral",]
plotQMultiline(matching_rows_a_1, exportpath = getwd(), sortind = "Cluster2", showindlab = FALSE, clustercol = c("orange","grey"),imgtype='pdf', grplabsize = 3, grplab = matching_b_1, ordergrp = T, height = 2, width = 10)



matching_rows_b <- which(b$V1 == "3.Bell")

# Step 2: Extract the corresponding rows in 'a'
matching_rows_a_1 <- a[[1]][matching_rows_b, , drop = FALSE]
matching_rows_a_1 <- list(matching_rows_a_1)
names(matching_rows_a_1) <- "Bell"
matching_b_1 <- b[b$V1=="3.Bell",]
plotQMultiline(matching_rows_a_1, exportpath = getwd(), sortind = "Cluster2", showindlab = FALSE, clustercol = c("orange","grey"),imgtype='pdf', grplabsize = 3, grplab = matching_b_1, ordergrp = T, height = 2, width = 10)




matching_rows_b <- which(b$V1 == "4.Norquay")

# Step 2: Extract the corresponding rows in 'a'
matching_rows_a_1 <- a[[1]][matching_rows_b, , drop = FALSE]
matching_rows_a_1 <- list(matching_rows_a_1)
names(matching_rows_a_1) <- "Norquay"
matching_b_1 <- b[b$V1=="4.Norquay",]
plotQMultiline(matching_rows_a_1, exportpath = getwd(), sortind = "Cluster2", showindlab = FALSE, clustercol = c("orange","grey"), imgtype='pdf',grplabsize = 3,grplab = matching_b_1, ordergrp = T, height = 2, width = 10)



matching_rows_b <- which(b$V1 == "5.Cascade")

# Step 2: Extract the corresponding rows in 'a'
matching_rows_a_1 <- a[[1]][matching_rows_b, , drop = FALSE]
matching_rows_a_1 <- list(matching_rows_a_1)
names(matching_rows_a_1) <- "Cascade"
matching_b_1 <- b[b$V1=="5.Cascade",]
plotQMultiline(matching_rows_a_1, exportpath = getwd(), sortind = "Cluster2", showindlab = FALSE, clustercol = c("orange","grey"),imgtype='pdf',grplabsize = 3, grplab = matching_b_1, ordergrp = T, height = 2, width = 10)



matching_rows_b <- which(b$V1 == "6.Canmore")

# Step 2: Extract the corresponding rows in 'a'
matching_rows_a_1 <- a[[1]][matching_rows_b, , drop = FALSE]
matching_rows_a_1 <- list(matching_rows_a_1)
names(matching_rows_a_1) <- "Canmore"
matching_b_1 <- b[b$V1=="6.Canmore",]
plotQMultiline(matching_rows_a_1, exportpath = getwd(), sortind = "Cluster2", showindlab = FALSE, clustercol = c("orange","grey"), grplabsize = 3,imgtype='pdf', grplab = matching_b_1, ordergrp = T, height = 2, width = 10)



matching_rows_b <- which(b$V1 == "7.Yates")

# Step 2: Extract the corresponding rows in 'a'
matching_rows_a_1 <- a[[1]][matching_rows_b, , drop = FALSE]
matching_rows_a_1 <- list(matching_rows_a_1)
names(matching_rows_a_1) <- "Yates"
matching_b_1 <- b[b$V1=="7.Yates",]
plotQMultiline(matching_rows_a_1, exportpath = getwd(), sortind = "Cluster2", showindlab = FALSE, clustercol = c("orange","grey"), grplabsize = 3,imgtype='pdf',grplab = matching_b_1, ordergrp = T, height = 2, width = 10)



matching_rows_b <- which(b$V1 == "8.Old Baldy")

# Step 2: Extract the corresponding rows in 'a'
matching_rows_a_1 <- a[[1]][matching_rows_b, , drop = FALSE]
matching_rows_a_1 <- list(matching_rows_a_1)
names(matching_rows_a_1) <- "Old Baldy"
matching_b_1 <- b[b$V1=="8.Old Baldy",]
plotQMultiline(matching_rows_a_1, exportpath = getwd(), sortind = "Cluster2", showindlab = FALSE, clustercol = c("orange","grey"), grplabsize = 3,imgtype='pdf', grplab = matching_b_1, ordergrp = T, height = 2, width = 10)
