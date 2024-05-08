library(data.table)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(dgconstraint)


# Format Pg- to Pg."
lg <- read.table("final_order.txt")


a <- readRDS("likelihood_ratio_results.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("0. Overall", length(a$SNP) )
a$index <- row.names(a)

merged_data_0 <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_0) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")



a <- readRDS("likelihood_ratio_results_1.Emerald.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("1. Emerald", length(a$SNP) )
a$index <- row.names(a)

merged_data_a <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_a) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")


a <- readRDS("likelihood_ratio_results_2.Cathedral.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("2. Cathedral", length(a$SNP) )
a$index <- row.names(a)

merged_data_b <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_b) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")

a <- readRDS("likelihood_ratio_results_3.Bell.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("3. Bell", length(a$SNP) )
a$index <- row.names(a)

merged_data_c <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_c) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")

a <- readRDS("likelihood_ratio_results_4.Norquay.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("4. Norquay", length(a$SNP) )
a$index <- row.names(a)

merged_data_d <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_d) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")

a <- readRDS("likelihood_ratio_results_5.Cascade.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("5. Cascade", length(a$SNP) )
a$index <- row.names(a)

merged_data_e <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_e) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")

a <- readRDS("likelihood_ratio_results_6.Canmore.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("6. Canmore", length(a$SNP) )
a$index <- row.names(a)

merged_data_f <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_f) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")

a <- readRDS("likelihood_ratio_results_7.Yates.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("7. Yates", length(a$SNP) )
a$index <- row.names(a)

merged_data_g <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_g) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")

a <- readRDS("likelihood_ratio_results_8.Old Baldy.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("8. Old Baldy", length(a$SNP) )
a$index <- row.names(a)

merged_data_h <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_h) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")


merged_data <- rbind(merged_data_0, merged_data_a,merged_data_b,merged_data_c,merged_data_d,merged_data_e,merged_data_f,merged_data_g,merged_data_h)

lg_groups <- split(merged_data, merged_data$LG)


trans_list <- list()
all_list <- list()

for (lg_name in names(lg_groups)) {
  
  lg_data <- lg_groups[[lg_name]]
  


# Create empty dataframe to store results
results_df <- data.frame(
  Vector1 = character(),
  Vector2 = character(),
  Score = numeric(),
  stringsAsFactors = FALSE
)

# Function to calculate score
calculate_score <- function(vector1, vector2) {
  # Assuming single_c_hyper is a function that calculates the score
  single_c_hyper(vector1, vector2, na.rm = FALSE)
}


temp_data <- lg_data[lg_data$Transect == "1. Emerald",]
emerald <- as.vector(ifelse(is.na(temp_data$p_value) | temp_data$p_value >= 0.0001, 0, 1))


temp_data  <- as.vector(lg_data[lg_data$Transect == "2. Cathedral",])
cathedral <- as.vector(ifelse(is.na(temp_data$p_value) | temp_data$p_value >= 0.0001, 0, 1))

temp_data  <- as.vector(lg_data[lg_data$Transect == "3. Bell",])
bell <- as.vector(ifelse(is.na(temp_data$p_value) | temp_data$p_value >= 0.0001, 0, 1))

temp_data  <- as.vector(lg_data[lg_data$Transect == "4. Norquay",])
norquay <- as.vector(ifelse(is.na(temp_data$p_value) | temp_data$p_value >= 0.0001, 0, 1))

temp_data  <- as.vector(lg_data[lg_data$Transect == "5. Cascade",])
cascade <- as.vector(ifelse(is.na(temp_data$p_value) | temp_data$p_value >= 0.0001, 0, 1))

temp_data  <- as.vector(lg_data[lg_data$Transect == "6. Canmore",])
canmore <- as.vector(ifelse(is.na(temp_data$p_value) | temp_data$p_value >= 0.0001, 0, 1))

temp_data  <- as.vector(lg_data[lg_data$Transect == "7. Yates",])
yates <- as.vector(ifelse(is.na(temp_data$p_value) | temp_data$p_value >= 0.0001, 0, 1))

temp_data  <- as.vector(lg_data[lg_data$Transect == "8. Old Baldy",])
old <- as.vector(ifelse(is.na(temp_data$p_value) | temp_data$p_value >= 0.0001, 0, 1))

temp_data  <- as.vector(lg_data[lg_data$Transect == "0. Overall",])
all <- as.vector(ifelse(is.na(temp_data$p_value) | temp_data$p_value >= 0.0000000001, 0, 1))


# List of vectors
vectors <- list(
  emerald = emerald,
  cathedral = cathedral,
  bell = bell,
  norquay = norquay,
  cascade = cascade,
  canmore = canmore,
  yates = yates,
  old = old,
  all = all
) 
 
# Iterate over all possible pairs of vectors
for (i in 1:(length(vectors) - 1)) {
  for (j in (i + 1):length(vectors)) {
    vector1_name <- names(vectors)[i]
    vector2_name <- names(vectors)[j]
    
    score <- calculate_score(vectors[[i]], vectors[[j]])
    
    # Store results in dataframe
    results_df <- rbind(results_df, data.frame(
      Vector1 = vector1_name,
      Vector2 = vector2_name,
      Score = score
    ))
  }
}

results_df <- results_df[order(results_df$Score, decreasing = TRUE), ]


all_results <- results_df[results_df$Vector2 == "all",]
transects_results <- results_df[results_df$Vector2 != "all",]


# Display the results
trans_list[[lg_name]] <- transects_results
all_list[[lg_name]] <- all_results
 
  }