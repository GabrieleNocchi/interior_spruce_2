library(lmtest)
library(data.table)

################# Adapt Tree #################

my_data <- fread("final_adapt_tree_banff_gea_input.txt", sep = "\t", header = TRUE)

# Obtain the number of columns
num_columns <- ncol(my_data)
my_data <- my_data[my_data$dataset == "adapt_tree",]

# set up outputs
results <- data.frame()
models <- list()

# Loop through columns
for (i in 8:num_columns) {
  
  # Fit the model
  my_model <- lm(my_data[[i]] ~ PC2, data = my_data)
  my_model_summary <- summary(my_model)
  
  my_null <- lm(my_data[[i]] ~ 1, data = my_data)
  my_null_summary <- summary(my_null)
  
  # Significance based on likelihood ratio test
  ratio_test <- lrtest(my_null, my_model)
  p_value <- ratio_test$`Pr(>Chisq)`[2]
  SNP <- colnames(my_data)[i]
  line <- cbind(SNP, p_value)
  results <- rbind(results, line)
  
  # Save all model coefficients, just in case
  models[[SNP]] <- my_model_summary$coefficients
 
}

# Save models
saveRDS(models, file = "adapt_tree_models_PC2.rds")
saveRDS(results, file = "adapt_tree_likelihood_ratio_results_PC2.rds")


################# Adapt Tree 2 #################

my_data <- fread("final_adapt_tree_banff_gea_input.txt", sep = "\t", header = TRUE)

# Obtain the number of columns
num_columns <- ncol(my_data)
my_data <- my_data[my_data$dataset == "adapt_tree",]

# set up outputs
results <- data.frame()
models <- list()

# Loop through columns
for (i in 8:num_columns) {
  
  # Fit the model
  my_model <- lm(my_data[[i]] ~ PC1, data = my_data)
  my_model_summary <- summary(my_model)
  
  my_null <- lm(my_data[[i]] ~ 1, data = my_data)
  my_null_summary <- summary(my_null)
  
  # Significance based on likelihood ratio test
  ratio_test <- lrtest(my_null, my_model)
  p_value <- ratio_test$`Pr(>Chisq)`[2]
  SNP <- colnames(my_data)[i]
  line <- cbind(SNP, p_value)
  results <- rbind(results, line)
  
  # Save all model coefficients, just in case
  models[[SNP]] <- my_model_summary$coefficients
 
}

# Save models
saveRDS(models, file = "adapt_tree_models_PC1.rds")
saveRDS(results, file = "adapt_tree_likelihood_ratio_results_PC1.rds")




################# Banff #################
my_data <- fread("final_adapt_tree_banff_gea_input.txt", sep = "\t", header = TRUE)

# Obtain the number of columns
num_columns <- ncol(my_data)
my_data <- my_data[my_data$dataset == "banff",]

# set up outputs
results <- data.frame()
models <- list()

# Loop through columns
for (i in 8:num_columns) {
  
  # Fit the model
  my_model <- lm(my_data[[i]] ~ Topography * H_Index, data = my_data)
  my_model_summary <- summary(my_model)
  
  my_null <- lm(my_data[[i]] ~ H_Index, data = my_data)
  my_null_summary <- summary(my_null)
  
  # Significance based on likelihood ratio test
  ratio_test <- lrtest(my_null, my_model)
  p_value <- ratio_test$`Pr(>Chisq)`[2]
  SNP <- colnames(my_data)[i]
  line <- cbind(SNP, p_value)
  results <- rbind(results, line)
  
  # Save all model coefficients, just in case
  models[[SNP]] <- my_model_summary$coefficients
  
}

# Save models
saveRDS(models, file = "banff_models.rds")
saveRDS(results, file = "banff_likelihood_ratio_results.rds")





################# Transects #################
my_data <- fread("final_adapt_tree_banff_gea_input.txt", sep = "\t", header = TRUE)

# Obtain the number of columns
num_columns <- ncol(my_data)
my_data <- my_data[my_data$dataset == "banff",]

# Split data by Transect
transect_groups <- split(my_data, my_data$Transect)

# Loop through transect groups
for (transect_name in names(transect_groups)) {
  
  transect_data <- transect_groups[[transect_name]]
  
  results <- data.frame()
  models <- list()
  
  # Loop through columns
  for (i in 8:num_columns) {
    
    # Fit the model
    my_model <- lm(transect_data[[i]] ~ Topography * H_Index, data = transect_data)
    my_model_summary <- summary(my_model)
    
    my_null <- lm(transect_data[[i]] ~ H_Index, data = transect_data)
    my_null_summary <- summary(my_null)
    
    # Significance based on likelihood ratio test
    ratio_test <- lrtest(my_null, my_model)
    p_value <- ratio_test$`Pr(>Chisq)`[2]
    SNP <- colnames(transect_data)[i]
    line <- cbind(SNP, p_value)
    results <- rbind(results, line)
    
    # Save all model coefficients, just in case
    models[[SNP]] <- my_model_summary$coefficients
    
  }
  
  # Save models
  saveRDS(models, file = paste0("models_", transect_name, ".rds"))
  saveRDS(results, file = paste0("likelihood_ratio_results_", transect_name, ".rds"))

}






######################### Plot









# Manhattan Plot by transect and LG
library(data.table)
library(gridExtra)
library(ggplot2)
library(dplyr)


# Format Pg- to Pg."
lg <- read.table("final_order.txt")


a <- readRDS("banff_likelihood_ratio_results.rds")
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




a <- readRDS("adapt_tree_likelihood_ratio_results_PC2.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("00. AdaptTree", length(a$SNP) )
a$index <- row.names(a)

merged_data_i <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_i) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")


merged_data <- rbind(merged_data_0, merged_data_a,merged_data_b,merged_data_c,merged_data_d,merged_data_e,merged_data_f,merged_data_g,merged_data_h, merged_data_i)



lg_groups <- split(merged_data, merged_data$LG)




for (lg_name in names(lg_groups)) {
  
  lg_data <- lg_groups[[lg_name]]
  lg_data$Transect <- factor(lg_data$Transect, levels = c("00. AdaptTree", "0. Overall", "1. Emerald", "5. Cascade","2. Cathedral",  "6. Canmore", "3. Bell",  "7. Yates", "4. Norquay", "8. Old Baldy" ))
# Load necessary libraries manhattan by LG
library(ggplot2)

# Create the Manhattan plot
manhattan_plot <- ggplot(lg_data, aes(x = as.numeric(index), y = -log10(p_value))) +
    geom_point(size = 1) +
    labs(x = "SNP Index", y = "-log10(p-value)") +
    theme_minimal() +
    facet_wrap(~ Transect, scales = "free", ncol = 2) + ggtitle(lg_name)

# Save the plot to a PDF file
pdf(paste(lg_name,".pdf"))
print(manhattan_plot)
dev.off()

}






# Just Overall vs AdaptTree

# Manhattan Plot by transect and LG
library(data.table)
library(gridExtra)
library(ggplot2)
library(dplyr)


# Format Pg- to Pg."
lg <- read.table("final_order.txt")


a <- readRDS("banff_likelihood_ratio_results.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("Banff", length(a$SNP) )
a$index <- row.names(a)

merged_data_00 <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_00) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")




a <- readRDS("adapt_tree_likelihood_ratio_results_PC2.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("AdaptTree PC2", length(a$SNP) )
a$index <- row.names(a)

merged_data_ii <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_ii) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")


a <- readRDS("adapt_tree_likelihood_ratio_results_PC1.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
a$site <- rep("AdaptTree PC1", length(a$SNP) )
a$index <- row.names(a)

merged_data_iii <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data_iii) <- c("SNP", "p_value","Transect","index","LG","cM","cDNA")


merged_data <- rbind(merged_data_00, merged_data_ii, merged_data_iii)



lg_groups <- split(merged_data, merged_data$LG)






for (lg_name in names(lg_groups)) {
  
  lg_data <- lg_groups[[lg_name]]
  lg_data$Transect <- factor(lg_data$Transect, levels = c("Banff","AdaptTree PC1", "AdaptTree PC2"))
# Load necessary libraries manhattan by LG
library(ggplot2)

# Create the Manhattan plot
manhattan_plot <- ggplot(lg_data, aes(x = as.numeric(index), y = -log10(p_value))) +
    geom_point(size = 1) +
    labs(x = "SNP Index", y = "-log10(p-value)") +
    theme_minimal() +
    facet_wrap(~ Transect, scales = "free", ncol = 1) + ggtitle(lg_name)

# Save the plot to a PDF file
pdf(paste(lg_name,"_banff_adapt_tree.pdf"))
print(manhattan_plot)
dev.off()

}

