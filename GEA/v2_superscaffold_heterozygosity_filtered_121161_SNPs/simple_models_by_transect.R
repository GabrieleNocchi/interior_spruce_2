### NO PC1, but interactions
#### Same as above, but data.table

lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}


library(lmtest)
library(data.table)

# Load data using fread
my_data <- fread("final_gea_input.txt", sep = "\t", header = TRUE)

# Obtain the number of columns
num_columns <- ncol(my_data)

# Split data by Transect
transect_groups <- split(my_data, my_data$Transect)

# Loop through transect groups
for (transect_name in names(transect_groups)) {
  
  transect_data <- transect_groups[[transect_name]]
  
  results <- data.frame()
  coefficient_results <- data.frame()
  model_pvals <- data.frame()
  models <- list()
  
  # Loop through columns
  for (i in 22:num_columns) {
    
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
    
    # Extracting p-values for Topography coefficients and interactions involving topography and taking the minimum
    p_values_coeff <- my_model_summary$coefficients[c(2:3), 4]
    line <- cbind(SNP, min(p_values_coeff))
    coefficient_results <- rbind(coefficient_results, line)
    model_p <- lmp(my_model)
    line_2 <- cbind(SNP,model_p)
    model_pvals <- rbind(model_pvals, line_2)
  }
  
  # Naming the columns of coefficient_results
  colnames(coefficient_results) <- c("SNP", "coeff_p_value")
  colnames(model_pvals) <- c("SNP", "p_value")
  
  # Save models
  saveRDS(models, file = paste0("models_", transect_name, ".rds"))
  
  # Save results
  saveRDS(results, file = paste0("likelihood_ratio_results_", transect_name, ".rds"))
  
  # Save coefficient_results
  saveRDS(coefficient_results, file = paste0("coefficient_results_", transect_name, ".rds"))
  
  saveRDS(model_pvals, file = paste0("model_pval_results_", transect_name, ".rds"))
}













# Manhattan Plot by transect and LG
library(data.table)
library(gridExtra)
library(ggplot2)
library(dplyr)


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




for (lg_name in names(lg_groups)) {
  
  lg_data <- lg_groups[[lg_name]]
  lg_data$Transect <- factor(lg_data$Transect, levels = c("0. Overall", "5. Cascade", "1. Emerald", "6. Canmore", "2. Cathedral", "7. Yates", "3. Bell", "8. Old Baldy", "4. Norquay"))
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