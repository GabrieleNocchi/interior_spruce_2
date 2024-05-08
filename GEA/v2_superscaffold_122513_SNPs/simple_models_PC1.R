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

# set up outputs
results <- data.frame()
coefficient_results <- data.frame()
model_pvals <- data.frame()
models <- list()

# Loop through columns
for (i in 22:num_columns) {
  
  # Fit the model
  my_model <- lm(my_data[[i]] ~ Topography * H_Index * PC1, data = my_data)
  my_model_summary <- summary(my_model)
  
  my_null <- lm(my_data[[i]] ~ H_Index * PC1, data = my_data)
  my_null_summary <- summary(my_null)
  
  # Significance based on likelihood ratio test
  ratio_test <- lrtest(my_null, my_model)
  p_value <- ratio_test$`Pr(>Chisq)`[2]
  SNP <- colnames(my_data)[i]
  line <- cbind(SNP, p_value)
  results <- rbind(results, line)
  
  # Save all model coefficients, just in case
  models[[SNP]] <- my_model_summary$coefficients
  
  # Extracting p-values for Topography coefficients and interactions involving topography and taking the minimum
  p_values_coeff <- my_model_summary$coefficients[c(2:3, 6:9, 11:12), 4]
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
saveRDS(models, file = "models_PC1.rds")

# Save results
saveRDS(results, file = "likelihood_ratio_results_PC1.rds")

# Save coefficient_results
saveRDS(coefficient_results, file = "coefficient_results_PC1.rds")

saveRDS(model_pvals, file = "model_pval_results_PC1.rds")

