######################### Likelihood ratio test with normal distribution using lm with interactions but no relatedness matrix
library(lmtest)

my_data <- read.table("final_gea_input.txt", sep = "\t", h = T)
num_columns <- ncol(my_data)

# set up outputs
results <- data.frame()
coefficient_results <- data.frame()
models <- list()

for (i in 22:num_columns) {

 my_model <- lm(my_data[,i] ~ Topography  * H_Index * PC1, data = my_data)
 my_model_summary <- summary(my_model)

 my_null <- lm(my_data[,i] ~ H_Index * PC1, data = my_data)
 my_null_summary <- summary(my_null)

 # Significance based on likelihood ratio test
 ratio_test <- lrtest(my_null,my_model)
 p_value=ratio_test$`Pr(>Chisq)`[2]
 SNP <- colnames(my_data)[i]
 line <- cbind(SNP, p_value)
 results <- rbind(results, line)
 # Save all model coefficients, just in case
 models[[SNP]] <- my_model_summary$coefficients

 # Extracting p-values for Tpography coeffiecients and interactions involving topography and ntaking minimum
 p_values_coeff <- my_model_summary$coefficients[c(2:3, 6:9, 11:12), 4]
 line <- cbind(SNP, min(p_values_coeff))
 coefficient_results <- rbind(coefficient_results, line)

}

colnames(coefficient_results) <- c("SNP","coeff_p_value")


# Explore significant ie. column 25 SNP
library(ggplot2)
df <- my_data
# Create a histogram
ggplot(df, aes(x = df[, 4], fill = as.factor(df[, 25]))) +
  geom_bar(position = "dodge") +
  labs(title = "Distribution of Column 25 Across Topography",
       x = paste("Column", 25),
       y = "Frequency") +
  theme_minimal()








#### Same as above, but data.table
library(lmtest)
library(data.table)

# Load data using fread
my_data <- fread("final_gea_input.txt", sep = "\t", header = TRUE)

# Obtain the number of columns
num_columns <- ncol(my_data)

# set up outputs
results <- data.frame()
coefficient_results <- data.frame()
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
}

# Naming the columns of coefficient_results
colnames(coefficient_results) <- c("SNP", "coeff_p_value")



# Explore significant ie. column 25 SNP
library(ggplot2)
df <- my_data
# Create a histogram
ggplot(df, aes(x = df[[4]], fill = as.factor(df[[25]]))) +
  geom_bar(position = "dodge") +
  labs(title = "Distribution of Column 25 Across Topography",
       x = paste("Column", 25),
       y = "Frequency") +
  theme_minimal()









######################### LMEKIN (no interactions) but random effects relatedness matrix

library(coxme)

my_data <- read.table("final_gea_input.txt", sep = "\t", h = T)
num_columns <- ncol(my_data)
my_mat <- read.table("relatedness.txt")
my_mat <- as.matrix(my_mat)

# set up outputs
coefficient_results <- data.frame()
models <- list()

extract_coxme_table <- function (mod){
    beta <- mod$coefficients$fixed
    nvar <- length(beta)
    nfrail <- nrow(mod$var) - nvar
    se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
    z<- round(beta/se, 2)
    p<- signif(1 - pchisq((beta/se)^2, 1), 2)
    table=data.frame(cbind(beta,se,z,p))
    return(table)
}


for (i in 22:num_columns) {

SNP <- colnames(my_data)[i]
 my_model <- lmekin(my_data[,i] ~ Topography  + H_Index + PC1 + (1|Name), my_data, varlist = my_mat,na.action=na.omit)
 # Take topography min p-value
 my_table <- extract_coxme_table(my_model)
 my_table <- my_table[2:3,]
 p_values_coeff <- my_table$p
 line <- cbind(SNP, min(p_values_coeff))
 coefficient_results <- rbind(coefficient_results, line)
 # Save all models coefficients just in case
 models[[SNP]] <- extract_coxme_table(my_model)

}
colnames(coefficient_results) <- c("SNP","coeff_p_value")











######################### GLMMKIN  but random effects relatedness matrix   ### IN progress, where the p-values at

library(GMMAT)

my_data <- read.table("final_gea_input.txt", sep = "\t", h = T)
num_columns <- ncol(my_data)
my_mat <- read.table("relatedness.txt")
my_mat <- as.matrix(my_mat)

# set up outputs
coefficient_results <- data.frame()
models <- list()




for (i in 22:num_columns) {

SNP <- colnames(my_data)[i]
 my_model <- glmmkin(my_data[,i] ~ Topography  * H_Index * PC1, my_data,id = "Name", kin = my_mat,family = gaussian(link = "identity"))
 # Take topography min p-value


}
