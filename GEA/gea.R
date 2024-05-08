my_data <- read.table("final_gea_input.txt", sep = "\t", h = T)


# Choose distribution based on smaller AIC
library(MASS)
df <- data.frame()
num_columns <- ncol(my_data)


for (i in 22:num_columns) {

 gau <- AIC(fitdistr(as.numeric(na.omit(my_data[,i])), "normal"))
 poi <- AIC(fitdistr(as.numeric(na.omit(my_data[,i])), "Poisson"))
 nbd <- AIC(fitdistr(as.numeric(na.omit(my_data[,i])), "negative binomial"))
 line <- cbind(gau,nbd,poi)
 df <- rbind(df,line)

}


# Checking correlation
library(ggcorrplot)
library(dplyr)
library(ggplot2)

a <- my_data[,4:5]
b <- my_data[,20:21]
z <- cbind(a,b)

model.matrix(~0+., data=z) %>%
cor(use="pairwise.complete.obs") %>%
ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2)





## Different method to test correlation more appropriate for when there are both numerical and categorical predictors
library(tidyverse)
library(rcompanion)


# Calculate a pairwise association between all variables in a data-frame. In particular nominal vs nominal with Chi-square, numeric vs numeric with Pearson correlation, and nominal vs numeric with ANOVA.
# Adopted from https://stackoverflow.com/a/52557631/590437
mixed_assoc = function(df, cor_method="spearman", adjust_cramersv_bias=TRUE){
 df_comb = expand.grid(names(df), names(df),  stringsAsFactors = F) %>% set_names("X1", "X2")

 is_nominal = function(x) class(x) %in% c("factor", "character")
 # https://community.rstudio.com/t/why-is-purr-is-numeric-deprecated/3559
 # https://github.com/r-lib/rlang/issues/781
 is_numeric <- function(x) { is.integer(x) || is_double(x)}

 f = function(xName,yName) {
  x =  pull(df, xName)
  y =  pull(df, yName)

  result = if(is_nominal(x) && is_nominal(y)){
   # use bias corrected cramersV as described in https://rdrr.io/cran/rcompanion/man/cramerV.html
   cv = cramerV(as.character(x), as.character(y), bias.correct = adjust_cramersv_bias)
   data.frame(xName, yName, assoc=cv, type="cramersV")

  }else if(is_numeric(x) && is_numeric(y)){
   correlation = cor(x, y, method=cor_method, use="complete.obs")
   data.frame(xName, yName, assoc=correlation, type="correlation")

  }else if(is_numeric(x) && is_nominal(y)){
   # from https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable/124618#124618
   r_squared = summary(lm(x ~ y))$r.squared
   data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")

  }else if(is_nominal(x) && is_numeric(y)){
   r_squared = summary(lm(y ~x))$r.squared
   data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")

  }else {
   warning(paste("unmatched column type combination: ", class(x), class(y)))
  }

  # finally add complete obs number and ratio to table
  result %>% mutate(complete_obs_pairs=sum(!is.na(x) & !is.na(y)), complete_obs_ratio=complete_obs_pairs/length(x)) %>% rename(x=xName, y=yName)
 }

 # apply function to each variable combination
 map2_df(df_comb$X1, df_comb$X2, f)
}

mixed_assoc(z)


require(corrr)
library(viridis)
z %>%
mixed_assoc() %>%
select(x, y, assoc) %>%
spread(y, assoc) %>%
column_to_rownames("x") %>%
as.matrix %>%
as_cordf %>%
network_plot() +
scale_color_viridis_c(limits = c(-1, 1))





##### MODELS TESTS #####


######################### Null vs Residual deviance rather than LRT to calculate model significance
num_columns <- ncol(my_data)
results <- data.frame()

for (i in 22:num_columns) {
 my_model <- glm(my_data[,i] ~ Topography  * H_Index * PC1, data = my_data, family = "gaussian")
 my_model_summary <- summary(my_model)

 # Significance based on chi square between null and residual deviance
 p_value <- 1 - pchisq(my_model_summary$null.deviance - my_model_summary$deviance, my_model_summary$df.null - my_model_summary$df.residual)
 SNP <- colnames(my_data)[i]
 line <- cbind(SNP, p_value)
 results <- rbind(results, line)
}





######################### Normal Distribution lm with builtin F statistic p-value and saving coefficients
num_columns <- ncol(my_data)
results <- data.frame()
models <- list()

for (i in 22:num_columns) {
 my_model <- lm(my_data[,i] ~ Topography  * H_Index * PC1, data = my_data)
 my_model_summary <- summary(my_model)
 f <- my_model_summary$fstatistic
 p_value <- pf(f[1],f[2],f[3],lower.tail=F)
 SNP <- colnames(my_data)[i]
 line <- cbind(SNP, p_value)
 results <- rbind(results, line)
 models[[SNP]] <- my_model_summary$coefficients
}





######################### Likelihood ratio test with normal distribution and saving coefficients
library(lmtest)
num_columns <- ncol(my_data)
results <- data.frame()
models <- list()

for (i in 22:num_columns) {
 my_model <- glm(my_data[,i] ~ Topography  * H_Index * PC1, data = my_data, family = "gaussian")
 my_model_summary <- summary(my_model)

 my_null <- glm(my_data[,i] ~ H_Index * PC1, data = my_data, family = "gaussian")
 my_null_summary <- summary(my_null)

 # Significance based on likelihood ratio test
 ratio_test <- lrtest(my_null,my_model)
 p_value=ratio_test$`Pr(>Chisq)`[2]
 SNP <- colnames(my_data)[i]
 line <- cbind(SNP, p_value)
 results <- rbind(results, line)
 models[[SNP]] <- my_model_summary$coefficients
}





######################### Likelihood ratio test with binomial distribution and saving coefficients
######################### Do not use unless you have 0 and 1 - do not attempt converting 0,1,2 to 0, 0.5 and 1 -- it works, but wrongly
library(lmtest)
num_columns <- ncol(my_data)
results <- data.frame()
models <- list()

for (i in 22:num_columns) {

 my_model <- glm(my_data[,i]~ Topography  * H_Index * PC1, data = my_data, family = "binomial")
 my_model_summary <- summary(my_model)

 my_null <- glm(my_data[,i] ~ H_Index * PC1, data = my_data, family = "binomial")
 my_null_summary <- summary(my_null)

 # Significance based on likelihood ratio test
 ratio_test <- lrtest(my_null,my_model)
 p_value=ratio_test$`Pr(>Chisq)`[2]
 SNP <- colnames(my_data)[i]
 line <- cbind(SNP, p_value)
 results <- rbind(results, line)
 models[[SNP]] <- my_model_summary$coefficients
}





######################### Likelihood ratio test with Poisson distribution and saving coefficients
library(lmtest)
num_columns <- ncol(my_data)
results <- data.frame()
models <- list()

for (i in 22:num_columns) {

 my_model <- glm(my_data[,i] ~ Topography  * H_Index * PC1, data = my_data, family = "poisson")
 my_model_summary <- summary(my_model)

 my_null <- glm(my_data[,i] ~ H_Index * PC1, data = my_data, family = "poisson")
 my_null_summary <- summary(my_null)

 # Significance based on likelihood ratio test
 ratio_test <- lrtest(my_null,my_model)
 p_value=ratio_test$`Pr(>Chisq)`[2]
 SNP <- colnames(my_data)[i]
 line <- cbind(SNP, p_value)
 results <- rbind(results, line)
 models[[SNP]] <- my_model_summary$coefficients

}




######################### Likelihood ratio test using negative binomial distribution and saving coefficients
library(lmtest)
library(MASS)

num_columns <- ncol(my_data)
results <- data.frame()
models <- list()

for (i in 22:num_columns) {

 my_model <- glm.nb(my_data[,i] ~ Topography  * H_Index * PC1, data = my_data)
 my_model_summary <- summary(my_model)

 my_null <- glm.nb(my_data[,i] ~ H_Index * PC1, data = my_data)
 my_null_summary <- summary(my_null)

 # Significance based on likelihood ratio test
 ratio_test <- lrtest(my_null,my_model)
 p_value=ratio_test$`Pr(>Chisq)`[2]
 SNP <- colnames(my_data)[i]
 line <- cbind(SNP, p_value)
 results <- rbind(results, line)
 models[[SNP]] <- my_model_summary$coefficients
}





######################### lme4qtl Maybe but no p-values, better suited to estimate heriability in GWAS
library(lmtest)
library(lme4)
library(lattice)
library(lme4qtl)

my_data <- read.table("final_gea_input.txt", sep = "\t", h = T)
num_columns <- ncol(my_data)
my_mat <- read.table("relatedness.txt")
my_mat <- as.matrix(my_mat)


# set up outputs
results <- data.frame()
coefficient_results <- data.frame()
models <- list()


for (i in 22:num_columns) {

 my_model <- relmatLmer(my_data[,i] ~ Topography  * H_Index * PC1 * (1|Name), my_data, relmat = list(Name = my_mat))
 my_model_summary <- summary(my_model)

 my_null <- relmatLmer(my_data[,i] ~ H_Index * PC1 * (1|Name), my_data, relmat = list(Name = my_mat))
 my_null_summary <- summary(my_null)

 # Significance based on likelihood ratio test
 ratio_test <- lrtest(my_null,my_model)
 p_value=ratio_test$`Pr(>Chisq)`[2]
 SNP <- colnames(my_data)[i]
 line <- cbind(SNP, p_value)
 results <- rbind(results, line)

}
