# Manhattan Plot 2 based on likelihood ratio test p values
library(data.table)
library(gridExtra)
library(ggplot2)
library(dplyr)

# Define the initial data frame with proper row and column specifications
my_c_input <- data.frame(matrix(NA, nrow = 43643, ncol = 8))

a <- readRDS("likelihood_ratio_results_1.Emerald.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
# Convert each value of column p_value to 0 and 1
a$emerald <- ifelse(is.na(a$p_value) | a$p_value >= 0.0001, 0, 1)
my_c_input[,1] <-  a$emerald
emerald <- as.vector(a$emerald)

a <- readRDS("likelihood_ratio_results_2.Cathedral.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
# Convert each value of column p_value to 0 and 1
a$cathedral <- ifelse(is.na(a$p_value) | a$p_value >= 0.0001, 0, 1)
my_c_input[,2] <-  a$cathedral
cathedral <- as.vector(a$cathedral)

a <- readRDS("likelihood_ratio_results_3.Bell.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
# Convert each value of column p_value to 0 and 1
a$bell <- ifelse(is.na(a$p_value) | a$p_value >= 0.0001, 0, 1)
my_c_input[,3] <-  a$bell
bell <- as.vector(a$bell)

a <- readRDS("likelihood_ratio_results_4.Norquay.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
# Convert each value of column p_value to 0 and 1
a$norquay <- ifelse(is.na(a$p_value) | a$p_value >= 0.0001, 0, 1)
my_c_input[,4] <-  a$norquay
norquay <- as.vector(a$norquay)

a <- readRDS("likelihood_ratio_results_5.Cascade.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
# Convert each value of column p_value to 0 and 1
a$cascade <- ifelse(is.na(a$p_value) | a$p_value >= 0.0001, 0, 1)
my_c_input[,5] <-  a$cascade
cascade <- as.vector(a$cascade)

a <- readRDS("likelihood_ratio_results_6.Canmore.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
# Convert each value of column p_value to 0 and 1
a$canmore <- ifelse(is.na(a$p_value) | a$p_value >= 0.0001, 0, 1)
my_c_input[,6] <-  a$canmore
canmore <- as.vector(a$canmore)

a <- readRDS("likelihood_ratio_results_7.Yates.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
# Convert each value of column p_value to 0 and 1
a$yates <- ifelse(is.na(a$p_value) | a$p_value >= 0.0001, 0, 1)
my_c_input[,7] <-  a$yates
yates <- as.vector(a$yates)

a <- readRDS("likelihood_ratio_results_8.Old Baldy.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
# Convert each value of column p_value to 0 and 1
a$old <- ifelse(is.na(a$p_value) | a$p_value >= 0.0001, 0, 1)
my_c_input[,8] <-  a$old
old <- as.vector(a$old)

################################################## Adjust threshold
a <- readRDS("banff_likelihood_ratio_results.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
# Convert each value of column p_value to 0 and 1
a$all <- ifelse(is.na(a$p_value) | a$p_value >= 0.0000000001, 0, 1)
all <- as.vector(a$all)


colnames(my_c_input) <- c("emerald","cathedral","bell","norquay","cascade","canmore","yates","oldbaldy")


library(dgconstraint)


m <- as.matrix(my_c_input)
c_score_all <-  pairwise_c_hyper(m, na.rm = F)


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


all_results_banff <- results_df[results_df$Vector2 == "all",]
transects_results_banff <- results_df[results_df$Vector2 != "all",]




#############################################################################################àà Adapt Tree -- adjust threshold
a <- readRDS("adapt_tree_likelihood_ratio_results_PC2.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
# Convert each value of column p_value to 0 and 1
a$adapt_pc2 <- ifelse(is.na(a$p_value) | a$p_value >= 0.000000000000001, 0, 1)
adapt_pc2 <- as.vector(a$adapt_pc2)


colnames(my_c_input) <- c("emerald","cathedral","bell","norquay","cascade","canmore","yates","oldbaldy")


library(dgconstraint)


m <- as.matrix(my_c_input)
c_score_all <-  pairwise_c_hyper(m, na.rm = F)


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
  adapt_pc2 = adapt_pc2
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


all_results_pc2 <- results_df[results_df$Vector2 == "adapt_pc2",]
transects_results_pc2  <- results_df[results_df$Vector2 != "adapt_pc2",]







#### PC1





#############################################################################################àà Adapt Tree -- adjust threshold
a <- readRDS("adapt_tree_likelihood_ratio_results_PC1.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)
# Convert each value of column p_value to 0 and 1
a$adapt_pc1 <- ifelse(is.na(a$p_value) | a$p_value >= 0.000000000000001, 0, 1)
adapt_pc1 <- as.vector(a$adapt_pc1)


colnames(my_c_input) <- c("emerald","cathedral","bell","norquay","cascade","canmore","yates","oldbaldy")


library(dgconstraint)


m <- as.matrix(my_c_input)
c_score_all <-  pairwise_c_hyper(m, na.rm = F)


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
  adapt_pc1 = adapt_pc1
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


all_results_pc1 <- results_df[results_df$Vector2 == "adapt_pc1",]
transects_results_pc1 <- results_df[results_df$Vector2 != "adapt_pc1",]







#### Banff vs adapt PC1 vs adapt pc2


# Display the results
print(all_results_banff)
#print(transects_results_banff)
print(all_results_pc1)
#print(transects_results_pc1)
print(all_results_pc2)
print(transects_results_pc2)
print(c_score_all)





#### Banff vs adapt PC1 vs adapt pc2
single_c_hyper(all, adapt_pc1, na.rm = FALSE)
single_c_hyper(all, adapt_pc2, na.rm = FALSE)
single_c_hyper(adapt_pc1, adapt_pc2, na.rm = FALSE)