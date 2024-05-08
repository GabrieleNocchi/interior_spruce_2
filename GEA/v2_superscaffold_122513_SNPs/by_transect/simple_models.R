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

# set up outputs
results <- data.frame()
coefficient_results <- data.frame()
model_pvals <- data.frame()
models <- list()

# Loop through columns
for (i in 22:num_columns) {
  
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
  
  # Extracting p-values for Topography coefficients and interactions involving topography and taking the minimum
  p_values_coeff <- my_model_summary$coefficients[c(2:3,5:6), 4]
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
saveRDS(models, file = "models.rds")

# Save results
saveRDS(results, file = "likelihood_ratio_results.rds")

# Save coefficient_results
saveRDS(coefficient_results, file = "coefficient_results.rds")

saveRDS(model_pvals, file = "model_pval_results.rds")




##### Plot Top SNPs


library(data.table)
library(gridExtra)
library(ggplot2)

a <- readRDS("coefficient_results.rds")
a$coeff_p_value <- as.numeric(as.character(a$coeff_p_value))
sorted_a <- a[order(a$coeff_p_value), ]

b <- readRDS("likelihood_ratio_results.rds")
b$p_value <- as.numeric(as.character(b$p_value))
sorted_b <- b[order(b$p_value), ]

c <- readRDS("model_pval_results.rds")
c$p_value <- as.numeric(as.character(c$p_value))
sorted_c <- c[order(c$p_value), ]




### Based on model like ratio p values
indices <- as.numeric(rownames(sorted_b))
indices <- indices + 21
top <- head(indices)

df <- fread("final_gea_input.txt", sep = "\t", header = TRUE)



pdf("top_SNPs_likelihood_ratio.pdf")
# Create a list to store individual plots
plot_list <- list()
 for (column in top) {
 category <- "0. Overall"
# Create a histogram
p <- ggplot(df, aes(x = df[[4]], fill = as.factor(df[[column]]))) +
  geom_bar(position = "dodge") +
  labs(title = paste("Distribution of ",names(df)[column], "Overall"),
       x = paste("Column", column),
       y = "Frequency") +
  theme_minimal() +
		theme(axis.text.x = element_text(size = 5), # Adjust x-axis label size
            axis.text.y = element_text(size = 5), # Adjust y-axis label size
            axis.title.x = element_text(size = 5), # Adjust x-axis title size
            axis.title.y = element_text(size = 5), # Adjust y-axis title size
            plot.title = element_text(size = 5)) +	
		    theme(legend.position = "right", # Placing legend at the bottom
            legend.text = element_text(size = 5), # Adjust legend text size
            legend.title = element_blank()) # Remove legend title

# Modifying the legend label
p <- p + scale_fill_discrete(name = "Allele")

# Adjusting the size of the legend key
p <- p + guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

# Adjusting the size of the legend
p <- p + theme(legend.box.margin = margin(0, 0, 0, 0))
    
    # Store the plot in the plot_list
    plot_list[[category]] <- p
 
# Loop through unique categories in column 3
for (category in unique(df[[3]])) {
    # Subset the data for the current category
    subset_df <- df[df[[3]] == category, ]
	
	
    
    # Create a histogram for the current category
    p <- ggplot(subset_df, aes_string(x = names(subset_df)[4], fill = as.factor(subset_df[[column]]))) +
        geom_bar(position = "dodge") +
        labs(title = paste("Distribution of ",names(subset_df)[column], category),
             x = "Topography",
             y = "Frequency") +
        theme_minimal() +
		theme(axis.text.x = element_text(size = 5), # Adjust x-axis label size
            axis.text.y = element_text(size = 5), # Adjust y-axis label size
            axis.title.x = element_text(size = 5), # Adjust x-axis title size
            axis.title.y = element_text(size = 5), # Adjust y-axis title size
            plot.title = element_text(size = 5)) +	
		    theme(legend.position = "right", # Placing legend at the bottom
            legend.text = element_text(size = 5), # Adjust legend text size
            legend.title = element_blank()) # Remove legend title

# Modifying the legend label
p <- p + scale_fill_discrete(name = "Allele")

# Adjusting the size of the legend key
p <- p + guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

# Adjusting the size of the legend
p <- p + theme(legend.box.margin = margin(0, 0, 0, 0))
    
    # Store the plot in the plot_list
    plot_list[[category]] <- p
 
}



# Arrange plots in alphanumerical order based on category value
sorted_plot_list <- plot_list[order(names(plot_list))]

# Arrange plots using grid.arrange
grid.arrange(grobs = sorted_plot_list, ncol = 1) 
}
dev.off()






# Manhattan Plot based on coefficients p values
library(data.table)
library(gridExtra)
library(ggplot2)
library(dplyr)

a <- readRDS("coefficient_results.rds")
a$coeff_p_value <- as.numeric(as.character(a$coeff_p_value))
a$SNP <- gsub("_.*", "", a$SNP)

# Format Pg- to Pg."
lg <- read.table("final_order.txt")


# Perform left join
merged_data <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data) <- c("SNP", "coeff_p_value","LG","cM","cDNA")

# Load the ggplot2 package
library(ggplot2)
custom_colors <- c("grey", "black","grey", "black","grey", "black","grey", "black","grey", "black","grey", "black")
# Create a Manhattan plot
pdf("manhattan_coefficient.pdf")
ggplot(merged_data, aes(x = as.numeric(row.names(merged_data)), y = -log10(coeff_p_value), color = LG)) +
    geom_point() +
    scale_color_manual(values = custom_colors) +  # Color points based on V1
    labs(x = "SNP Index", y = "-log10(p-value)", title = "Manhattan Plot") +
    theme_minimal()

# Print the plot
dev.off()

# Load necessary libraries manhattan by LG
library(ggplot2)

# Create the Manhattan plot
manhattan_plot <- ggplot(merged_data, aes(x = as.numeric(row.names(merged_data)), y = -log10(coeff_p_value))) +
    geom_point(size = 1) +
    labs(x = "SNP Index", y = "-log10(p-value)", title = "Manhattan Plot") +
    theme_minimal() +
    facet_wrap(~ LG, scales = "free", ncol = 2)

# Save the plot to a PDF file
pdf("manhattan_coefficient_by_LG.pdf")
print(manhattan_plot)
dev.off()




# Manhattan Plot 2 based on likelihood ratio test p values
library(data.table)
library(gridExtra)
library(ggplot2)
library(dplyr)

a <- readRDS("likelihood_ratio_results.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)

# Format Pg- to Pg."
lg <- read.table("final_order.txt")


# Perform left join
merged_data <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data) <- c("SNP", "p_value","LG","cM","cDNA")

# Load the ggplot2 package
library(ggplot2)
custom_colors <- c("grey", "black","grey", "black","grey", "black","grey", "black","grey", "black","grey", "black")
# Create a Manhattan plot
pdf("manhattan_likelihood_ratio.pdf")
ggplot(merged_data, aes(x = as.numeric(row.names(merged_data)), y = -log10(p_value), color = LG)) +
    geom_point() +
    scale_color_manual(values = custom_colors) +  # Color points based on V1
    labs(x = "SNP Index", y = "-log10(p-value)", title = "Manhattan Plot") +
    theme_minimal()

# Print the plot
dev.off()


# Load necessary libraries manhattan by LG
library(ggplot2)

# Create the Manhattan plot
manhattan_plot <- ggplot(merged_data, aes(x = as.numeric(row.names(merged_data)), y = -log10(p_value))) +
    geom_point(size = 1) +
    labs(x = "SNP Index", y = "-log10(p-value)", title = "Manhattan Plot") +
    theme_minimal() +
    facet_wrap(~ LG, scales = "free", ncol = 2)

# Save the plot to a PDF file
pdf("manhattan_likelihood_ratio_by_LG.pdf")
print(manhattan_plot)
dev.off()