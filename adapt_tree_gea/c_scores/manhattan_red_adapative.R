# Manhattan Plot 2 based on likelihood ratio test p values
library(data.table)
library(gridExtra)
library(ggplot2)
library(dplyr)

a <- readRDS("likelihood_ratio_results_1.Emerald.rds")
a$p_value <- as.numeric(as.character(a$p_value))
a$SNP <- gsub("_.*", "", a$SNP)

# Format Pg- to Pg."
lg <- read.table("final_order.txt")


# Perform left join
merged_data <- left_join(a, lg, by = c("SNP" = "V4"))
colnames(merged_data) <- c("SNP", "p_value","LG","cM","cDNA")



# Load necessary libraries manhattan by LG
library(ggplot2)

# Create the Manhattan plot
manhattan_plot <- ggplot(merged_data, aes(x = as.numeric(row.names(merged_data)), y = -log10(p_value))) +
    geom_point(aes(color = -log10(p_value) > 4), size = 1) +
    scale_color_manual(values = c("black", "red"), guide = FALSE) +
    labs(x = "SNP Index", y = "-log10(p-value)", title = "Manhattan Plot") +
    theme_minimal() +
    facet_wrap(~ LG, scales = "free", ncol = 2)

# Save the plot to a PDF file
pdf("emerald_manhattan_likelihood_ratio_by_LG.pdf")
print(manhattan_plot)
dev.off()


# -log10 p above 10 for banff
# -log10 p above 15 adapt_tree