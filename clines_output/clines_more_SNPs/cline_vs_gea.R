a <- readRDS("likelihood_ratio_results_het_filtered.rds")
# Convert 'SNP' column to character type if it's not already
a$SNP <- as.character(a$SNP)

# Remove the final 2 characters from the 'SNP' values
a$SNP <- substring(a$SNP, 1, nchar(a$SNP) - 2)





library("ClineHelpR")
library(dplyr)
plotDIR <- "./plots"
genesDIR <- "./results/"
prefix <- "super_scaffolds_thinned_3_snp_heterozygosity_filtered"
admixPop <- "admixed"
popMap <- "./data/clines_095.txt"
genes.loci.file <- genes.loci.file <- "data/super_scaffolds_thinned_3_snp_heterozygosity_filtered_loci.txt"


bgc.genes <-
  combine_bgc_output(results.dir = genesDIR,
                     prefix = prefix)
					 

				
gene.outliers <-
  get_bgc_outliers(
    df.list = bgc.genes, 
    admix.pop = admixPop, 
    popmap =  popMap,
    loci.file = genes.loci.file, 
    qn = 0.975)
	
	


my_res <- gene.outliers[[1]]



your_data_frame <- my_res
# Create a new column 'mylabel' by joining 'X.CHROM' and 'POS' with an underscore
your_data_frame$mylabel <- paste(your_data_frame$X.CHROM, your_data_frame$POS, sep = "_")

# Substitute the dash after "Pg" with a dot in the mylabel column
your_data_frame$mylabel <- gsub("Pg-", "Pg.", your_data_frame$mylabel)



# Merge the data frames based on common columns SNP and mylabel
merged_data <- merge(your_data_frame, a, by.x = "mylabel", by.y = "SNP")

# Print the merged data
print(merged_data)


# Replace NA with "neut" in alpha.excess column
merged_data$alpha.excess <- ifelse(is.na(merged_data$alpha.excess), "neut", merged_data$alpha.excess)

# Convert alpha.excess to factor with "neut" included
merged_data$alpha.excess <- factor(merged_data$alpha.excess, levels = c("neg", "pos", "neut"))





pdf("gea_vs_clines.pdf")
library(ggplot2)

	


ggplot(merged_data, aes(x = as.numeric(as.character(alpha)), y = as.numeric(as.character(p_value)), shape = as.factor(crazy.a))) +
    geom_point(aes(color = alpha.excess)) +
    # Manipulate y-axis ticks
    scale_y_continuous(breaks = seq(0, max(as.numeric(as.character(merged_data$p_value))), by = 0.1)) +
    scale_color_manual(values = c("darkgreen", "deepskyblue1", "grey"), 
                       labels = c("Negative", "Positive", "Neutral"),
                       guide = guide_legend(title = "Alpha Outlier")) +
    scale_shape_manual(values = c(16, 17), 
                        labels = c("FALSE", "TRUE"),
                        guide = guide_legend(title = "Both Outlier Tests"))+
  labs(x = "BGC alpha", y = "GEA likelihood-ratio test p-val")


dev.off()

pdf("gea_vs_clines_log.pdf")
library(ggplot2)

	


ggplot(merged_data, aes(x = as.numeric(as.character(alpha)), y = -log10(as.numeric(as.character(p_value))), shape = as.factor(crazy.a))) +
    geom_point(aes(color = alpha.excess)) +
    # Manipulate y-axis ticks
    scale_y_continuous(breaks = seq(0, max(-log10(as.numeric(as.character(merged_data$p_value)))), by = 1)) +
    scale_color_manual(values = c("darkgreen", "deepskyblue1", "grey"), 
                       labels = c("Negative", "Positive", "Neutral"),
                       guide = guide_legend(title = "Alpha Outlier")) +
    scale_shape_manual(values = c(16, 17), 
                        labels = c("FALSE", "TRUE"),
                        guide = guide_legend(title = "Both Outlier Tests")) +
  labs(x = "BGC alpha", y = "-log10(GEA likelihood-ratio test p-val)")


dev.off()