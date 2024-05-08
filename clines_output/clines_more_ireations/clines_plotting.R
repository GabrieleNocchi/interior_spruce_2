library("ClineHelpR")
library(dplyr)
plotDIR <- "./plots"
genesDIR <- "./results/"
prefix <- "super_scaffolds_thinned"
admixPop <- "admixed"
popMap <- "./data/clines_095.txt"
genes.loci.file <- genes.loci.file <- "data/super_scaffolds_thinned_loci.txt"


bgc.genes <-
  combine_bgc_output(results.dir = genesDIR,
                     prefix = prefix)
					 
plot_traces(df.list = bgc.genes,
            prefix = prefix,
            plotDIR = plotDIR, showPLOTS=TRUE)
				
gene.outliers <-
  get_bgc_outliers(
    df.list = bgc.genes, 
    admix.pop = admixPop, 
    popmap =  popMap,
    loci.file = genes.loci.file, 
    qn = 0.975)
	
	


my_res <- gene.outliers[[1]]

# Format Pg- to Pg."
lg <- read.table("final_order.txt")


# Perform left join
merged_data <- left_join(my_res, lg, by = c("X.CHROM" = "V4"))



# Replace NA with "neut" in alpha.excess column
merged_data$alpha.excess <- ifelse(is.na(merged_data$alpha.excess), "neut", merged_data$alpha.excess)

# Convert alpha.excess to factor with "neut" included
merged_data$alpha.excess <- factor(merged_data$alpha.excess, levels = c("neg", "pos", "neut"))





# Convert alpha.excess to factor with NA as a category
my_pos <- subset(merged_data, !is.na(alpha.excess) & alpha.excess == "pos")
my_neg <- subset(merged_data, !is.na(alpha.excess) & alpha.excess == "neg")
my_double <- subset(merged_data, crazy.a == "TRUE")
merged_data <- merged_data %>%
  rename(LG = V1)
  
 pdf("clines_man.pdf")
# Load the ggplot2 package
library(ggplot2)
custom_colors <- c("grey", "black","grey", "black","grey", "black","grey", "black","grey", "black","grey", "black")
# Create a Manhattan plot
ggplot(merged_data, aes(x = as.numeric(row.names(merged_data)), y = alpha, color = LG)) +
    geom_point() +
    scale_color_manual(values = custom_colors) +  # Color points based on LG
    labs(x = "SNP Index", y = "BGC alpha", title = "Manhattan Plot") +
    theme_minimal() +
    geom_point(data = my_pos, 
               aes(x = as.numeric(row.names(my_pos)), y = alpha), color ="deepskyblue1") +
    geom_point(data = my_neg, 
               aes(x = as.numeric(row.names(my_neg)), y = alpha), color ="darkgreen")+
    geom_point(data = my_double, 
               aes(x = as.numeric(row.names(my_double)), y = alpha), color ="red")
			   
dev.off()
			   

pdf("clines_man_LG.pdf")
# Create the Manhattan plot with different colors for alpha.excess
ggplot(merged_data, aes(x = as.numeric(row.names(merged_data)), y = alpha, color = alpha.excess, size = 0.5, shape = factor(crazy.a))) +
    geom_point(size = 1) +
    scale_color_manual(values = c("darkgreen", "deepskyblue1", "grey"), 
                       labels = c("Negative", "Positive", "Neutral"),
                       guide = guide_legend(title = "Alpha Outlier")) +
    scale_shape_manual(values = c(16, 17), 
                        labels = c("FALSE", "TRUE"),
                        guide = guide_legend(title = "Both Outlier Tests")) +
    labs(x = "SNP Index", y = "BGC alpha", title = "Manhattan Plot") +
    theme_minimal() +
    facet_wrap(~ LG, scales = "free", ncol = 2)

dev.off()






# one test only necessary to pass

phiPlot(
  outlier.list = gene.outliers,
  popname = paste0(admixPop, " Genes"),
  line.size = 0.35,
  saveToFile = paste0(prefix, "_genes"),
  plotDIR = plotDIR,
  showPLOTS = TRUE,
  both.outlier.tests = FALSE,
  neutral.color = "gray60",
  alpha.color = "cornflowerblue",
  beta.color = "firebrick",
  both.color = "purple",
  hist.y.origin = 1.2,
  hist.height = 1.8,
  margins = c(160.0, 5.5, 5.5, 5.5),
  hist.binwidth = 0.05
)


alphaBetaPlot(
  gene.outliers,
  alpha.color = "cornflowerblue",
  both.outlier.tests = FALSE,
  beta.color = "orange",
  neutral.color = "gray60",
  saveToFile = prefix,
  plotDIR = plotDIR,
  showPLOTS = TRUE,
  padding = 0.2,
)