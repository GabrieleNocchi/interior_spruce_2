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

merged_data_2 <- rbind(merged_data_a,merged_data_b,merged_data_c,merged_data_d,merged_data_e,merged_data_f,merged_data_g,merged_data_h)


lg_groups <- split(merged_data, merged_data$LG)


lg_groups_2 <- split(merged_data_2, merged_data_2$LG)


for (lg_name in names(lg_groups)) {
  
  lg_data <- lg_groups[[lg_name]]
  lg_data_2 <- lg_groups_2[[lg_name]]
  lg_data$Transect <- factor(lg_data$Transect, levels = c("00. AdaptTree", "0. Overall", "1. Emerald", "5. Cascade","2. Cathedral",  "6. Canmore", "3. Bell",  "7. Yates", "4. Norquay", "8. Old Baldy"))
  
  plot_list <- list()
  
  for (tra_name in unique(lg_data$Transect)) {
    lg_temp <- lg_data[lg_data$Transect == tra_name,]
    lg_temp$Transect <- factor(lg_temp$Transect, levels = c("00. AdaptTree", "0. Overall", "1. Emerald", "5. Cascade","2. Cathedral",  "6. Canmore", "3. Bell",  "7. Yates", "4. Norquay", "8. Old Baldy"))
    
# Set different y-axis limits based on the transect value
if (tra_name %in% c("0. Overall", "00. AdaptTree")) {
  y_limit <- max(-log10(lg_temp$p_value), na.rm = TRUE)
} else {
  max_p_value <- max(-log10(lg_data_2$p_value), na.rm = TRUE)
  if (max_p_value > 50) {
    y_limit <- 15
  } else {
    y_limit <- max_p_value
  }
}
    
    # Create the Manhattan plot
    manhattan_plot <- ggplot(lg_temp, aes(x = as.numeric(index), y = -log10(p_value))) +
      geom_point(size = 1) +
      labs(x = "SNP Index", y = "-log10(p-value)") +
      theme_minimal() +
      ggtitle(tra_name) + 
      ylim(0, y_limit)
    
    plot_list[[tra_name]] <- manhattan_plot
  }
  
  # Save the plot to a PDF file
  pdf(paste(lg_name,"_standardized_axis.pdf"))
  plot_list_ordered <- plot_list[levels(lg_data$Transect)]
  print(do.call(grid.arrange, c(plot_list_ordered, ncol = 2)))
  dev.off()
}