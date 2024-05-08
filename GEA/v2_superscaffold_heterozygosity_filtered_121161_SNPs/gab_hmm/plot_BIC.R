# Get list of files in working directory
files <- list.files(pattern = "4.Norquay_ReOrder_window_1_the_crits_list__seed_[0-9]*.rds")

# Iterate over each file
for (file in files) {
  # Extract number from file name
  num <- gsub("[^0-9]", "", file)
  num <- gsub("41", "", num)

a <- readRDS(file)

# Function to extract column 2 from each entry
get_column_2 <- function(entry) {
  entry[, 2]
}

# Extract column 2 from each entry in the list
column_2_list <- lapply(a, get_column_2)


pdf(paste("Norquay_",num, ".pdf"))

plots_list <- list()
# Plot each entry in column 2 and highlight minimum point
library(ggplot2)
library(gridExtra)

# Set global parameter for axis label size
theme_set(theme_bw(base_size = 8)) # Adjust the size as needed

# Plot each entry in column 2 and highlight the minimum point
for (i in 1:length(column_2_list)) {
  df <- data.frame(Index = 1:length(column_2_list[[i]]), Value = column_2_list[[i]])
  
  p <- ggplot(df, aes(x = Index, y = Value)) +
    geom_point() +
    labs(title = paste("Chromosome", i), x = "Index", y = "Value") +
    geom_point(data = df[which.min(df$Value), , drop = FALSE], aes(color = "Minimum"), size = 2) +
    scale_color_manual(values = c("Minimum" = "red")) +
    theme_bw(base_size = 8) + 
      scale_x_continuous(breaks = seq(1, length(a[[i]]), by = 1)) + theme(legend.position = "none")

    
  
    plots_list[[i]] <- p
}
grid_arrange = do.call("grid.arrange", c(plots_list, ncol = 3))
print(grid_arrange)
dev.off()
}