library(ggplot2)
library(gridExtra)
library(viridis)

# Load your data
gabriele <- readRDS("1.Emerald_ReOrder_window_1_window_stats_seed_1001.rds")

# Vector indicating which column to use for color filling
my_states <- c(20, 28, 24, 24, 24, 24, 24, 28, 24, 20, 24, 24)

# Create a function to generate Manhattan plots
manhattan_plot <- function(data, color_col, plot_title) {
  # Convert data to data frame
  df <- as.data.frame(data)
  # Create Manhattan plot
  p <- ggplot(df, aes(x = V1, y = V7, color = factor(df[, color_col]))) +
    geom_point() +
    labs(title = plot_title,
         x = "Chromosome position",
         y = "-log10 p-value",
         color = "States")  + 
         scale_color_viridis(discrete= TRUE, direction = -1) +	 # Change legend title to "States"
    theme(legend.key.size = unit(0.2, "cm")) # Adjust the size of the legend
  
  # Convert ggplot object to grob
  grob <- ggplotGrob(p)
  return(grob)
}

# Create a list to store plots
plots <- list()

# Generate Manhattan plots for each entry in the list
for (i in 1:length(gabriele)) {
    a <- manhattan_plot(gabriele[[i]], my_states[i], paste("Chromosome", i))
    plots[[i]] <- a
}

# Arrange and print the plots on the same page
grid.arrange(grobs = plots, ncol = 2)
