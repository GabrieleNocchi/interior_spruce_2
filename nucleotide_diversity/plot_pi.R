# Load required libraries
library(ggplot2)
library(gridExtra)

# Function to read file and create boxplot
create_boxplot <- function(file_path, color) {
    # Read data
    data <- read.table(file_path, header = TRUE)
    
    # Create boxplot
    p <- ggplot(data, aes(x = "", y = PI)) +
        geom_boxplot(fill = color) +  # Specify fill color
        labs(title = basename(file_path)) +  # Set title to file name
        theme_minimal() +
        theme(axis.title.x=element_blank(), 
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
    
    return(p)
}

# List of file names with suffix ".pi"
file_names <- list.files(pattern = "\\.pi$")

# Generate a vector of colors
colors <- rainbow(length(file_names))

# Create a list of plots with corresponding colors
plots <- Map(create_boxplot, file_names, colors)

# Arrange plots in a grid
grid.arrange(grobs = plots, ncol = 8)  # Change ncol as per your preference
