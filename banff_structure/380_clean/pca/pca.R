library(ggplot2)

# Load PCA data without header
pca_data <- read.table("380_pca.eigenvec", header = FALSE)

# Load metadata without header
metadata <- read.table("metadata.txt", header = FALSE, sep = "\t")

# Create PCA plot
p1 <- ggplot() +
    geom_point(data = pca_data, aes(x = V2, y = V3, color = metadata$V4)) +
    labs(title = "PC1 vs PC2",
         x = "PC1",
         y = "PC2",
         color = "Aspect") + theme_bw() +
  theme(legend.position = "none") 
		 
		 
# Create PCA plot
p2 <- ggplot() +
    geom_point(data = pca_data, aes(x = V2, y = V4, color = metadata$V4)) +
    labs(title = "PC1 vs PC3",
         x = "PC1",
         y = "PC3",
         color = "Aspect")+ theme_bw() +
  theme(legend.position = "none") 
		 
		 
# Create PCA plot
p3 <- ggplot() +
    geom_point(data = pca_data, aes(x = V3, y = V4, color = metadata$V4)) +
    labs(title = "PC2 vs PC3",
         x = "PC2",
         y = "PC3",
         color = "Aspect")+ theme_bw()
		 
		 
		 
		 
		 
# Create PCA plot
p4 <- ggplot() +
    geom_point(data = pca_data, aes(x = V2, y = V3, color = metadata$V3)) +
    labs(title = "PC1 vs PC2",
         x = "PC1",
         y = "PC2",
         color = "Transect")+ theme_bw() +
  theme(legend.position = "none") 
		 
		 
# Create PCA plot
p5 <- ggplot() +
    geom_point(data = pca_data, aes(x = V2, y = V4, color = metadata$V3)) +
    labs(title = "PC1 vs PC3",
         x = "PC1",
         y = "PC3",
         color = "Transect")+ theme_bw() +
  theme(legend.position = "none") 
		 
		 
# Create PCA plot
p6 <- ggplot() +
    geom_point(data = pca_data, aes(x = V3, y = V4, color = metadata$V3)) +
    labs(title = "PC2 vs PC3",
         x = "PC2",
         y = "PC3",
         color = "Transect")+ theme_bw()
		 
		 
		 # Read the data from the file
data <- read.table("380_pca.eigenval", header = FALSE)

# Create a line plot
p7 <- ggplot(data, aes(x = 1:10, y = V1)) +
  geom_line(color = "red") +
  labs(title = "PCA Eigenvalues",
       x = "Principal component",
       y = "Eigenvalue") + theme_bw()+geom_point(color = "black", size = 3) +  # Set y-axis limits using coord_cartesian
  scale_x_continuous(breaks = seq(0, 10, by = 1)) 
	   
	   
	   library(gridExtra)
	   library(grid)
	   
	   
	   
	   # Specify the widths and heights for the arrangement
widths <- c(1, 1, 1.3)
heights <- c(1, 1, 1)

# Create the arranged plot
arranged_plot <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7,
                             ncol = 3,
                             layout_matrix = rbind(c(1, 2, 3),
                                                  c(4, 5, 6),
                                                  c(7, 7, 7)),
                             widths = widths,
                             heights = heights)

# Draw the arranged plot
grid.draw(arranged_plot)