library(corrr)
library(ggcorrplot)
library(factoextra)
env_data <- read.table("metadata_ordered_header.txt", h = T)
env_data <- env_data[,4:25]
colSums(is.na(env_data))
data_normalized <- scale(env_data)

pdf("pca_summary.pdf")

data.pca <- princomp(data_normalized)
summary(data.pca)
data.pca$loadings[, 1:2]
fviz_eig(data.pca, addlabels = TRUE)
# Graph of the variables
fviz_pca_var(data.pca, col.var = "black")
fviz_cos2(data.pca, choice = "var", axes = 1:2)
fviz_pca_var(data.pca, col.var = "cos2",
            gradient.cols = c("black", "orange", "green"),
            repel = TRUE)
			
species <- read.table("species_assignment_0.9.txt")



library(ggplot2)

# Read species assignment data
species <- read.table("species_assignment_0.9.txt")

# Define custom colors
custom_colors <- c("engelmann" = "darkgreen", "glauca" = "lightgrey", "hybrid" = "yellow")

# Plot with ggplot
ggplot(as.data.frame(data.pca$scores), aes(x = Comp.1, y = Comp.2, color = factor(species$V1))) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +
  labs(color = "Species") +
  theme_minimal()
 
 
dev.off()
names <- read.table("names.txt")
write.table(cbind(names$V1,data.pca$scores[,1:2]), file = "pca_eigenvectors.txt", quote = F, row.names = F, col.names = F, sep = "\t")