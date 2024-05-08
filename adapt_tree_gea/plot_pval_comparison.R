a <- readRDS("banff_likelihood_ratio_results.rds")
b <- readRDS("adapt_tree_likelihood_ratio_results_PC1.rds")
c <- readRDS("adapt_tree_likelihood_ratio_results_PC2.rds")
a$p_value <- as.numeric(as.character(a$p_value))
b$p_value <- as.numeric(as.character(b$p_value))
c$p_value <- as.numeric(as.character(c$p_value))
z <- cbind(a$p_value,b$p_value,c$p_value)



colnames(z)<- c("banff","adapt_pc1", "adapt_pc2")
z <- as.data.frame(z)


library(gridExtra)
library(ggplot2)

# Create individual ggplot objects
plot1 <- ggplot(z, aes(x = -log10(banff), y = -log10(adapt_pc1))) +
  geom_point() +
  labs(title = "banff vs adapt_pc1")

plot2 <- ggplot(z, aes(x = -log10(banff), y = -log10(adapt_pc2))) +
  geom_point() +
  labs(title = "banff vs adapt_pc2")

plot3 <- ggplot(z, aes(x = -log10(adapt_pc1), y = -log10(adapt_pc2))) +
  geom_point() +
  labs(title = "adapt_pc1 vs adapt_pc2")

# Arrange plots using grid.arrange
grid.arrange(plot1, plot2, plot3, ncol = 1)