a <- read.table("relatedness.txt")
a <-as.matrix(a)
custom_order <- read.table("names_ordered_by_pop.txt")
custom_order <- as.vector(custom_order$V1)
a <- a[custom_order, custom_order, drop = FALSE]


library(ggplot2)
# Create a ggplot heatmap with custom settings
ggplot(data = reshape2::melt(a), aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(width = 1) +
    scale_fill_fermenter(palette = "YlGnBu",breaks = seq(0, 1, by = 0.2),  limits = c(-0.2, 1.2), na.value = "grey", direction = 1) +
    theme_minimal() +
    labs(title = "Relatedness", fill = "Relatedness") +
    theme(axis.text = element_blank(),
          axis.title = element_blank())
