library(ggplot2)
library(viridis)


data <- read.table("metadata.txt",sep = "\t")


colnames(data) <- c("ID", "Code", "Category", "Type")

ggplot(data, aes(x = Category, fill = Type)) +
  geom_bar(position = "dodge",color = "black") +
  scale_fill_viridis(discrete = TRUE) +  # Use viridis color palette
  labs(title = "Spruce Sampling",
       x = "Transect",
       y = "Sampled individuals") +
  coord_cartesian(ylim = c(0, 20)) +  # Set y-axis limits using coord_cartesian
  scale_y_continuous(breaks = seq(0, 20, by = 2)) +  # Set y-axis ticks every 2
  theme_minimal()