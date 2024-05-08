library(ggplot2)
library(gridExtra)

a <- read.table("forR_best_hit.txt")
b <- read.table("forR_second_best.txt")


# Create histograms using ggplot2
histogram_a <- ggplot(a, aes(x = V11)) +
  geom_histogram(bins = 200, fill = "blue", alpha = 0.5) +
  ggtitle("Histogram of best hit identity")

histogram_b <- ggplot(b, aes(x = V11)) +
  geom_histogram(bins = 200, fill = "red", alpha = 0.5) +
  ggtitle("Histogram of second best hit identity")


  
# Create histograms using ggplot2
histogram_aa <- ggplot(a, aes(x = V10)) +
  geom_histogram(bins = 200, fill = "blue", alpha = 0.5) +
  ggtitle("Histogram of best hit coverage")

histogram_bb <- ggplot(b, aes(x = V10)) +
  geom_histogram(bins = 200, fill = "red", alpha = 0.5) +
  ggtitle("Histogram of second best hit coverage")


  
  
  
  

# Arrange histograms using grid.arrange
grid.arrange(histogram_a,histogram_aa, histogram_b,histogram_bb, ncol = 2)