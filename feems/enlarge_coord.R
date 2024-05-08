a <- read.table("*coord")

points <- as.matrix(a)

# Calculate center of the points
center <- colMeans(points)

# Calculate distances from the center
distances <- sqrt((points[,1] - center[1])^2 + (points[,2] - center[2])^2)

# Determine scale factor
scale_factor <- 12  # Change this to desired scale factor

# Scale distances
scaled_distances <- distances * scale_factor

# Calculate new coordinates based on scaled distances
new_points <- matrix(NA, nrow = nrow(points), ncol = 2)
for (i in 1:nrow(points)) {
  direction_vector <- points[i,] - center
  direction_vector <- direction_vector / sqrt(sum(direction_vector^2))  # Normalize direction vector
  new_points[i,] <- center + scaled_distances[i] * direction_vector
}

new_points <- as.data.frame(new_points)

write.table(new_points, file = "enlarged.coord", col.names = F, quote = F, row.names = F)