a <- readRDS("1.Emerald_ReOrder_window_1_the_crits_list__seed_1001.rds")
b <- readRDS("cohen/1.Emerald_ReOrder_window_1_bumper_0_opt_cW.rds")


# Define a function to extract the row number of the minimum value in the second column
get_row_of_min_second_column <- function(mat) {
  which.min(mat[, 2])
}

# Apply the function to each entry in the list 'a'
min_row_numbers <- lapply(a, get_row_of_min_second_column)


# Subtract 1 from each value in the list
adjusted_values <- sapply(min_row_numbers, function(x) x - 1)


# Extract the corresponding row number of column 1 for each entry in list 'b'
extracted_rows <- lapply(seq_along(adjusted_values), function(i) {
  b[[i]][adjusted_values[i], 1]
})



# Unlist extracted rows
unlisted_rows <- unlist(extracted_rows)

# Create a dataframe with a column containing the unlisted rows
cohen_df <- data.frame(emerald = unlisted_rows)






a <- readRDS("2.Cathedral_ReOrder_window_1_the_crits_list__seed_1001.rds")
b <- readRDS("cohen/2.Cathedral_ReOrder_window_1_bumper_0_opt_cW.rds")


# Define a function to extract the row number of the minimum value in the second column
get_row_of_min_second_column <- function(mat) {
  which.min(mat[, 2])
}

# Apply the function to each entry in the list 'a'
min_row_numbers <- lapply(a, get_row_of_min_second_column)


# Subtract 1 from each value in the list
adjusted_values <- sapply(min_row_numbers, function(x) x - 1)


# Extract the corresponding row number of column 1 for each entry in list 'b'
extracted_rows <- lapply(seq_along(adjusted_values), function(i) {
  b[[i]][adjusted_values[i], 1]
})



# Unlist extracted rows
unlisted_rows <- unlist(extracted_rows)

# Create a dataframe with a column containing the unlisted rows
df <- data.frame(cathedral = unlisted_rows)
cohen_df <- cbind(cohen_df,df)











a <- readRDS("3.Bell_ReOrder_window_1_the_crits_list__seed_1001.rds")
b <- readRDS("cohen/3.Bell_ReOrder_window_1_bumper_0_opt_cW.rds")


# Define a function to extract the row number of the minimum value in the second column
get_row_of_min_second_column <- function(mat) {
  which.min(mat[, 2])
}

# Apply the function to each entry in the list 'a'
min_row_numbers <- lapply(a, get_row_of_min_second_column)



# Subtract 1 from each value in the list
adjusted_values <- sapply(min_row_numbers, function(x) x - 1)


# Extract the corresponding row number of column 1 for each entry in list 'b'
extracted_rows <- lapply(seq_along(adjusted_values), function(i) {
  b[[i]][adjusted_values[i], 1]
})



# Unlist extracted rows
unlisted_rows <- unlist(extracted_rows)

# Create a dataframe with a column containing the unlisted rows
df <- data.frame(bell = unlisted_rows)
cohen_df <- cbind(cohen_df,df)




a <- readRDS("4.Norquay_ReOrder_window_1_the_crits_list__seed_1001.rds")
b <- readRDS("cohen/4.Norquay_ReOrder_window_1_bumper_0_opt_cW.rds")


# Define a function to extract the row number of the minimum value in the second column
get_row_of_min_second_column <- function(mat) {
  which.min(mat[, 2])
}

# Apply the function to each entry in the list 'a'
min_row_numbers <- lapply(a, get_row_of_min_second_column)


# Subtract 1 from each value in the list
adjusted_values <- sapply(min_row_numbers, function(x) x - 1)


# Extract the corresponding row number of column 1 for each entry in list 'b'
extracted_rows <- lapply(seq_along(adjusted_values), function(i) {
  b[[i]][adjusted_values[i], 1]
})



# Unlist extracted rows
unlisted_rows <- unlist(extracted_rows)

# Create a dataframe with a column containing the unlisted rows
df <- data.frame(norquay = unlisted_rows)
cohen_df <- cbind(cohen_df,df)





a <- readRDS("5.Cascade_ReOrder_window_1_the_crits_list__seed_1001.rds")
b <- readRDS("cohen/5.Cascade_ReOrder_window_1_bumper_0_opt_cW.rds")


# Define a function to extract the row number of the minimum value in the second column
get_row_of_min_second_column <- function(mat) {
  which.min(mat[, 2])
}

# Apply the function to each entry in the list 'a'
min_row_numbers <- lapply(a, get_row_of_min_second_column)


# Subtract 1 from each value in the list
adjusted_values <- sapply(min_row_numbers, function(x) x - 1)


# Extract the corresponding row number of column 1 for each entry in list 'b'
extracted_rows <- lapply(seq_along(adjusted_values), function(i) {
  b[[i]][adjusted_values[i], 1]
})



# Unlist extracted rows
unlisted_rows <- unlist(extracted_rows)

# Create a dataframe with a column containing the unlisted rows
df <- data.frame(cascade = unlisted_rows)
cohen_df <- cbind(cohen_df,df)












a <- readRDS("6.Canmore_ReOrder_window_1_the_crits_list__seed_1001.rds")
b <- readRDS("cohen/6.Canmore_ReOrder_window_1_bumper_0_opt_cW.rds")


# Define a function to extract the row number of the minimum value in the second column
get_row_of_min_second_column <- function(mat) {
  which.min(mat[, 2])
}

# Apply the function to each entry in the list 'a'
min_row_numbers <- lapply(a, get_row_of_min_second_column)



# Subtract 1 from each value in the list
adjusted_values <- sapply(min_row_numbers, function(x) x - 1)


# Extract the corresponding row number of column 1 for each entry in list 'b'
extracted_rows <- lapply(seq_along(adjusted_values), function(i) {
  b[[i]][adjusted_values[i], 1]
})



# Unlist extracted rows
unlisted_rows <- unlist(extracted_rows)

# Create a dataframe with a column containing the unlisted rows
df <- data.frame(canmore = unlisted_rows)
cohen_df <- cbind(cohen_df,df)







a <- readRDS("7.Yates_ReOrder_window_1_the_crits_list__seed_1001.rds")
b <- readRDS("cohen/7.Yates_ReOrder_window_1_bumper_0_opt_cW.rds")


# Define a function to extract the row number of the minimum value in the second column
get_row_of_min_second_column <- function(mat) {
  which.min(mat[, 2])
}

# Apply the function to each entry in the list 'a'
min_row_numbers <- lapply(a, get_row_of_min_second_column)


# Subtract 1 from each value in the list
adjusted_values <- sapply(min_row_numbers, function(x) x - 1)


# Extract the corresponding row number of column 1 for each entry in list 'b'
extracted_rows <- lapply(seq_along(adjusted_values), function(i) {
  b[[i]][adjusted_values[i], 1]
})



# Unlist extracted rows
unlisted_rows <- unlist(extracted_rows)

# Create a dataframe with a column containing the unlisted rows
df <- data.frame(yates = unlisted_rows)
cohen_df <- cbind(cohen_df,df)








a <- readRDS("8.Old_Baldy_ReOrder_window_1_the_crits_list__seed_1001.rds")
b <- readRDS("cohen/8.Old_Baldy_ReOrder_window_1_bumper_0_opt_cW.rds")


# Define a function to extract the row number of the minimum value in the second column
get_row_of_min_second_column <- function(mat) {
  which.min(mat[, 2])
}

# Apply the function to each entry in the list 'a'
min_row_numbers <- lapply(a, get_row_of_min_second_column)


# Subtract 1 from each value in the list
adjusted_values <- sapply(min_row_numbers, function(x) x - 1)


# Extract the corresponding row number of column 1 for each entry in list 'b'
extracted_rows <- lapply(seq_along(adjusted_values), function(i) {
  b[[i]][adjusted_values[i], 1]
})



# Unlist extracted rows
unlisted_rows <- unlist(extracted_rows)

# Create a dataframe with a column containing the unlisted rows
df <- data.frame(old = unlisted_rows)
cohen_df <- cbind(cohen_df,df)






a <- readRDS("0.All_ReOrder_window_1_the_crits_list__seed_1001.rds")
b <- readRDS("cohen/0.All_ReOrder_window_1_bumper_0_opt_cW.rds")


# Define a function to extract the row number of the minimum value in the second column
get_row_of_min_second_column <- function(mat) {
  which.min(mat[, 2])
}

# Apply the function to each entry in the list 'a'
min_row_numbers <- lapply(a, get_row_of_min_second_column)


# Subtract 1 from each value in the list
adjusted_values <- sapply(min_row_numbers, function(x) x - 1)


# Extract the corresponding row number of column 1 for each entry in list 'b'
extracted_rows <- lapply(seq_along(adjusted_values), function(i) {
  b[[i]][adjusted_values[i], 1]
})



# Unlist extracted rows
unlisted_rows <- unlist(extracted_rows)

# Create a dataframe with a column containing the unlisted rows
df <- data.frame(all = unlisted_rows)
cohen_df <- cbind(cohen_df,df)



print(cohen_df)



library(ggplot2)
library(tidyr)

# Define the order of columns
column_order <- names(cohen_df)

# Reshape the dataframe to long format
cohen_df_long <- gather(cohen_df, key = "Variable", value = "Value")

# Convert Variable to factor with desired order
cohen_df_long$Variable <- factor(cohen_df_long$Variable, levels = column_order)

# Plot
# Convert row names to factor with desired order
row_names_factor <- factor(row.names(cohen_df), levels = row.names(cohen_df))

# Plot
ggplot(cohen_df_long, aes(x = Variable, y = Value)) +
  geom_point(aes(color = rep(row_names_factor, times = 9)), size = 1) +  
  labs(title = "Cohen's W", x = "Transects", y = "Cohen's W", color = "Chromosome") +
  scale_color_discrete(labels = row.names(cohen_df)) +  # Set legend labels
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability
  
  
 # With lines
  ggplot(cohen_df_long, aes(x = Variable, y = Value)) +
  geom_point(aes(color = rep(row_names_factor, times = 9)), size = 1) +
  geom_line(aes(group = rep(row_names_factor, times = 9), color = rep(row_names_factor, times = 9)), size = 0.5) + # Add lines connecting points of the same color
  labs(title = "Cohen's W", x = "Transects", y = "Cohen's W", color = "Chromosome") +
  scale_color_discrete(labels = row.names(cohen_df)) +  # Set legend labels
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability




gab_column <- rep(row_names_factor, 9)
gab_df <- cbind(cohen_df_long,gab_column)


# Create separate plots for each factor with 2 columns
ggplot(gab_df, aes(x = Variable, y = Value)) +
    geom_point(size = 1) +  
    geom_line(aes(group = 1), size = 0.5) +  # Add line connecting points
    labs(title = "Cohen's W", x = "Transects", y = "Cohen's W") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Rotate x-axis labels for better readability
    facet_wrap(~ gab_column, ncol = 3, scales = "free")
  
# Assuming cohen_df is your dataframe
mean_values <- colMeans(cohen_df, na.rm = TRUE)








#####################################################################################################################################################################################################
















#################3D Plots  with C scores -- all Chrom together with genome wide c scores and mean Ws
library(plot3D)

### Read genome wide c scores first....--> transects_results and all_results
# FORMATTING

transects_results_backup <- transects_results
all_results_backup <- all_results

transects_results$W1 <- mean_values[transects_results$Vector1]
transects_results$W2 <- mean_values[transects_results$Vector2]

all_results$W1 <- mean_values[all_results$Vector1]
all_results$W2 <- mean_values[all_results$Vector2]


# Load the plotly package
############################## PLOTS START HERE
library(plotly)

# Create the plotly 3D scatter plot for transect vs transect
plot_ly(data = transects_results, 
        x = ~W1, y = ~W2, z = ~Score, 
        type = "scatter3d", mode = "markers",
                  text = ~paste(Vector1, Vector2),
                  marker = list(symbol = "circle", size = 5)) %>%
    layout(scene = list(xaxis = list(title = "mean-W1"),
                        yaxis = list(title = "mean-W2"),
                        zaxis = list(title = "genome-wide C-hyper")),
           title = "3D Plot of Transect vs Transect C-hyper vs mean-W1 and mean-W2")

# Create the plotly 3D scatter plot for transect vs all
plot_ly(data = all_results, 
        x = ~W1, y = ~W2, z = ~Score, 
        type = "scatter3d", mode = "markers",
                  text = ~paste(Vector1, Vector2),
                  marker = list(symbol = "circle", size = 5)) %>%
    layout(scene = list(xaxis = list(title = "mean-W1"),
                        yaxis = list(title = "mean-W2"),
                        zaxis = list(title = "genome-wide C-hyper")),
           title = "3D Plot of Transect vs All C-hyper vs mean-W1 and mean-W2")





















#################3D Plots by Chrom -- read c scores by crhom first --> trans_list and all_list.  C scores and cohen Ws per chromosome
### FORMATTING
trans_list_backup <- trans_list
all_list_backup <- all_list

for (i in 1:nrow(cohen_df)) {

  row <- cohen_df[i, ]
  row <-as.vector(row)
  row <- unlist(row)  
  trans_list[[i]]$W1 <- row[trans_list[[i]]$Vector1]
  trans_list[[i]]$W2 <- row[trans_list[[i]]$Vector2]

}


for (i in 1:nrow(cohen_df)) {

  row <- cohen_df[i, ]
  row <-as.vector(row)
  row <- unlist(row)  
  all_list[[i]]$W1 <- row[all_list[[i]]$Vector1]
  all_list[[i]]$W2 <- row[all_list[[i]]$Vector2]
 

}




########### PLOTS START HERE

#########################  Plot transect results in 3D by Chrom - no repeated values
library(plotly)

# Function to create the plot for each dataframe
create_plot <- function(df, name) {
  plot <- plot_ly(data = df, 
                  x = ~W1, y = ~W2, z = ~Score,
                  type = "scatter3d", mode = "markers",
                  text = ~paste(Vector1, Vector2),
                  marker = list(symbol = "circle", size = 5)) %>%
    layout(scene = list(xaxis = list(title = "W1"),
                        yaxis = list(title = "W2"),
                        zaxis = list(title = "C-hyper")),
           title = paste("3D Plot of C-hyper vs W1 and W2 for chromosome", name))
  
  return(plot)
}


# Save Plots
for (i in seq_along(trans_list)) {
    plot <- create_plot(trans_list[[i]], i)
    plot_name <- paste("trans_vs_trans_", i, ".html", sep = "")
    htmlwidgets::saveWidget(plot, plot_name)
}





################################# Plot in 3D but all chromosomes together and including both lower and upper diagonal of the data
library(plotly)

# Combine all dataframes into one dataframe and add a categorical variable
combined_df <- do.call(rbind, Map(cbind, trans_list, name = seq_along(trans_list)))
combined_df$ID <- paste(combined_df$Vector1,combined_df$Vector2)



# Add duplicate points to plot 

# Duplicate rows and invert values of Vector1, Vector2, W1, and W2
inverted_df <- combined_df
inverted_df[c("Vector1", "Vector2")] <- combined_df[c("Vector2", "Vector1")]
inverted_df[c("W1", "W2")] <- combined_df[c("W2", "W1")]

# Combine original and inverted dataframes
final_df <- rbind(combined_df, inverted_df)

# Reset row names
rownames(final_df) <- NULL

# Show the resulting dataframe
head(final_df)


# Function to create the combined plot 3D
create_combined_plot <- function(df) {
  plot <- plot_ly(data = df, 
                  x = ~W1, y = ~W2, z = ~Score,
                  type = "scatter3d", mode = "markers",
                  text = ~paste(Vector1, Vector2),
                  color = ~as.factor(ID), # Set color based on the categorical variable
                  colors = "Set1", # Set color scheme
                  marker = list(symbol = "circle", size = 5)) %>%
    layout(scene = list(xaxis = list(title = "W1"),
                        yaxis = list(title = "W2"),
                        zaxis = list(title = "C-hyper")),
           title = "Combined 3D Plot of C-hyper vs W1 and W2 for all chromosomes")
  
  return(plot)
}

# Create the combined plot
combined_plot <- create_combined_plot(final_df)
print(combined_plot)




######################  By chrmo 3D plot with upper and lower diagonal
#######  Plot transect results in 3D
library(plotly)

# Function to create the plot for each dataframe
create_plot <- function(df, name) {
  plot <- plot_ly(data = df, 
                  x = ~W1, y = ~W2, z = ~Score,
                  type = "scatter3d", mode = "markers",
                  text = ~paste(Vector1, Vector2),
                  marker = list(symbol = "circle", size = 5)) %>%
    layout(scene = list(xaxis = list(title = "W1"),
                        yaxis = list(title = "W2"),
                        zaxis = list(title = "C-hyper")),
           title = paste("3D Plot of C-hyper vs W1 and W2 for chromosome", name))
  
  return(plot)
}


# Save Plots
for (i in seq_along(trans_list)) {
    combined_df <- trans_list[[i]]
    combined_df$ID <- paste(combined_df$Vector1,combined_df$Vector2)
	# Duplicate rows and invert values of Vector1, Vector2, W1, and W2
    inverted_df <- combined_df
    inverted_df[c("Vector1", "Vector2")] <- combined_df[c("Vector2", "Vector1")]
    inverted_df[c("W1", "W2")] <- combined_df[c("W2", "W1")]

    # Combine original and inverted dataframes
     final_df <- rbind(combined_df, inverted_df)

    # Reset row names
    rownames(final_df) <- NULL
    plot <- create_plot(final_df, i)
    plot_name <- paste("duplicated_trans_vs_trans_", i, ".html", sep = "")
    htmlwidgets::saveWidget(plot, plot_name)
}





















############ Plot all results
############################### 3D Plot
create_plot <- function(df, name) {
  plot <- plot_ly(data = df, 
                  x = ~W1, y = ~W2, z = ~Score,
                  type = "scatter3d", mode = "markers",
                  text = ~paste(Vector1, Vector2),
                  marker = list(symbol = "circle", size = 5)) %>%
    layout(scene = list(xaxis = list(title = "W1"),
                        yaxis = list(title = "W2"),
                        zaxis = list(title = "C-hyper")),
           title = paste("3D Plot of C-hyper vs W1 and W2 for chromosome", name))
  
  return(plot)
}

# Save Plots
for (i in seq_along(all_list)) {
    plot <- create_plot(all_list[[i]], i)
    plot_name <- paste("trans_vs_all_", i, ".html", sep = "")
    htmlwidgets::saveWidget(plot, plot_name)
}



############################## 2D Plots omitting all(banff) cohen W
# Combine all dataframes into one -- all chromosomes together
all_data <- do.call(rbind, all_list)

# Create the plot
ggplot(all_data, aes(x = W1, y = Score, color = Vector1)) +
  geom_point() +
  labs(title = "C-hyper vs cW", x = "Trasect cW", y = "C-hyper") +
  theme_minimal()





################################# 2D plot omitting all cohen W, each chromosopme in a different quadrant
library(ggplot2)
library(purrr)

# Function to create plot for each entry in all_list
create_plot <- function(data, title) {
  ggplot(data, aes(x = W1, y = Score, color = Vector1)) +
    geom_point() +
    labs(title = title, x = "Transect cW", y = "C-hyper", color = "Transect") +
    theme_minimal() +
    theme(legend.key.size = unit(0.5, "lines"))  # Adjust the size as per your preference
}

# Create separate plots for each entry in all_list and arrange them in a grid
plots <- map2(all_list, names(all_list), ~ create_plot(.x, .y))

# Arrange the plots in a grid
gridExtra::grid.arrange(grobs = plots, ncol = 3)  # Adjust ncol as per your preference


