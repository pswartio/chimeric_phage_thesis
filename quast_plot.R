library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(tibble)
library(tidyverse)
library(ggbreak) 
library(scales)


# Get the list of files in the directory
file_list <- list.files(pattern = "metaquast_p\\d+_withref_report kopie.tsv")
file_list

# Create a list to store the dataframes for each file
dataframes <- list()

# Iterate over each file and read it as a dataframe
for (file in file_list) {
  df <- read.delim(file, sep = "\t", header = TRUE, row.names = 1)
  dataframes[[file]] <- df
}

# Define the statistics of interest
statistics <- c("# contigs", "Total length", "# misassemblies", "Unaligned length", "Total aligned length", "Duplication ratio", "Genome fraction (%)")

# Create an empty dataframe to store the combined values
combined_df <- data.frame(assembler = character(),
                          statistic = character(),
                          value = numeric(),
                          stringsAsFactors = FALSE)

# Combine the values from each dataframe into a single dataframe
for (i in 1:length(dataframes)) {
  df <- dataframes[[i]]
  print(df)
  for (statistic in statistics) {
    assembler_values <- unlist(df[statistic, ])
    assembler <- colnames(df)
    combined_df <- rbind(combined_df, data.frame(assembler = assembler, statistic = statistic, value = assembler_values))
  }
}

# Convert the value column to numeric
combined_df$value <- as.numeric(combined_df$value)

# Calculate the average for each statistic and assembler
average_df <- combined_df %>%
  group_by(assembler, statistic) %>%
  summarize(average_value = mean(value, na.rm = TRUE))

# Convert the assembler column to a factor to control the order in the plot
average_df$assembler <- factor(average_df$assembler, levels = unique(average_df$assembler))

# Filter out the "Canu_regular_old" assembler from the dataframe
average_df <- average_df %>%
  filter(assembler != "Canu_regular_old")

# Convert the 'assembler' column to character
average_df <- average_df %>%
  mutate(assembler = as.character(assembler))

# Rename the data for "Canu_regular_new" to "Canu_regular" in the dataframe
average_df <- average_df %>%
  mutate(assembler = ifelse(assembler == "Canu_regular_new", "Canu_regular", assembler))

# Remove decimal places in the 'average_value' column of average_df
average_df$average_value <- round(average_df$average_value, 0)

plot_width <- 15
plot_height <- 14

# Create the bar chart plot
plot <- ggplot(average_df, aes(x = assembler, y = average_value, fill = statistic)) +
  geom_col(position = position_dodge()) +
  labs(x = "Assembler", y = "Average Value", fill = "Statistic") +
  ggtitle("Average Statistics for Assemblies by Assembler") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_break(c(10000,50000), scales=1.5)  + 
  scale_y_break(c(10000000,70000000), scales=1.5) 

plot

# Save the resized plot as an image file
ggsave("resized_plot.png", plot, width = plot_width, height = plot_height)