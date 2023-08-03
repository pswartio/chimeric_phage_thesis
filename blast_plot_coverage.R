library(reshape)
library(ggplot2)
library(dplyr)
library(magrittr)

# Read the count data from the TSV file
count_df <- read.table("/Users/pimswart/count_file_manuallyadded kopie.tsv", header = TRUE, sep = "\t")

# Specify the folder path where the coverage files are stored
coverage_folder <- "/Users/pimswart/coverage_blast_test/"

# Create a blank list to store the coverage plots
coverage_plots <- list()

# Iterate over each row in the count_df DataFrame
for (i in 1:nrow(count_df)) {
  # Extract the sseqid and count from each row
  sseqid <- count_df$sseqid[i]
  count <- count_df$count[i]
  scinames <- count_df$scinames[i]
  
  # Construct the coverage file path
  coverage_file <- paste0(coverage_folder, sseqid, ".txt")
  
  # Read the coverage data from the file
  coverage <- read.table(coverage_file, sep = "\t", header = FALSE)
  
  # Rename the columns of the coverage data
  colnames(coverage) <- c("locus", "depth")
  
  # Create the ggplot coverage plot
  p <- ggplot(coverage, aes(x = locus, y = depth)) +
    geom_point(colour = "blue", size = 1, shape = 20, alpha = 1/3) +
    labs(y = "Depth", x = "Base position") +
    ggtitle(paste("Coverage of P", count, "assemblies per basepair of", sseqid))
  
  # Add the coverage plot to the list
  coverage_plots[[sseqid]] <- p
}

#REMOVE
# Plot and save the coverage plots
for (sseqid in names(coverage_plots)) {
  # Print the plot
  print(coverage_plots[[sseqid]] + labs(y = "Depth", x = "Base position"))

  # Save the plot as a PNG file
  filename <- paste0("/Users/pimswart/", sseqid, "_coverage_plot.png")
  ggsave(filename = filename, plot = coverage_plots[[sseqid]], device = "png", width = 8, height = 4)
}
#REMOVE



# Sort the filtered count_df by count in descending order
ordered_count_df <- count_df %>%
  arrange(desc(count))

# Find the index of the row where count is 1 for the second time
index <- min(which(cumsum(ordered_count_df$count == 1) == 2))

# Remove the lines after the index
ordered_count_df <- ordered_count_df[seq_len(index), ]

# Create a new column with duplicated sciname and a newline separator
ordered_count_df$label <- paste(ordered_count_df$sciname, ordered_count_df$sseqid, sep = "\n")

# Print the resulting dataframe
print(ordered_count_df)

# Create the horizontal bar chart plot with sorted bars
bar_chart <- ggplot(ordered_count_df, aes(x = count, y = reorder(label, count))) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.8, position = position_dodge(width = 0.5)) +
  labs(x = "sseqid", y = "Count") +
  ggtitle("Count of sseqid") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 7))

bar_chart
# Save the bar chart plot as a PDF file
output_file <- "/Users/pimswart/bar_chart.png"
ggsave(filename = output_file, plot = bar_chart, device = "png", width = 8, height = 15)

