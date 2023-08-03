library(reshape)
library(ggplot2)
library(cowplot)

# Read the count data from the TSV file
count_df <- read.table("/Users/pimswart/tester_count_file.tsv", header = TRUE, sep = "\t")

# Specify the folder path where the coverage files are stored
coverage_folder <- "/Users/pimswart/coverage_blast_test/"

# Create a blank list to store the coverage plots
coverage_plots <- list()

# Iterate over each row in the count_df DataFrame
for (i in 1:nrow(count_df)) {
  # Extract the sseqid and count from each row
  sseqid <- count_df$sseqid[i]
  count <- count_df$count[i]
  
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

# Create the bar chart plot with horizontal bars
bar_chart <- ggplot(count_df, aes(x = reorder(sseqid, count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "sseqid", y = "Count") +
  coord_flip() +
  ggtitle("Count of sseqid")

# Combine the bar chart and coverage plots
combined_plot <- plot_grid(bar_chart, plotlist = coverage_plots, nrow = 1, align = "v", labels = "AUTO")

# Save the combined plot as a PNG file
output_file <- "/Users/pimswart/combined_plot.png"
ggsave(filename = output_file, plot = combined_plot, device = "png")