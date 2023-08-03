library(gplots)
library(dplyr)
library(stringr)

# Read the count data from the CSV file
data <- read.csv("/Users/pimswart/go_terms_counts.csv")
# Remove leading and trailing whitespaces from the GO.Term column
data$GO.Term <- trimws(data$GO.Term)

# Split the data into three data frames for BP, MF, and CC
bp_data <- subset(data, GO.Term.Type == "Biological Process", select = c("Propagation", "GO.Term", "Occurrences"))
mf_data <- subset(data, GO.Term.Type == "Molecular Function", select = c("Propagation", "GO.Term", "Occurrences"))
cc_data <- subset(data, GO.Term.Type == "Cellular Component", select = c("Propagation", "GO.Term", "Occurrences"))

# Function to create a heatmap plot given the data frame
create_heatmap <- function(data, title) {
  # Create a matrix for heatmap
  heatmap_matrix <- reshape2::dcast(data, GO.Term ~ Propagation, value.var = "Occurrences", fill = 0)
  
  # Set row names to be GO terms
  rownames(heatmap_matrix) <- heatmap_matrix$GO.Term
  heatmap_matrix <- heatmap_matrix[, -1]
  
  # Order the columns based on P1 to P10
  propagation_order <- paste0("P", 1:10)
  heatmap_matrix <- heatmap_matrix[, order(factor(colnames(heatmap_matrix), levels = propagation_order))]
  
  # Calculate the relative abundance per propagation (normalize by the total count of genes per propagation)
  relative_abundance <- t(t(heatmap_matrix) / colSums(heatmap_matrix))
  
  # Create the heatmap plot using heatmap.2 without log transformation
  # Adjust the margins and size as needed
  heatmap.2(relative_abundance, Rowv = FALSE, Colv = FALSE, dendrogram = "none", trace = "none", 
            col = colorRampPalette(c("blue", "white", "red"))(100), key = TRUE, keysize = 1.5, density.info = "none", 
            main = title, cexCol = 1, cexRow = 1.5, margins = c(10, 10), labCol = c(10, 3))
}

# Create the three heatmaps
bp <- create_heatmap(bp_data, "Biological Process Heatmap (Relative Abundance)")
mf <- create_heatmap(mf_data, "Molecular Function Heatmap (Relative Abundance)")
cc <- create_heatmap(cc_data, "Cellular Component Heatmap (Relative Abundance)")