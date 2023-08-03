# Load the required packages
library(ggplot2)
library(dplyr)
library(scales)
# Set the file names
lengths_file <- "/Users/pimswart/Downloads/lengths_sourmash_compare_assembly_allreads_regularsettings_sept2021_k31"
heatmap_file <- "/Users/pimswart/Downloads/sourmash_compare_assembly_allreads_regularsettings_sept2021_k31_dna_ordered.csv"
# IMPORTANT: THIS ORDERED CSV IS 'LEFT TO RIGHT' EQUAL TO THE ORDER OF THE SOURMASH PLOT
# THIS MEANS THE FIRST/LEFT SEQUENCE NAME IN THE CSV IS THE FIRST/TOP SEQUENCE NAME IN THE SOURMASH PLOT 

# Read in the sequence lengths file
lengths_df <- read.delim(lengths_file, header = FALSE, col.names = c("sequence_name", "length"))
# Read in the heatmap file
heatmap_df <- read.csv(heatmap_file)

# Extract the sequence names from the column headers of the heatmap file
seq_names <- sub("\\..*", "", colnames(heatmap_df ))
seq_names

######## LENGTH PLOT
# Find the lengths corresponding to the sequence names and add them to the heatmap data frame
heatmap_df$length <- lengths_df$length[match(seq_names, lengths_df$sequence_name)]
heatmap_df$length
# Create a data frame with the sequence names and lengths
lengths_plot <- data.frame(seq_namess = seq_names, length = heatmap_df$length)
lengths_plot$seq_names <- factor(lengths_plot$seq_names, levels = unique(lengths_plot$seq_names))

# Create the plot (THIS IS ORIENTED CORRECTLY)
seq_lengths_plot <- ggplot(lengths_plot, aes(x = length, y = seq_names)) +
  geom_point(size = 0.1) +
  scale_x_log10(expand = c(0, 0), limits = c(2000, max(lengths_plot$length+50000)), labels = label_number(suffix="K",scale=0.001)) +
  xlab("") +
  #ylab("Sequence ID") +
  theme(axis.ticks.y = element_blank(),
        #axis.text.y = element_text(size = 8, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(),
        plot.margin = margin(0,0,0,0),
        axis.text.x = element_text(angle=270,hjust=1,vjust=0.5,size=3))
#Display plot
seq_lengths_plot
ggsave(filename='primitivepool_additional_sourmash_protein.png',plot=seq_lengths_plot, height=10,width=1.5,units='cm',dpi=1000)


######## PROP PLOT
prop_df <- lengths_plot
# Add new column 'prop' containing the first part of the seq_names column
prop_df$prop <- sub("\\..*", "", prop_df$seq_names)
# Drop the length column
prop_df <- prop_df[, -2]

#TESTING
prop_df$prop[1:100] <- "P10"

# create a named vector of colors
color_palette <- c("#000041", "#0000A9", "#0049FF", "#00A4DE", "#03E070", "#5DFC21", "#F6F905", "#FFD701", "#FF9500", "#FF3300")
# Create a named vector of the colors
colors <- setNames(color_palette, paste0("P", 1:10))
# Create a duplicate of the prop_df dataframe for plotting
prop_df_plotting <- prop_df

# For testing purposes: sample the dataframe
#prop_df_plotting <- tail(prop_df,n=50)

# Map colors to prop column in prop_df
prop_df_plotting$color <- colors[prop_df_plotting$prop]
prop_df_plotting$seq_names <- paste0("Row", seq_len(nrow(prop_df_plotting)))

prop_df_plotting$seq_names <- factor(prop_df_plotting$seq_names, levels=unique(prop_df_plotting$seq_names))

# Create plot (COORD_FLIP ROTATES IT TO PLOT LINES HORIZONTALLY, BUT DOES PLOT UPSIDE DOWN!!!!)
ggplot(prop_df_plotting, aes(x = seq_names, y = 1)) + 
  geom_tile(aes(fill = color), width = 1, height = 1) +
  scale_fill_identity() + coord_flip() + theme(aspect.ratio=1.5)
  #theme_void()

# Create makeshift legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#000041", "#0000A9", "#0049FF", "#00A4DE", "#03E070", "#5DFC21", "#F6F905", "#FFD701", "#FF9500", "#FF3300"))
mtext("Propagations", at=0.25, cex=2)





######## MAPPING PLOT
best_mappings_file = "/Users/pimswart/2_extracted_best_mappings_minimap_assembly_allreads_regularsettings_sept2021.contigs_cat_references_Ecolibl21DE_LE392_S155L_S236s_S255s1.txt"
best_mappings_df <- read.csv(best_mappings_file, header = FALSE, sep='\t', col.names = c("seq_names", "best_mapping"))
# Reorder best_mappings_df to 'original' order of lengths_plot df
# Reorder rows of best_mappings_df based on seq_names in lengths_plot
best_mappings_df <- best_mappings_df[match(lengths_plot$seq_names, best_mappings_df$seq_names), ]

# Print the reordered dataframe
#print(reordered_best_mappings)
# Define the base colors
color1 <- "#1f77b4"
color2 <- "#ff7f0e"
color3 <- "#2ca02c"
color4 <- "#d62728"
color5 <- "#8c564b"

color_palette_2 <- c(color1, color2, color3, color4, color5)

# Create a named vector of the colors
colors_2 <- setNames(color_palette_2, c("S1_55L","S2_55s1","S2_36s","Ecolibl21DE3","LE392"))

# Create a duplicate of the best_mappings_df dataframe for plotting
best_mappings_df_plotting <- best_mappings_df

# Map colors to prop column in prop_df
best_mappings_df_plotting$color <- colors_2[best_mappings_df_plotting$best_mapping]
best_mappings_df_plotting$seq_names <- paste0("Row", seq_len(nrow(best_mappings_df_plotting)))

best_mappings_df_plotting$seq_names <- factor(best_mappings_df_plotting$seq_names, levels=unique(best_mappings_df_plotting$seq_names))

# Create plot (COORD_FLIP ROTATES IT TO PLOT LINES HORIZONTALLY, BUT DOES PLOT UPSIDE DOWN!!!!)
ggplot(best_mappings_df_plotting, aes(x = seq_names, y = 1)) + 
  geom_tile(aes(fill = color), width = 1, height = 1) +
  scale_fill_identity() + coord_flip() + theme(aspect.ratio=1.5 )
#theme_void()


# Create makeshift legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("S155L","S255s1","S236s","Ecolibl21DE3","LE392"), 
       pch=16, pt.cex=3, cex=1.5, bty='n',
       col = color_palette_2)
mtext("Legend B. Best mapping", at=0.25, cex=2)

# Create makeshift legend
color_palette_2
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("S155L","S255s1","S236s","Ecolibl21DE3","LE392"), 
       pch=16, pt.cex=2, cex=1, bty='n',
       col = color_palette_2)
mtext("Legend A/D. Heatmap scale", at=0.5, cex=2)



colors_2
fade_colors4
fade_colors5
