# Load the required packages
library(ggplot2)
library(dplyr)
library(scales)
# Set the file names
lengths_file <- "/Users/pimswart/P_all_pharokka_of_representatives_manualremovedbacterial_fakelengths.csv"
heatmap_file <- "/Users/pimswart/P_all_pharokka_of_representatives_manualremovedbacterial_fasta_sigs_compare.csv"
# IMPORTANT: THIS ORDERED CSV IS 'LEFT TO RIGHT' EQUAL TO THE ORDER OF THE SOURMASH PLOT
# THIS MEANS THE FIRST/LEFT SEQUENCE NAME IN THE CSV IS THE FIRST/TOP SEQUENCE NAME IN THE SOURMASH PLOT 

# Read in the sequence lengths file
lengths_df <- read.delim(lengths_file, header = FALSE, col.names = c("sequence_name", "length"))
# Read in the heatmap file
heatmap_df <- read.csv(heatmap_file)

# Extract the sequence names from the column headers of the heatmap file
seq_names <- colnames(heatmap_df)
#seq_names <- sub("^[^.]+\\.", "", colnames(heatmap_df))
#is.list(seq_names)
seq_names <- sub("\\.fa$", "", seq_names)
seq_names <- gsub("\\.", "\\-", seq_names)
seq_names <- sub("\\-", ".", seq_names)
seq_names
######## LENGTH PLOT
# Find the lengths corresponding to the sequence names and add them to the heatmap data frame
heatmap_df$length <- lengths_df$length[match(seq_names, lengths_df$sequence_name)]
heatmap_df$length
# Create a data frame with the sequence names and lengths
lengths_plot <- data.frame(seq_names = seq_names, length = heatmap_df$length)
lengths_plot <- lengths_plot
# Create the plot (THIS IS ORIENTED CORRECTLY)
seq_lengths_plot <- ggplot(lengths_plot, aes(x = length, y = seq_names)) +
  geom_point(size = 0.1) +
  scale_x_continuous(expand = c(0, 0), limits = c(15000, max(lengths_plot$length+10000)), labels = label_number(suffix="K",scale=0.001)) +
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
ggsave(filename='test',plot=seq_lengths_plot, height=10,width=1.5,units='cm',dpi=1000)


######## PROP PLOT
prop_df <- lengths_plot
# Add new column 'prop' containing the first part of the seq_names column
prop_df$prop <- sub("\\..*", "", prop_df$seq_names)
# Drop the length column
prop_df <- prop_df[, -2]

#TESTING
#prop_df$prop[1:100] <- "P10"

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
plot <- ggplot(prop_df_plotting, aes(x = seq_names, y = 1)) + 
  geom_tile(aes(fill = color), width = 1, height = 1) +
  scale_fill_identity() + coord_flip() + theme(aspect.ratio=1.5)
  #theme_void()
plot
# Create makeshift legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#000041", "#0000A9", "#0049FF", "#00A4DE", "#03E070", "#5DFC21", "#F6F905", "#FFD701", "#FF9500", "#FF3300"))
mtext("Propagations", at=0.25, cex=2)





######## MAPPING PLOT
best_mappings_file = "/Users/pimswart/2_extracted_best_mappings_minimap_P_all_representatives_manualremovedbacterial_fasta_catted_cat_references_namecharactersfixedforcircoletto_allprophageregionsmaskednotjustregion2_Ecolibl21DE3_LE392_and_maskedprophageregions_and_parentalphages.txt"
best_mappings_df <- read.csv(best_mappings_file, header = FALSE, sep='\t', col.names = c("seq_names", "best_mapping"))
# Reorder best_mappings_df to 'original' order of lengths_plot df
# Reorder rows of best_mappings_df based on seq_names in lengths_plot
best_mappings_df <- best_mappings_df[match(lengths_plot$seq_names, best_mappings_df$seq_names), ]

# Print the reordered dataframe
print(reordered_best_mappings)
# Define the base colors
color1 <- "#1f77b4"
color2 <- "#ff7f0e"
color3 <- "#2ca02c"
color4 <- "#d62728"
color5 <- "#8c564b"

# Generate the fading versions of color4
n_fade4 <- 6
fade_colors4 <- colorRampPalette(c(color4, "#FFFFFF"))(n_fade4 + 1)

# Generate the fading versions of color5
n_fade5 <- 11
fade_colors5 <- colorRampPalette(c(color5, "#FFFFFF"))(n_fade5 + 1)

# Combine the colors for the five categories
color_palette_2 <- c(color1, color2, color3, fade_colors4, fade_colors5)

# Create a named vector of the colors
colors_2 <- setNames(color_palette_2, c("S155L","S255s1","S236s","Ecolibl21DE3","Ecolibl21DEppregion1","Ecolibl21DE3ppregion2","Ecolibl21DE3ppregion3","Ecolibl21DE3ppregion4","Ecolibl21DE3ppregion5","Ecolibl21DE3ppregion6","LE392","LE392ppregion1","LE392ppregion2","LE392ppregion3","LE392ppregion4","LE392ppregion5","LE392ppregion6","LE392ppregion7","LE392ppregion8","LE392ppregion9","LE392ppregion10","LE392ppregion11"))

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
legend("topleft", legend =c("S155L","S255s1","S236s","Ecolibl21DE3","Ecolibl21DEppregion1","Ecolibl21DE3ppregion2","Ecolibl21DE3ppregion3","Ecolibl21DE3ppregion4","Ecolibl21DE3ppregion5","Ecolibl21DE3ppregion6","LE392","LE392ppregion1","LE392ppregion2","LE392ppregion3","LE392ppregion4","LE392ppregion5","LE392ppregion6","LE392ppregion7","LE392ppregion8","LE392ppregion9","LE392ppregion10","LE392ppregion11"), 
       pch=16, pt.cex=3, cex=1.5, bty='n',
       col = color_palette_2)
mtext("Reference genome", at=0.25, cex=2)

# Create makeshift legend
color_palette_2
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("S155L","S255s1","S236s","Ecolibl21DE3","Ecolibl21DEppregion1","Ecolibl21DE3ppregion2","Ecolibl21DE3ppregion3","Ecolibl21DE3ppregion4","Ecolibl21DE3ppregion5","Ecolibl21DE3ppregion6","LE392","LE392ppregion1","LE392ppregion2","LE392ppregion3","LE392ppregion4","LE392ppregion5","LE392ppregion6","LE392ppregion7","LE392ppregion8","LE392ppregion9","LE392ppregion10","LE392ppregion11"), 
       pch=16, pt.cex=2, cex=1, bty='n',
       col = color_palette_2)
mtext("Legend B. Propagation", at=0.25, cex=2)



colors_2
fade_colors4
fade_colors5


# Plot length distribution histogram for reads and assemblies
library(Biostrings)
fasta_data <- readDNAStringSet("/Users/pimswart/assembly_allreads_regularsettings_sept2021.contigs.fasta")
# Extract read lengths
read_lengths <- width(fasta_data)
bin_size <- 1000
# Create a data frame with the read lengths
data <- data.frame(Lengths = read_lengths)

# Plot the histogram with axis break
ggplot(data, aes(x = Lengths)) +
  geom_histogram(binwidth = bin_size, boundary = 0, fill = "steelblue", color = "white") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  labs(x = "Read Length", y = "Frequency") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))
