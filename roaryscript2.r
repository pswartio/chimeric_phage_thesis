library(gplots)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
#/Users/pimswart/p_all_roary_pharokka_t1_i70_e/gene_presence_absence.csv
#/Users/pimswart/P_all_separate_representatives_pharokka_gffs_roary_pharokka_t1_i70_e/gene_presence_absence.csv
#/Users/pimswart/P_all_separate_representatives_manualremovedbacterial_pharokka_gffs_roary_pharokka_t1_i70_e/gene_presence_absence.csv
table_input <- read.csv("/Users/pimswart/P_all_separate_representatives_pharokka_gffs_roary_pharokka_t1_i70_e/gene_presence_absence.csv", sep=",", na.strings=c("","NA"))
table_input <- as.data.frame(table_input)
table_values <- within(table_input, rm("Gene","Non.unique.Gene.name","No..isolates","No..sequences","Avg.sequences.per.isolate","Genome.Fragment","Order.within.Fragment","Accessory.Fragment","Accessory.Order.with.Fragment","QC","Min.group.size.nuc","Max.group.size.nuc","Avg.group.size.nuc"))
abscence_presence <- as.matrix(table_values[,-1])
rownames(abscence_presence) <- table_values[,1]
abscence_presence[is.na(abscence_presence)] <- 0
abscence_presence[which(abscence_presence!=0)] <- 1
a_p_matrix <- mapply(abscence_presence, FUN=as.numeric)
a_p_matrix <- matrix(data=a_p_matrix, ncol=length(colnames(abscence_presence)), nrow=length(row.names(abscence_presence)))
row.names(a_p_matrix) <- row.names(abscence_presence)
colnames(a_p_matrix) <- colnames(abscence_presence)

# Define color palette for column labels
label_colors <- c("#000041", "#0000A9", "#0049FF", "#00A4DE", "#03E070", "#5DFC21", "#F6F905", "#FFD701", "#FF9500", "#FF3300")

# Create a function to retrieve the prefix number from a label
get_prefix_number <- function(label) {
  as.numeric(gsub("P", "", label))
}

# Sort the column labels based on prefix number
sorted_colnames <- colnames(a_p_matrix)[order(sapply(colnames(a_p_matrix), get_prefix_number))]

# Create a vector of label colors
col_label_colors <- sapply(colnames(a_p_matrix), function(label) {
  if (startsWith(label, "P10")) {
    "#FF0000"  # Change the color for labels starting with "P10"
  } else {
    label_colors[which(get_prefix(label) == paste0("P", 1:10))]
  }
})
# Plot the heatmap
pdf(file="coloured_roary_separate_representatives_pharokka_heatmap.pdf", width = 50, height = 50)
heatmap.2(a_p_matrix, col = c("#FFE986","#FF736E"), main = "Absence/Presence of genes",
          trace = "none", labRow = FALSE, cexRow = 1, cexCol = 1, margins = c(20, 20),
          colCol = col_label_colors)  # Set the column label colors
dev.off()
#plot







genomes_count <- length(colnames(a_p_matrix))
abscence_presence <- cbind(a_p_matrix, rowSums(a_p_matrix))
summary_table <- matrix(data=NA, nrow=3, ncol=length(colnames(abscence_presence)))
colnames(summary_table) <- colnames(abscence_presence)
rownames(summary_table) <- c("Total_genes","Unique_genes","Core_genes")
summary_table[1,] <- colSums(abscence_presence)
summary_table[2,] <- colSums(abscence_presence[which(abscence_presence[,ncol(abscence_presence)] == 1),])
summary_table[3,] <- colSums(abscence_presence[which(abscence_presence[,ncol(abscence_presence)] >= (genomes_count*0.95)),])
summary_table <- summary_table[,-ncol(summary_table)]
average_table <- data.frame(x=1:6, y=1:6, z=1:6)
average_table[,1] <- c("Total genes analyzed","Orthologous groups","Average gene count","Average core genes","Average unique genes","Total unique genes")
average_table[1,2] <- sum(summary_table[1,])
average_table[2,2] <- length(rownames(abscence_presence))
average_table[3,2] <- median(summary_table[1,])
average_table[4,2] <- length(rownames(abscence_presence[which(abscence_presence[,ncol(abscence_presence)] >= (genomes_count*0.95)),]))
average_table[5,2] <- round(length(rownames(abscence_presence[which(abscence_presence[,ncol(abscence_presence)] == 1),]))/length(colnames(abscence_presence)))
average_table[6,2] <- length(rownames(abscence_presence[which(abscence_presence[,ncol(abscence_presence)] == 1),]))
melt_summary_table <- melt(summary_table)
melt_summary_table <- melt_summary_table[order(melt_summary_table$value),]
p1 <- ggplot(melt_summary_table, aes(x = reorder(Var2, value), y = value)) + 
          geom_bar( stat = 'identity') + 
          facet_grid(. ~ Var1, scales = "free_x") + 
          xlab("Genomes") +
          ylab("Count") +
          coord_flip()
p2 <- ggplot(data=average_table[-c(1,2,6),], aes(x=x, y=y))+
          geom_bar(stat = 'identity') +
          theme (axis.text.x=element_text(angle=90,hjust=1,vjust=0.3),
          axis.title.x = element_blank()) +
          geom_text(aes(y = 10, label = paste("N =" ,y),vjust = 0), colour = "white", show.legend=FALSE) +
          ylab("Count")
t1  <- textGrob(paste(c("Total number of genomes:\n",
                        length(colnames(summary_table)),
                        "\n\nNumber of analyzed genes:\n",
                        as.numeric(average_table[1,2]),
                        "\n\nTotal orthologous groups\n",
                        as.numeric(average_table[2,2]),
                        "\n\nTotal unique genes\n",
                        as.numeric(average_table[6,2])), collapse = " ")) 
lay <- rbind(c(1,1,2),
             c(1,1,2),
             c(1,1,3),
             c(1,1,3))
grid.arrange(p1,p2,t1, layout_matrix = lay)
