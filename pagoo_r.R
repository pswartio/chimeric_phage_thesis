devtools::install_github("iferres/pagoo")
library(pagoo)
# Load Roary output and Prokka .gff output
gffs <- list.files(path = "/Users/pimswart/Downloads/P_all_collected_gff_new_prokka_notcorrupt_removed", pattern = "[.]gff$", full.names = TRUE)
gpa_csv <- "/Users/pimswart/ontvangstdir/gene_presence_absence.csv"
#p <- roary_2_pagoo(gene_presence_absence_csv = gpa_csv, gffs = gffs)

#orgs_file <- "/Users/macbookpro/Documents/uni/Msc/P18.tsv"
#orgs_meta <- read.table(orgs_file, header = TRUE, sep = "\t", quote = "")
#head(orgs_meta)
library(pagoo)
# Load the Roary output into Pagoo object
p <- roary_2_pagoo(gene_presence_absence_csv = gpa_csv, gffs = gffs)

# Convert the column with Propagation numbering into a vector to feed into dataframe to add to Pagoo object
orgs_file <- "/Users/pimswart/Downloads/P_all.tsv"
data=read.csv(orgs_file,sep="\t",header=TRUE)
data=data$propagation
vectordata=as.vector(data)
length(vectordata)
length(p$organisms$org)
propagation_df <- data.frame(org=p$organisms$org, propagation=vectordata)
p$add_metadata(map="org", propagation_df)
p$organisms
table(p$organisms$propagation)

colnames(propagation_df)

p$summary_stats
p$core_genes
p$core_clusters
p$pan_matrix[1:10, 1:10]
p$clusters
pca <- p$pan_pca()

p$gg_pca(colour = "propagation", size = 4) + theme_bw(base_size = 15) + scale_color_brewer(palette = "Set3")
ggsave(filename="P_all_PCA",device=png)
p$gg_curves(size = 2) + ggtitle("Pangenome curves") + geom_point(alpha = 0.1, size = 4) + theme_bw(base_size = 15) + ylim(0, 5000) + scale_color_brewer(palette = "Accent")
ggsave(filename="P_all_Pangenomecurves",device=png)
p$gg_barplot() + theme_bw(base_size = 15)+ theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12)) + geom_bar(stat = "identity", color = "black", fill = "black")
ggsave(filename="P_all_Barchart",device=png)
p$gg_dist() + theme_bw(base_size = 15) #maybe remove this theme thing if it errors
ggsave(filename="P_all_distanceheatmap",device=png)
p$gg_pie() + theme_bw(base_size = 15) + scale_fill_brewer(palette = "Blues") + scale_x_discrete(breaks = c(0, 25, 50, 75)) + theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 10), legend.margin = margin(0, 0, 13, 0), legend.box.margin = margin(0, 0, 5, 0), axis.title = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank())
ggsave(filename="P_all_piechart",device=png)