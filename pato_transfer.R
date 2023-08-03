library(pato)
library(tidyverse)
library(ggtree)
library(ewrefxn)
library(tidyr)
library(ape)
library(Matrix)
library(dplyr)
find.package('pato')
sessionInfo()
getwd()
setwd("/Users/pimswart/Downloads/P_all_collected_gff_new_prokka_notcorrupt")

gff_files <- dir("/Users/pimswart/Downloads/P_all_collected_gff_new_prokka_notcorrupt", pattern = ".gff", full.names = T)
gffs <- load_gff_list(gff_files)
print("kanker")
strain_names <- read.table("/Users/pimswart/Downloads/namess.txt")%>%
  #Rename the columns
  rename(Genome = V1, Name = V2) %>%
  #delete the '-' character
  #mutate(Name = gsub("-","",Name)) %>%
  #Extract the first 4 character as Sample name
  #mutate(Sample = str_sub(Name,1,4)) %>%
  #Extract the final characters as Source
  mutate(Source = gsub("[.][^ ]*","",Name)) %>%y
  #remove the final part of the filename
  mutate(Genome = str_replace(Genome,".gff",""))

mash <- mash(gffs, type ="wgs", n_cores = 20)
#DOUBLE
my_mmseq <- mmseqs(gffs, type = "prot",cov_mode=2,coverage=0.8,identity=0.8,evalue=1e-6)
#my_accnet <- accnet(my_mmseq,threshold=0.001,singles=TRUE)
nr_list <- non_redundant(mash, distance = 0.01)
#nr_list
#typeof(nr_list)
representative_df = subset(nr_list, select = -c(Source) )
# Add the row names as a column to the data frame
representative_df$row_name <- row.names(representative_df)
# Group the data by "cluster"
representative_df <- split(representative_df, representative_df$cluster)
# Keep only the row with the highest "centrality" value for each group
representative_df <- lapply(representative_df, function(x) {
  x[x$centrality == max(x$centrality), ]
})
# Combine the groups back into a single data frame
representative_df <- do.call(rbind, representative_df)
# Set the row names to the values in the "row_name" column
row.names(representative_df) <- representative_df$row_name
# Remove the "row_name" column
representative_df <- representative_df[, -which(names(representative_df) == "row_name")]
# Get the row names from the dataframe
row_names <- rownames(representative_df)
# Remove the "_out.fna" string from the row names
row_names <- sub("_out.fna", "", row_names)
# Create a vector of file paths using the row names
file_paths <- paste0("/Users/pimswart/Downloads/P_all_collected_gff_new_prokka_notcorrupt/", row_names, "_out.gff")
# Bind the row names and file paths into a dataframe
output_df <- data.frame(row_names, file_paths)
# Write the output dataframe to a tab-separated file
write.table(output_df, file = "outputt_P_all_representatives.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)