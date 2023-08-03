library(tidyverse)
library(ewrefxn)
#library(sparseMatrixUtils)
#library(Matrix)
library(readsparse)
library(ape)
library(phangorn)
library(slam)
#library(futile.matrix)
library(tidyr)
library(dplyr)
library(phylobase)
library(phytools)
#sessionInfo()

#MY WORKING CODE
#reading file into df, filtering out duplicates
a <- read.table("/Users/pimswart/Downloads/odgipaths_test_notgrouped", sep="\t", header = TRUE) #maybe remove header
naamgeving_df <- data.frame(path.a=a$path.a)
names(naamgeving_df)[names(naamgeving_df) == "path.a"] <- "Naam"
naamgeving_df <- naamgeving_df %>% distinct(Naam, .keep_all = TRUE)
#ordering the df
naamgeving_df$Naam <- naamgeving_df[order(naamgeving_df$Naam),]
#print(naamgeving_df)
#now add column to the dataframe, write the first characters before . in there (so P1.asdfb -- P1)
naamgeving_df$Propagation <- substr(naamgeving_df$Naam, start = 1, stop = 3)
naamgeving_df$Propagation <- gsub('[^[:alnum:] ]', '', naamgeving_df$Propagation)
#construction of tree matrix starts here
nameVals <- sort(unique(unlist(a[1:1])))
#print(nameVals)
#construct 0 matrix of correct dimensions with row and column names
myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
# fill in the matrix with matrix indexing on row and column names
#print(myMat)
myMat[as.matrix(a[c("path.a", "path.b")])] <- a[["euclidean"]]
#print(myMat)
#nj on the matrix
my_nj <- ape::nj(myMat)

##PLOTTING VERSION THREE (EDGE COLOURING) USING PHYLO
colIDs <- colnames(myMat)
names(colIDs) <- colIDs
colIDs
trait<-as.factor(colIDs)
trait
factor(substring(colIDs,1,2))
colIDs
levels <- substr(colIDs,1,2)
levels
colnames<- as.factor(colIDs)
colnames
TYPE = gsub("[.][^ ]*","",colIDs)
xfactor <- factor(TYPE)
xfactor


#nj on the matrix
my_nj <- ape::nj(myMat)
plot(my_nj)

p1 <- names(xfactor)[xfactor=="P1"]
p2 <- names(xfactor)[xfactor=="P2"]
p3 <- names(xfactor)[xfactor=="P3"]
p4 <- names(xfactor)[xfactor=="P4"]
p5 <- names(xfactor)[xfactor=="P5"]
p6 <- names(xfactor)[xfactor=="P6"]
p7 <- names(xfactor)[xfactor=="P7"]
p8 <- names(xfactor)[xfactor=="P8"]
p9 <- names(xfactor)[xfactor=="P9"]
p10 <- names(xfactor)[xfactor=="P10"]
pdoesnexist <- names(xfactor)[xfactor=="doesntexist"]

#maybe just copy this logic for all P's?
tt<-paintBranches(my_nj,edge=sapply(p1,match,my_nj$tip.label), state="P1",anc="doesntexist")
tt<-paintBranches(tt,edge=sapply(p2,match,my_nj$tip.label), state="P2")
tt<-paintBranches(tt,edge=sapply(p3,match,my_nj$tip.label), state="P3")
tt<-paintBranches(tt,edge=sapply(p4,match,my_nj$tip.label), state="P4")
tt<-paintBranches(tt,edge=sapply(p5,match,my_nj$tip.label), state="P5")
tt<-paintBranches(tt,edge=sapply(p6,match,my_nj$tip.label), state="P6")
tt<-paintBranches(tt,edge=sapply(p7,match,my_nj$tip.label), state="P7")
tt<-paintBranches(tt,edge=sapply(p8,match,my_nj$tip.label), state="P8")
tt<-paintBranches(tt,edge=sapply(p9,match,my_nj$tip.label), state="P9")
tt<-paintBranches(tt,edge=sapply(p10,match,my_nj$tip.label), state="P10")

cols<-setNames(c("#FFFFCC", "#FFFF99", "#FFFF66", "#FFFF33", "#FFFF00", "#FFCC00", "#FF9900", "#FF6600", "#FF3300", "#FF0000","grey"),c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","doesntexist"))
png(file="odgi_output_sample.png",res=200,width=3000,height=3000)
plot(tt,type="fan",colors=cols,lwd=0.5,split.vertical=TRUE,ftype="i",fsize=0.2)
dev.off()
#can change ftype to "off" to remove labels, otherwise keep as "i"



##PLOTTING VERSION ONE (TIPLABEL COLOURING) AND NOT-WORKING VERSION TWO) (EDGE COLOURING) USING GGTREE AND APE
rich10equal = c("#000041", "#0000A9", "#0049FF", "#00A4DE", "#03E070", "#5DFC21", "#F6F905", "#FFD701", "#FF9500", "#FF3300")
#https://stackoverflow.com/questions/58554889/in-r-how-can-i-color-the-labels-from-my-phylogenetic-tree-using-bionj-from-ape
TYPE = gsub("[.][^ ]*","",my_nj$tip.label)
col_assignment = rich10equal[1:length(unique(TYPE))]
names( col_assignment) = unique(TYPE)
COLS = col_assignment[TYPE]
TYPE
col_assignment
COLS
#plotting that does work: tip.color, colours correctly
ape::plot.phylo(unroot(my_nj), cex=0.3, tip.color=COLS)
#plotting that doesn't work:
ape::plot.phylo(unroot(my_nj), cex=0.3, edge.color=COLS)
#plotting that doesn't work with the above way of assigning colours from COLS:
ggtree(my_nj, layout="unrooted")
#Old plotting options that don't work
#ggtree(my_nj) %<+% data.frame(node=1:nrow(my_nj$edge), group.name=factor(c(as.character(a$group.name),rep("internal",nrow(y$edge)-nrow(y))), levels=c(levels(y$group.name),"internal") )) + aes(color=group.name) + geom_tree() + scale_color_manual("passage",values=c(phage.colors, "black"))
#ggtree(my_nj) %<+% data.frame(node=1:nrow(my_nj$edge), group.name=factor(c(as.character(myMat$group.name),rep("internal",nrow(my_nj$edge)-nrow(y))), levels=c(levels(myMat$group.name),"internal") )) + aes(color=group.name) + geom_tree() + scale_color_manual("passage",values=c(phage.colors, "black"))
#ggsave(paste(output, "ggtree.passage.pdf", sep="."), height=40, width=9)




#mine+eriks code
library(tidyverse)
library(ewrefxn)
#library(sparseMatrixUtils)
#library(Matrix)
library(readsparse)
library(ape)
library(phangorn)
library(slam)
#library(futile.matrix)
library(tidyr)
library(dplyr)
library(phylobase)
library(phytools)
#sessionInfo()

a <- read.table("/Users/pimswart/Downloads/odgipaths_test_grouped", sep="\t", header = TRUE) #maybe remove header
aa <- as.matrix(a)
aaa<-unique(aa[c('path.a','path.b')])
#print(aa)
typeof(aa)
nameVals <- sort(unique(unlist(a[1:1])))
print(nameVals)
# construct 0 matrix of correct dimensions with row and column names
myMat <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
# fill in the matrix with matrix indexing on row and column names
print(myMat)
myMat[as.matrix(a[c("group.a", "group.b")])] <- a[["euclidean"]]
print(myMat)
my_nj <- ape::nj(myMat)
plot(my_nj,cex=0.5)

#PCA stuff
colfunc <- colorRampPalette(c("red", "yellow"))
phage.colors=c(colfunc(10), rainbow(8)[3:7])
y.pca <- prcomp(myMat)
y.pca.df <- as.data.frame(y.pca$x)
y.pca.df$group.name <- y$group.name
ggplot(y.pca.df, aes(x=PC1, y=PC2, color=group.name)) + geom_point() + scale_color_manual("passage",values=phage.colors)
ggsave(paste(output, "pca.PC1.PC2.pdf", sep="."), height=8, width=9)
ggplot(y.pca.df, aes(x=PC2, y=PC3, color=group.name)) + geom_point() + scale_color_manual("passage",values=phage.colors)
ggsave(paste(output, "pca.PC2.PC3.pdf", sep="."), height=8, width=9)
ggplot(y.pca.df, aes(x=PC3, y=PC4, color=group.name)) + geom_point() + scale_color_manual("passage",values=phage.colors)
ggsave(paste(output, "pca.PC3.PC4.pdf", sep="."), height=8, width=9)
ggplot(y.pca.df, aes(x=PC4, y=PC5, color=group.name)) + geom_point() + scale_color_manual("passage",values=phage.colors)
ggsave(paste(output, "pca.PC4.PC5.pdf", sep="."), height=8, width=9)

