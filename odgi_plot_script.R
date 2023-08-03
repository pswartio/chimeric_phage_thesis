#http://blog.phytools.org/2017/01/plotting-terminal-edges-of-tree.html
library(phytools)

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

# Filter branches with negative length
#tt$edge.length[tt$edge.length < 0] <- 0

cols<-setNames(c("#000041", "#0000A9", "#0049FF", "#00A4DE", "#03E070", "#5DFC21", "#F6F905", "#FFD701", "#FF9500", "#FF3300","grey"),c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","doesntexist"))
plot(tt,type="fan",colors=cols,lwd=1,split.vertical=TRUE,ftype="i",fsize=0.2)
