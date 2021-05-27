
##
#  This first section is to make the xover dataframe from the crimap files. By including this the script as a whole can stand alone
#  rather than rely on another script to pull the data in.
##
rm(list = ls())
library(GenABEL)
library(crimaptools)
library(dplyr)
setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Scripts/Model")
load("sparrowABEL.RData")
sparrow.famped <- read.table("FamPed_20200414.txt", header = TRUE)
sparrow.famped <- sparrow.famped[which(sparrow.famped$ANIMAL %in% sparrow.abel@phdata$id),]
sparrow.famped <- sparrow.famped[sparrow.famped$Family %in% names(which(table(sparrow.famped$Family)>=5)),]
sparrow.famped$Family <- as.character(sparrow.famped$Family)
sparrow.famped$Family <- as.factor(sparrow.famped$Family)
setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Data/Crimap Runs/Final/crimap")
setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Data/Crimap Runs/Final/crimap")
#Make Maps
product <- list()
counter <- 0
for (i in (c(5:15, 17:28))) {
  
  sparrow.map <- parse_map(mapfile = paste("chr", i,"a.map", sep = ""))
  physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == i], 
                        Position = map(sparrow.abel)[chromosome(sparrow.abel) == i])
  counter <- counter +1
  physmap <- left_join(sparrow.map,physmap, by = "SNP.Name")
  product[[counter]] <- physmap[,c("SNP.Name", "Position", "Order", "analysisID")]
  
}

#Make Doubles
setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Data/Crimap Runs/Final/crimap")
sparrow.doubles <- list()
vec <- c(5:15, 17:28)
counter <- 0
for (i in (c(5:15, 17:28))) {
  sparrow.xovers <- parse_crossovers(chrompicfile = paste("chr",i,"a.cmp", sep = ""), familyPedigree = sparrow.famped)
  counter <- counter + 1
  sparrow.doubles[[counter]] <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = as.data.frame(product[which(vec == i)]))
  
}
#Make xovers
xover.list <- list()
counter <- 0
for (i in (c(5:15, 17:28))) {
  sparrow.xovers <- parse_crossovers(chrompicfile = paste("chr",i,"a.cmp", sep = ""), familyPedigree = sparrow.famped)
  counter <- counter + 1
  sparrow.xovers <- sparrow.xovers[,c(1:2,4:5,11:13)]
  xover.list[[counter]] <- sparrow.xovers
  
}
##
##
# Generating xovers.clean from doubles and parsed xover objects.
##
##
x <- do.call(rbind, sparrow.doubles)
sp.remove <- subset(x, Singleton == "yes")
xovers.clean.smlchr <- revise_double_crossovers(parsed.xovers = xovers, removesections = sp.remove)
xovers.clean.smlchr
#Do same for Large Chromosomes.
## Make Maps
##
##
## Make Maps
##
X <- list()
vec <- c(1:3)
counter <- 0
for (i in vec) {
  sparrow.mapa <- parse_map(mapfile = paste("chr", i,"a.map", sep = ""))
  sparrow.mapb <- parse_map(mapfile = paste("chr", i,"b.map", sep = ""))
  sparrow.mapc <- parse_map(mapfile = paste("chr", i,"c.map", sep = ""))
  physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == i], 
                        Position = map(sparrow.abel)[chromosome(sparrow.abel) == i])
  physmap$SNP.Name <- as.character(physmap$SNP.Name)
  counter <- counter + 1
  
  sparrow.map <- rbind(sparrow.mapa, sparrow.mapb,sparrow.mapc)
  sparrow.map <- sparrow.map[!duplicated(sparrow.map$SNP.Name),]
  physmap <- left_join(sparrow.map,physmap, by = "SNP.Name")
  test <- physmap[,c("SNP.Name", "Position", "Order", "analysisID")]
  test <- test[order(test$Position),]
  test$Order <- 0:(length(test$Position)-1)
  test$Order <- as.numeric(test$Order)
  test$analysisID <- paste(i,"a",sep = "")
  row.names(test) <- NULL
  X[[counter]] <- test
}
vec <- c(4,29)
counter <- 3
for (i in vec) {
  sparrow.mapa <- parse_map(mapfile = paste("chr", i,"a.map", sep = ""))
  sparrow.mapb <- parse_map(mapfile = paste("chr", i,"b.map", sep = ""))
  physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == i], 
                        Position = map(sparrow.abel)[chromosome(sparrow.abel) == i])
  physmap$SNP.Name <- as.character(physmap$SNP.Name)
  counter <- counter + 1
  sparrow.map <- rbind(sparrow.mapa, sparrow.mapb)
  sparrow.map <- sparrow.map[!duplicated(sparrow.map$SNP.Name),]
  physmap <- left_join(sparrow.map,physmap, by = "SNP.Name")            #Double check this function, might be misordering and/or doubling values; leftjoin alternative.
  test <- physmap[,c("SNP.Name", "Position", "Order", "analysisID")]
  test <- test[order(test$Position),]
  test$Order <- 0:(length(test$Position)-1)
  test$Order <- as.numeric(test$Order)
  test$analysisID <- paste(i,"a",sep = "")
  row.names(test) <- NULL
  X[[counter]] <- test
}
##
##Parse xovers and combine the seperate parts of the larger Chrs
##
##

vec<- c(1:3)
counter <- 23
for (i in vec) {
  a <- parse_crossovers(chrompicfile = paste("chr", i,"a.cmp", sep = ""), familyPedigree = sparrow.famped)
  b <- parse_crossovers(chrompicfile = paste("chr", i,"b.cmp", sep = ""), familyPedigree = sparrow.famped)
  c <- parse_crossovers(chrompicfile = paste("chr", i,"c.cmp", sep = ""), familyPedigree = sparrow.famped)
  #Remove overlaping portions in Data 
  a$data <- substr(a$data, 1, nchar(a$data)-101)
  b$data <- substr(b$data, 1, nchar(b$data)-101)
  z <- left_join(a, b, by=c("ANIMAL","parent","Family", "RRID"))
  z<- z[,c(1:5,11,14)]
  z <- left_join(z,c, by=c("ANIMAL","parent","Family", "RRID"))
  z<- z[,c(1:2,4:8,15:16)]
  z$UniqueID <- substr(z$UniqueID, 4, nchar(z$UniqueID))
  z$analysisID <- paste(i,"a",sep = "")
  z$data <- paste( paste(z$data.x, z$data.y, z$data, sep = ""))
  z <- z[,c(2:5,7:9)]
  counter <- counter + 1
  z <- z[,c(7,1:6)]
  xover.list[[counter]] <- z
}
vec<- c(4,29)
for (i in vec) {
  a <- parse_crossovers(chrompicfile = paste("chr", i,"a.cmp", sep = ""), familyPedigree = sparrow.famped)
  b <- parse_crossovers(chrompicfile = paste("chr", i,"b.cmp", sep = ""), familyPedigree = sparrow.famped)
  
  #Remove overlaping portions in Data 
  a$data <- substr(a$data, 1, nchar(a$data)-101)
  z <- left_join(a, b, by=c("ANIMAL","parent","Family", "RRID"))
  z<- z[,c(1:5,11,13:14)]
  z$UniqueID <- substr(z$UniqueID.x, 4, nchar(z$UniqueID.x))
  z$analysisID <- paste(i,"a",sep = "")
  z$data <- paste( paste(z$data.x, z$data.y, sep = ""))
  z <- z[,c(2,4:6,9:11)]
  counter <- counter + 1
  z <- z[,c(7,1:6)]
  xover.list[[counter]] <- z
}
xover.parsed <- do.call(rbind, xover.list)
##
##
##Perform check.double_crossovers() function and add to list with smaller chrs
##
##
totalvec <- c(5:15,17:28,1:4,29)
vec <- c(1:4,29)
counter <- 23
for (i in vec) {
  counter <- counter + 1
  sparrow.doubles[[counter]] <- check_double_crossovers(parsed.xovers = as.data.frame(xover.list[which(totalvec == i)]), physical.map = as.data.frame(X[which(vec == i)]))
}
sparrow.xovers.clean <- list()
counter<-0
vec <- c(5:15,17:28,1:4,29)
for(i in vec){
  counter <- counter + 1
  sparrow.remove <- subset(as.data.frame(sparrow.doubles[which(vec == i)]), Singleton == "yes")
  sparrow.xovers.clean[[counter]] <- revise_double_crossovers(parsed.xovers = as.data.frame(xover.list[which(vec == i)]), removesections = sparrow.remove)
}
all.xovers <- do.call(rbind,sparrow.xovers.clean)
doubles <- do.call(rbind, sparrow.doubles)
rm(sp.remove)
rm(sparrow.remove)
rm(a)
rm(b)
rm(c)
rm(physmap)
rm(sparrow.map)
rm(sparrow.mapa)
rm(sparrow.mapb)
rm(sparrow.mapc)
rm(X)
rm(test)
rm(sparrow.xovers)
rm(xover.list)
rm(xover.parsed)
rm(xovers.clean.smlchr)
rm(z)
rm(product)







##
##
##
#Take xovers from above and investigate crossovers in chosen bins.
##
##
#Split all xover for all chrs into 100snp chunks. 
#This method will keep the original whole snp data column with associated with each of the chunks
# it also assigns a chunk id


all.xovers.chunk <- list()
counter <- 1
row.names(all.xovers) <- NULL
for (i in 1:length(all.xovers$ANIMAL)) {
  #split the data column into the i'th row
x <- all.xovers$data[i]
snpchunks<-  sapply(seq(from=1, to=nchar(x), by=100), function(i) substr(x, i, i+99))
snpchunks <- as.data.frame(snpchunks)
snpchunks$snpchunks <- as.character(snpchunks$snpchunks)
row.names(snpchunks) <- NULL
snpchunks$chunkID <- 1:length(snpchunks$snpchunks)
snpchunks$chunkIDnumeric <- 1:length(snpchunks$snpchunks)
#take the variables for given i
pieces <- all.xovers[i,]
snpchunks$chunkID <- paste(pieces$analysisID,snpchunks$chunkID, sep = "_")
#combine the i'th variables with the chunked snps from data column
counter <- counter +1
all.xovers.chunk[[counter]] <- cbind(pieces, snpchunks)
}
all.xovers.chunk <- do.call(rbind, all.xovers.chunk)
#
##Check if chunks contain a crossover.##
#
# grepl() base r function searchs a character for specified values.
# WARNING: This method may cause me to miss crossovers if the grandparental origin switch occurs between chunks
# rather than within... This might be rare and therefore okay to use grepl() within method. However important to be aware of
# this limitation and consider writing a script which checks the adjecent chunks in some fashion for possible crossovers.
# If that proves difficult to write, another option could be to estimate the probability of this occuring and adjust
# overall counts accordingly. Issue might be rare enough to not justify doing this.

#Example; this line will check if a string "x" contains 0/o and 1/i, the way in which 'OR' is used allows a TRUE to pass if either a 'o'
# or a '0' is present and likewise for i/1. Then is both sides of the 'AND' operation are TRUE the output will be TRUE. Thus a vector or within
# crossover occuring TRUE/FALSE is created.
(grepl("o", x, fixed = FALSE) | grepl("0", x, fixed = FALSE))  &  (grepl("1", x, fixed = FALSE) | grepl("i", x, fixed = FALSE))
#Double check with Susie that o and i are interchangable with 0 and 1. 
#If not simply remove the 'OR' operations along with the arguments with o/i's.
#Loop it
crossover <- list()
counter <- 0
for (i in 1:length(all.xovers.chunk$snpchunks)) {
  x <- all.xovers.chunk$snpchunks[i]
  counter <- counter + 1
  crossover[[counter]] <- (grepl("o", x, fixed = FALSE) | grepl("0", x, fixed = FALSE))  &  (grepl("1", x, fixed = FALSE) | grepl("i", x, fixed = FALSE))
}
crossover <- do.call(rbind, crossover)
crossover <- as.vector(crossover)
all.xovers.chunk$chnkXover <- crossover
#Having these values as 1s and 0s may make data visualization easier than True/False
crossover <- as.integer(as.logical(crossover))
all.xovers.chunk$chnkXoverNum <- crossover


#Plot a chromosome
library(ggplot2)
test <- all.xovers.chunk[which(all.xovers.chunk$analysisID == "5a"),]

p <- ggplot(test, aes(x = chunkIDnumeric, y = chnkXoverNum)) +
  geom_bar(stat = "identity") +
  ggtitle(paste("Chromosome",gsub("a","",test$analysisID[1]), sep = " "))+
  xlab("Bin") + ylab("Number of Crossovers")
p












#Not a very good or efficient way to run on only small chromsomes but it does work.... 

#You can pick which chromosomes will be included in the loop by picking which are included in the smlchr vector.
#The other part: sapply(seq(from=1, to=nchar(x), by=10), function(i) substr(x, i, i+9))  ; this line determines the chunk lengths
#The length of the chunk is determined with the following two arguments: 
#by= "total length you want" ; i + "total length you chose MINUS 1"

smlchr <- c("6a","7a","9a","10a","11a","12a","13a","14a","15a","17a","18a","19a","20a","21a","22a","23a","24a","25a","26a","27a","28a")
sml.xovers.chunk <- list()
counter <- 1
row.names(all.xovers) <- NULL
for(j in 1:length(smlchr)){
  snglchrxover <- all.xovers[which(all.xovers$analysisID == smlchr[j]),]
for (i in 1:length(snglchrxover$ANIMAL)) {
  #split the data column into the i'th row
  x <- snglchrxover$data[i]
  snpchunks<-  sapply(seq(from=1, to=nchar(x), by=10), function(i) substr(x, i, i+9))
  snpchunks <- as.data.frame(snpchunks)
  snpchunks$snpchunks <- as.character(snpchunks$snpchunks)
  row.names(snpchunks) <- NULL
  snpchunks$chunkID <- 1:length(snpchunks$snpchunks)
  snpchunks$chunkIDnumeric <- 1:length(snpchunks$snpchunks)
  
  #take the variables for given i
  pieces <- snglchrxover[i,]
  snpchunks$chunkID <- paste(pieces$analysisID,snpchunks$chunkID, sep = "_")
  #combine the i'th variables with the chunked snps from data column
  counter <- counter +1
  sml.xovers.chunk[[counter]] <- cbind(pieces, snpchunks)
}
}
sml.xovers.chunk <- do.call(rbind, sml.xovers.chunk)
#Loop it
crossover <- list()
counter <- 0
for (i in 1:length(sml.xovers.chunk$snpchunks)) {
  x <- sml.xovers.chunk$snpchunks[i]
  counter <- counter + 1
  crossover[[counter]] <- (grepl("o", x, fixed = FALSE) | grepl("0", x, fixed = FALSE))  &  (grepl("1", x, fixed = FALSE) | grepl("i", x, fixed = FALSE))
}
crossover <- do.call(rbind, crossover)
crossover <- as.vector(crossover)
sml.xovers.chunk$chnkXover <- crossover
#Having these values as 1s and 0s may make data visualization easier than True/False
crossover <- as.integer(as.logical(crossover))
sml.xovers.chunk$chnkXoverNum <- crossover


test <- sml.xovers.chunk[which(sml.xovers.chunk$analysisID == "9a"),]
p <- ggplot(test, aes(x = chunkIDnumeric, y = chnkXoverNum)) +
  geom_bar(stat = "identity") +
  ggtitle(paste("Chromosome",gsub("a","",test$analysisID[1]), sep = " "))+
  xlab("Bin") + ylab("Number of Crossovers")
p
