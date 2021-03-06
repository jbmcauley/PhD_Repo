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