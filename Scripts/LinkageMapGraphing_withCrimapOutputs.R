#Graphing Linkage Maps;uses genabel obj
library(GenABEL)
library(GenABEL.data)
library(ggplot2)
library(dplyr)
library(crimaptools)

#Load in the genABEL object used to run crimap, in following script the genABEL object is named "sparrow.abel"
#Set working directory to the one containing the crimap outputs
#The script parses the maps and combines the longer chromosomes which were run in chunks. 
setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Scripts/Model")
load("sparrowABEL.RData")
setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Data/Crimap Runs/Final/crimap")


#Parse Maps for smaller chrs
product <- list()
linkagemaps <- list()
counter <- 0
counter2 <- 0
for (i in (c(5:15, 17:28))) {
  
  sparrow.map <- parse_map(mapfile = paste("chr", i,"a.map", sep = ""))
  physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == i], 
                        Position = map(sparrow.abel)[chromosome(sparrow.abel) == i])
  counter <- counter +1
  physmap <- left_join(sparrow.map,physmap, by = "SNP.Name")
  product[[counter]] <- physmap[,c("SNP.Name", "Position", "Order", "analysisID")]
  
  counter2 <- counter2 + 1
  linkagemaps[[counter2]] <- physmap
}




#Parse maps for larger chrs

X <- list()
vec <- c(1:3)
counter <- 0
for (i in vec) {
  sparrow.mapa <- parse_map(mapfile = paste("chr", i,"a.map", sep = ""))
  sparrow.mapb <- parse_map(mapfile = paste("chr", i,"b.map", sep = ""))
  sparrow.mapc <- parse_map(mapfile = paste("chr", i,"c.map", sep = ""))
  
  
  sparrow.mapb$cMPosition <- sparrow.mapb$cMPosition + (sparrow.mapa[which(sparrow.mapb$SNP.Name[100] == sparrow.mapa$SNP.Name),]$cMPosition)
  sparrow.mapc$cMPosition <- sparrow.mapc$cMPosition + (sparrow.mapb[which(sparrow.mapc$SNP.Name[100] == sparrow.mapb$SNP.Name),]$cMPosition)
  
  #adding M&F cM pos within Chunks
  sparrow.mapb$cMPosition.Female <- sparrow.mapb$cMPosition.Female + (sparrow.mapa[which(sparrow.mapb$SNP.Name[100] == sparrow.mapa$SNP.Name),]$cMPosition.Female)
  sparrow.mapc$cMPosition.Female <- sparrow.mapc$cMPosition.Female + (sparrow.mapb[which(sparrow.mapc$SNP.Name[100] == sparrow.mapb$SNP.Name),]$cMPosition.Female)
  sparrow.mapb$cMPosition.Male <- sparrow.mapb$cMPosition.Male + (sparrow.mapa[which(sparrow.mapb$SNP.Name[100] == sparrow.mapa$SNP.Name),]$cMPosition.Male)
  sparrow.mapc$cMPosition.Male <- sparrow.mapc$cMPosition.Male + (sparrow.mapb[which(sparrow.mapc$SNP.Name[100] == sparrow.mapb$SNP.Name),]$cMPosition.Male)
  #Adding order within chunks
  sparrow.mapb$Order <- sparrow.mapb$Order + (sparrow.mapa[which(sparrow.mapb$SNP.Name[1] == sparrow.mapa$SNP.Name),]$Order)
  sparrow.mapc$Order <- sparrow.mapc$Order + (sparrow.mapb[which(sparrow.mapc$SNP.Name[1] == sparrow.mapb$SNP.Name),]$Order)
  
  
  
  physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == i], 
                        Position = map(sparrow.abel)[chromosome(sparrow.abel) == i])
  physmap$SNP.Name <- as.character(physmap$SNP.Name)
  counter <- counter + 1
  
  sparrow.map <- rbind(sparrow.mapa, sparrow.mapb,sparrow.mapc)
  sparrow.map <- sparrow.map[!duplicated(sparrow.map$SNP.Name),]
  physmap <- left_join(sparrow.map,physmap, by = "SNP.Name")
  
  counter2 <- counter2 + 1
  physmap$analysisID <- paste(i,"a",sep = "")
  linkagemaps[[counter2]] <- physmap
  
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
  
  #adding cM pos within chunks
  sparrow.mapb$cMPosition <- sparrow.mapb$cMPosition + max(sparrow.mapa$cMPosition)
  #adding M&F cM pos within Chunks
  sparrow.mapb$cMPosition.Female <- sparrow.mapb$cMPosition.Female + (sparrow.mapa[which(sparrow.mapb$SNP.Name[100] == sparrow.mapa$SNP.Name),]$cMPosition.Female)
  sparrow.mapb$cMPosition.Male <- sparrow.mapb$cMPosition.Male + (sparrow.mapa[which(sparrow.mapb$SNP.Name[100] == sparrow.mapa$SNP.Name),]$cMPosition.Male)
  #Adding order within chunks
  sparrow.mapb$Order <- sparrow.mapb$Order + (sparrow.mapa[which(sparrow.mapb$SNP.Name[1] == sparrow.mapa$SNP.Name),]$Order)
  
  physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == i], 
                        Position = map(sparrow.abel)[chromosome(sparrow.abel) == i])
  physmap$SNP.Name <- as.character(physmap$SNP.Name)
  counter <- counter + 1
  sparrow.map <- rbind(sparrow.mapa, sparrow.mapb)
  sparrow.map <- sparrow.map[!duplicated(sparrow.map$SNP.Name),]
  physmap <- left_join(sparrow.map,physmap, by = "SNP.Name")            #Double check this function, might be misordering and/or doubling values; leftjoin alternative.
  
  counter2 <- counter2 + 1
  physmap$analysisID <- paste(i,"a",sep = "")
  linkagemaps[[counter2]] <- physmap
  
  test <- physmap[,c("SNP.Name", "Position", "Order", "analysisID")]
  test <- test[order(test$Position),]
  test$Order <- 0:(length(test$Position)-1)
  test$Order <- as.numeric(test$Order)
  test$analysisID <- paste(i,"a",sep = "")
  row.names(test) <- NULL
  X[[counter]] <- test
}


test <- linkagemaps[[28]]

linkagemaps <- linkagemaps[c(24,25,26,27,1:23,28)]

linkmaps <- do.call(rbind, linkagemaps)
linkmaps$analysisID <- factor(linkmaps$analysisID, 
                             levels = c("1a","2a","3a","4a","5a","6a","7a",
                              "8a","9a","10a","11a","12a","13a","14a",
                              "15a","17a","18a","19a","20a","21a",
                              "22a","23a","24a","25a","26a","27a","28a","29a"))

levels(linkmaps$analysisID) <- c("Chr 1","Chr 2","Chr 3","Chr 4","Chr 5","Chr 6","Chr 7",
              "Chr 8","Chr 9","Chr 10","Chr 11","Chr 12","Chr 13","Chr 14",
              "Chr 15","Chr 17","Chr 18","Chr 19","Chr 20","Chr 21",
              "Chr 22","Chr 23","Chr 24","Chr 25","Chr 26","Chr 27","Chr 28","Chr 29")

p <-ggplot(data = linkmaps, aes(x = Position/1000000, y = cMPosition)) +
  geom_point(size = 1) +
  labs(y= "Linkage Map Length (cM)", x ="Chr Length (Mb)")
p + facet_wrap(~analysisID, scales = "free")


p <-ggplot(data = linkmaps) +
  geom_point(aes(x = Position/1000000, y = cMPosition.Male), color = "blue", size = 1) +
  geom_point(aes(x = Position/1000000, y = cMPosition.Female), color = "red", size = 1) +
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")
p + facet_wrap(~analysisID, scales = "free")



p <-ggplot(subset(linkmaps, analysisID %in% "1a")) +
  geom_point(aes(x = Position/1000000, y = cMPosition),size = 1) +
  labs(y= "Linkage Map Length (cM)", x ="Chr Length (Mb)")
p

chr.labs <- c("Chr 1","Chr 2","Chr 3","Chr 4")
p <-ggplot(subset(linkmaps, analysisID %in% c("Chr 1", "Chr 2", "Chr 3","Chr 4"))) +
  geom_point(aes(x = Position/1000000, y = cMPosition.Male), color = "blue", size = 1) +
  geom_point(aes(x = Position/1000000, y = cMPosition.Female), color = "red", size = 1) +
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)") +
  theme(axis.text = element_text( size = 32), axis.title = element_text(size = 32,face = "bold"), 
        plot.title = element_text(size = 32)
p + facet_wrap(~analysisID, scales = "free") + theme(strip.text.x = element_text(size = 32))
