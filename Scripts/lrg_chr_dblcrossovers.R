rm(list = ls())
library(dplyr)
library(crimaptools)
library(GenABEL)

setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/Dblxoversremoved/")
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/crimap/crimap/crimap")
sparrow.famped <- read.table("FamPed_20200414.txt", header = TRUE)
badids_HighRecomb <- read.table("badids_HiRecomb.txt")
badids_HighRecomb <- as.integer(badids_HighRecomb$V1)
sparrow.famped <- sparrow.famped[which(!sparrow.famped$MOTHER %in% badids_HighRecomb),]
sparrow.famped <- sparrow.famped[which(!sparrow.famped$FATHER %in% badids_HighRecomb),]
sparrow.famped <- sparrow.famped[which(!sparrow.famped$ANIMAL %in% badids_HighRecomb),]
sparrow.famped <- sparrow.famped[sparrow.famped$Family %in% names(which(table(sparrow.famped$Family)>=5)),]
sparrow.famped$Family <- as.character(sparrow.famped$Family)
sparrow.famped$Family <- as.factor(sparrow.famped$Family)

setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Ped_Map_files")

goodsnps <- read.table("GoodSNPlist.txt")
goodsnps <- as.character(goodsnps$V1)
master_snpremove <- read.table("master_snpRemove.txt")
master_snpremove <- as.character(master_snpremove$V1)
load("sparrowgen.RData")
sparrowgen <- sparrowgen[,snpnames(sparrowgen) %in% goodsnps]
sparrowgen <- sparrowgen[,!snpnames(sparrowgen) %in% master_snpremove]
sparrow.abel <- sparrowgen

setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/crimap/crimap/")

sparrow.xover1a <- parse_crossovers(chrompicfile = "crimap/chr1a.cmp", familyPedigree = sparrow.famped)
sparrow.xover1b <- parse_crossovers(chrompicfile = "crimap/chr1b.cmp", familyPedigree = sparrow.famped)
sparrow.xover1c <- parse_crossovers(chrompicfile = "crimap/chr1c.cmp", familyPedigree = sparrow.famped)
sparrow.xover1d <- parse_crossovers(chrompicfile = "crimap/chr1d.cmp", familyPedigree = sparrow.famped)
sparrow.xover1e <- parse_crossovers(chrompicfile = "crimap/chr1e.cmp", familyPedigree = sparrow.famped)
sparrow.xover2a <- parse_crossovers(chrompicfile = "crimap/chr2a.cmp", familyPedigree = sparrow.famped)
sparrow.xover2b <- parse_crossovers(chrompicfile = "crimap/chr2b.cmp", familyPedigree = sparrow.famped)
sparrow.xover2c <- parse_crossovers(chrompicfile = "crimap/chr2c.cmp", familyPedigree = sparrow.famped)
sparrow.xover2d <- parse_crossovers(chrompicfile = "crimap/chr2d.cmp", familyPedigree = sparrow.famped)
sparrow.xover2e <- parse_crossovers(chrompicfile = "crimap/chr2e.cmp", familyPedigree = sparrow.famped)
sparrow.xover2f <- parse_crossovers(chrompicfile = "crimap/chr2f.cmp", familyPedigree = sparrow.famped)
#sparrow.xover2g <- parse_crossovers(chrompicfile = "crimap/chr2g.cmp", familyPedigree = sparrow.famped)
sparrow.xover3a <- parse_crossovers(chrompicfile = "crimap/chr3a.cmp", familyPedigree = sparrow.famped)
sparrow.xover3b <- parse_crossovers(chrompicfile = "crimap/chr3b.cmp", familyPedigree = sparrow.famped)
sparrow.xover3c <- parse_crossovers(chrompicfile = "crimap/chr3c.cmp", familyPedigree = sparrow.famped)
sparrow.xover3d <- parse_crossovers(chrompicfile = "crimap/chr3d.cmp", familyPedigree = sparrow.famped)
sparrow.xover3e <- parse_crossovers(chrompicfile = "crimap/chr3e.cmp", familyPedigree = sparrow.famped)
sparrow.xover4a <- parse_crossovers(chrompicfile = "crimap/chr4a.cmp", familyPedigree = sparrow.famped)
sparrow.xover4b <- parse_crossovers(chrompicfile = "crimap/chr4b.cmp", familyPedigree = sparrow.famped)
sparrow.xover4c <- parse_crossovers(chrompicfile = "crimap/chr4c.cmp", familyPedigree = sparrow.famped)
sparrow.xover5a <- parse_crossovers(chrompicfile = "crimap/chr5a.cmp", familyPedigree = sparrow.famped)
sparrow.xover5b <- parse_crossovers(chrompicfile = "crimap/chr5b.cmp", familyPedigree = sparrow.famped)
sparrow.xover5c <- parse_crossovers(chrompicfile = "crimap/chr5c.cmp", familyPedigree = sparrow.famped)
sparrow.xover6a <- parse_crossovers(chrompicfile = "crimap/chr6a.cmp", familyPedigree = sparrow.famped)
sparrow.xover6b <- parse_crossovers(chrompicfile = "crimap/chr6b.cmp", familyPedigree = sparrow.famped)
sparrow.xover7a <- parse_crossovers(chrompicfile = "crimap/chr7a.cmp", familyPedigree = sparrow.famped)
sparrow.xover7b <- parse_crossovers(chrompicfile = "crimap/chr7b.cmp", familyPedigree = sparrow.famped)
sparrow.xover8a <- parse_crossovers(chrompicfile = "crimap/chr8a.cmp", familyPedigree = sparrow.famped)
sparrow.xover8b <- parse_crossovers(chrompicfile = "crimap/chr8b.cmp", familyPedigree = sparrow.famped)
sparrow.xover8c <- parse_crossovers(chrompicfile = "crimap/chr8c.cmp", familyPedigree = sparrow.famped)
sparrow.xover29a <- parse_crossovers(chrompicfile = "crimap/chr29a.cmp", familyPedigree = sparrow.famped)
sparrow.xover29b <- parse_crossovers(chrompicfile = "crimap/chr29b.cmp", familyPedigree = sparrow.famped)
sparrow.xover29c <- parse_crossovers(chrompicfile = "crimap/chr29c.cmp", familyPedigree = sparrow.famped)



#Family, RRID, Parent, Unique ID, data, analysis ID)
#Only Data needs to be shaved/edited. Only 

library(dplyr)

sparrow.xover1a$data <- substr(sparrow.xover1a$data, 1, nchar(sparrow.xover1a$data)-101)
sparrow.xover1b$data <- substr(sparrow.xover1b$data, 1, nchar(sparrow.xover1b$data)-101)
sparrow.xover1c$data <- substr(sparrow.xover1c$data, 1, nchar(sparrow.xover1c$data)-101)
#sparrow.xover1d$data <- substr(sparrow.xover1d$data, 1, nchar(sparrow.xover1d$data)-101)


z <- merge(sparrow.xover1a, sparrow.xover1b, by=c("ANIMAL","parent","Family", "RRID"))

z$RecombCount.x <- NULL
z$No.Inf.Loci.x<- NULL
z$FATHER.x<- NULL
z$MOTHER.x<- NULL

z$First.Inf.Order.x<- NULL
z$Last.Inf.Order.x<- NULL

z$RecombCount.y<- NULL
z$No.Inf.Loci.y<- NULL
z$FATHER.y<- NULL
z$MOTHER.y<- NULL

z$First.Inf.Order.y<- NULL
z$Last.Inf.Order.y<- NULL

z <- merge(z,sparrow.xover1c, by=c("ANIMAL","parent","Family", "RRID"))
z$RecombCount<- NULL
z$No.Inf.Loci<- NULL
z$FATHER<- NULL
z$MOTHER<- NULL

z$First.Inf.Order<- NULL
z$Last.Inf.Order<- NULL

names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c")


z <- merge(z,sparrow.xover1d, by=c("ANIMAL","parent","Family", "RRID"))
z$RecombCount<- NULL
z$No.Inf.Loci<- NULL
z$FATHER<- NULL
z$MOTHER<- NULL

z$First.Inf.Order<- NULL
z$Last.Inf.Order<- NULL
names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c","data.d","analysisID.d","UniqueID.d")

#z <- merge(z,sparrow.xover1e, by=c("ANIMAL","parent","Family", "RRID"))
#z$RecombCount<- NULL
#z$No.Inf.Loci<- NULL
#z$FATHER<- NULL
#z$MOTHER<- NULL

#z$First.Inf.Order<- NULL
#z$Last.Inf.Order<- NULL
#names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c","data.d","analysisID.d","UniqueID.d","data.e","analysisID.e","UniqueID.e")

z$data <- paste(z$data.a, z$data.b, z$data.c, z$data.d, sep = "")
#z$data <- paste(z$data.a, z$data.b, z$data.c, z$data.d, z$data.e, sep = "")
nchar(z$data[1])
z$data.a <- NULL
z$data.b<- NULL
z$data.c<- NULL
z$data.d<- NULL
z$data.e<- NULL
z$analysisID.a <-NULL
z$analysisID.b <-NULL
z$analysisID.c <-NULL
z$analysisID.d <-NULL
z$analysisID.e <-NULL
z$analysisID <- "Chr1" 
z$UniqueID<- substr(z$UniqueID.a, 4, nchar(z$UniqueID.a))
z$UniqueID.a <-NULL
z$UniqueID.b <-NULL
z$UniqueID.c <-NULL
z$UniqueID.d <-NULL
z$UniqueID.e <-NULL

sparrow.doubles <- check_double_crossovers(parsed.xovers = z)

physmap <- data.frame(SNP.Name = snpnames(sparrowgen)[chromosome(sparrowgen) == 1], 
                      Position = map(sparrowgen)[chromosome(sparrowgen) == 1], 
                      Order = 1:length(which(chromosome(sparrowgen) == 1)), 
                      analysisID = "Chr1")
sparrow.doubles <- check_double_crossovers(parsed.xovers = z, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean1 <- revise_double_crossovers(parsed.xovers = z, removesections = sparrow.remove)
library(dplyr)
library(ggplot2)



write.table(sparrow.doubles,file = "sparrow.doubles.Chr1.txt", row.names = FALSE, col.names = TRUE)
write.table(sparrow.xovers.clean1, file = "sparrow.xovers.clean.Chr1.txt", row.names = FALSE, col.names = TRUE)
rm(sparrow.doubles)
rm(sparrow.xover1a)
rm(sparrow.xover1b)
rm(sparrow.xover1c)
rm(sparrow.xover1d)
rm(sparrow.xover1e)




length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)]) # 23366

sparrow.xover2a$data <- substr(sparrow.xover2a$data, 1, nchar(sparrow.xover2a$data)-101)
sparrow.xover2b$data <- substr(sparrow.xover2b$data, 1, nchar(sparrow.xover2b$data)-101)
sparrow.xover2c$data <- substr(sparrow.xover2c$data, 1, nchar(sparrow.xover2c$data)-101)
sparrow.xover2d$data <- substr(sparrow.xover2d$data, 1, nchar(sparrow.xover2d$data)-101)
sparrow.xover2e$data <- substr(sparrow.xover2e$data, 1, nchar(sparrow.xover2e$data)-101)
#sparrow.xover2f$data <- substr(sparrow.xover2f$data, 1, nchar(sparrow.xover2f$data)-101)


z <- merge(sparrow.xover2a, sparrow.xover2b, by=c("ANIMAL","parent","Family", "RRID"))

z$RecombCount.x <- NULL
z$No.Inf.Loci.x<- NULL
z$FATHER.x<- NULL
z$MOTHER.x<- NULL

z$First.Inf.Order.x<- NULL
z$Last.Inf.Order.x<- NULL

z$RecombCount.y<- NULL
z$No.Inf.Loci.y<- NULL
z$FATHER.y<- NULL
z$MOTHER.y<- NULL

z$First.Inf.Order.y<- NULL
z$Last.Inf.Order.y<- NULL

z <- merge(z,sparrow.xover2c, by=c("ANIMAL","parent","Family", "RRID"))
z$RecombCount<- NULL
z$No.Inf.Loci<- NULL
z$FATHER<- NULL
z$MOTHER<- NULL

z$First.Inf.Order<- NULL
z$Last.Inf.Order<- NULL

names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c")


z <- merge(z,sparrow.xover2d, by=c("ANIMAL","parent","Family", "RRID"))
z$RecombCount<- NULL
z$No.Inf.Loci<- NULL
z$FATHER<- NULL
z$MOTHER<- NULL

z$First.Inf.Order<- NULL
z$Last.Inf.Order<- NULL
names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c","data.d","analysisID.d","UniqueID.d")

z <- merge(z,sparrow.xover2e, by=c("ANIMAL","parent","Family", "RRID"))
z$RecombCount<- NULL
z$No.Inf.Loci<- NULL
z$FATHER<- NULL
z$MOTHER<- NULL

z$First.Inf.Order<- NULL
z$Last.Inf.Order<- NULL
names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c","data.d","analysisID.d","UniqueID.d","data.e","analysisID.e","UniqueID.e")


z <- merge(z,sparrow.xover2f, by=c("ANIMAL","parent","Family", "RRID"))
z$RecombCount<- NULL
z$No.Inf.Loci<- NULL
z$FATHER<- NULL
z$MOTHER<- NULL

z$First.Inf.Order<- NULL
z$Last.Inf.Order<- NULL
names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c","data.d","analysisID.d","UniqueID.d","data.e","analysisID.e","UniqueID.e","data.f","analysisID.f","UniqueID.f")

#z <- merge(z,sparrow.xover2g, by=c("ANIMAL","parent","Family", "RRID"))
#z$RecombCount<- NULL
#z$No.Inf.Loci<- NULL
#z$FATHER<- NULL
#z$MOTHER<- NULL

#z$First.Inf.Order<- NULL
#z$Last.Inf.Order<- NULL
#names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c","data.d","analysisID.d","UniqueID.d","data.e","analysisID.e","UniqueID.e","data.f","analysisID.f","UniqueID.f","data.g","analysisID.g","UniqueID.g")






z$data <- paste(z$data.a, z$data.b, z$data.c, z$data.d, z$data.e, z$data.f, sep = "")
#z$data <- paste(z$data.a, z$data.b, z$data.c, z$data.d, z$data.e, z$data.f, z$data.g, sep = "")
nchar(z$data[1])

z$data.a <- NULL
z$data.b<- NULL
z$data.c<- NULL
z$data.d<- NULL
z$data.e<- NULL
z$data.f<- NULL
#z$data.g<- NULL
z$analysisID.a <-NULL
z$analysisID.b <-NULL
z$analysisID.c <-NULL
z$analysisID.d <-NULL
z$analysisID.e <-NULL
z$analysisID.f <-NULL
#z$analysisID.g <-NULL

z$analysisID <- "Chr2"

z$UniqueID<- substr(z$UniqueID.a, 4, nchar(z$UniqueID.a))
z$UniqueID.a <-NULL
z$UniqueID.b <-NULL
z$UniqueID.c <-NULL
z$UniqueID.d <-NULL
z$UniqueID.e <-NULL
z$UniqueID.f <-NULL
#z$UniqueID.g <-NULL






sparrow.doubles <- check_double_crossovers(parsed.xovers = z)

physmap <- data.frame(SNP.Name = snpnames(sparrowgen)[chromosome(sparrowgen) == 2], 
                      Position = map(sparrowgen)[chromosome(sparrowgen) == 2], 
                      Order = 1:length(which(chromosome(sparrowgen) == 2)), 
                      analysisID = "Chr2")
sparrow.doubles <- check_double_crossovers(parsed.xovers = z, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean2 <- revise_double_crossovers(parsed.xovers = z, removesections = sparrow.remove)



write.table(sparrow.doubles,file = "sparrow.doubles.Chr2.txt", row.names = FALSE, col.names = TRUE)
write.table(sparrow.xovers.clean2, file = "sparrow.xovers.clean.Chr2.txt", row.names = FALSE, col.names = TRUE)
rm(sparrow.doubles)
rm(sparrow.xover2a)
rm(sparrow.xover2b)
rm(sparrow.xover2c)
rm(sparrow.xover2d)
rm(sparrow.xover2e)
rm(sparrow.xover2f)
rm(sparrow.xover2g)



length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)]) 


sparrow.xover3a$data <- substr(sparrow.xover3a$data, 1, nchar(sparrow.xover3a$data)-101)
sparrow.xover3b$data <- substr(sparrow.xover3b$data, 1, nchar(sparrow.xover3b$data)-101)
sparrow.xover3c$data <- substr(sparrow.xover3c$data, 1, nchar(sparrow.xover3c$data)-101)
#sparrow.xover3d$data <- substr(sparrow.xover3d$data, 1, nchar(sparrow.xover3d$data)-101)

z <- merge(sparrow.xover3a, sparrow.xover3b, by=c("ANIMAL","parent","Family", "RRID"))

z$RecombCount.x <- NULL
z$No.Inf.Loci.x<- NULL
z$FATHER.x<- NULL
z$MOTHER.x<- NULL

z$First.Inf.Order.x<- NULL
z$Last.Inf.Order.x<- NULL

z$RecombCount.y<- NULL
z$No.Inf.Loci.y<- NULL
z$FATHER.y<- NULL
z$MOTHER.y<- NULL

z$First.Inf.Order.y<- NULL
z$Last.Inf.Order.y<- NULL

z <- merge(z,sparrow.xover3c, by=c("ANIMAL","parent","Family", "RRID"))
z$RecombCount<- NULL
z$No.Inf.Loci<- NULL
z$FATHER<- NULL
z$MOTHER<- NULL

z$First.Inf.Order<- NULL
z$Last.Inf.Order<- NULL

names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c")


z <- merge(z,sparrow.xover3d, by=c("ANIMAL","parent","Family", "RRID"))
z$RecombCount<- NULL
z$No.Inf.Loci<- NULL
z$FATHER<- NULL
z$MOTHER<- NULL

z$First.Inf.Order<- NULL
z$Last.Inf.Order<- NULL
names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c","data.d","analysisID.d","UniqueID.d")

#z <- merge(z,sparrow.xover3e, by=c("ANIMAL","parent","Family", "RRID"))
#z$RecombCount<- NULL
#z$No.Inf.Loci<- NULL
#z$FATHER<- NULL
#z$MOTHER<- NULL

#z$First.Inf.Order<- NULL
#z$Last.Inf.Order<- NULL
#names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c","data.d","analysisID.d","UniqueID.d","data.e","analysisID.e","UniqueID.e")





z$data <- paste(z$data.a, z$data.b, z$data.c, z$data.d, sep = "")
#z$data <- paste(z$data.a, z$data.b, z$data.c, z$data.d, z$data.e, sep = "")
nchar(z$data[1])
z$data.a <- NULL
z$data.b<- NULL
z$data.c<- NULL
z$data.d<- NULL
z$data.e<- NULL
z$analysisID.a <-NULL
z$analysisID.b <-NULL
z$analysisID.c <-NULL
z$analysisID.d <-NULL
z$analysisID.e <-NULL
z$analysisID <- "Chr3"
z$UniqueID<- substr(z$UniqueID.a, 4, nchar(z$UniqueID.a))
z$UniqueID.a <-NULL
z$UniqueID.b <-NULL
z$UniqueID.c <-NULL
z$UniqueID.d <-NULL
z$UniqueID.e <-NULL


sparrow.doubles <- check_double_crossovers(parsed.xovers = z)

physmap <- data.frame(SNP.Name = snpnames(sparrowgen)[chromosome(sparrowgen) == 3], 
                      Position = map(sparrowgen)[chromosome(sparrowgen) == 3], 
                      Order = 1:length(which(chromosome(sparrowgen) == 3)), 
                      analysisID = "Chr3")
sparrow.doubles <- check_double_crossovers(parsed.xovers = z, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean3 <- revise_double_crossovers(parsed.xovers = z, removesections = sparrow.remove)


write.table(sparrow.doubles,file = "sparrow.doubles.Chr3.txt", row.names = FALSE, col.names = TRUE)
write.table(sparrow.xovers.clean3, file = "sparrow.xovers.clean.Chr3.txt", row.names = FALSE, col.names = TRUE)
rm(sparrow.doubles)
rm(sparrow.xover3a)
rm(sparrow.xover3b)
rm(sparrow.xover3c)
rm(sparrow.xover3d)
rm(sparrow.xover3e)
rm(sparrow.xovers.clean)
rm(sparrow.remove)








sparrow.xover4a$data <- substr(sparrow.xover4a$data, 1, nchar(sparrow.xover4a$data)-101)
sparrow.xover4b$data <- substr(sparrow.xover4b$data, 1, nchar(sparrow.xover4b$data)-101)


z <- merge(sparrow.xover4a, sparrow.xover4b, by=c("ANIMAL","parent","Family", "RRID"))

z$RecombCount.x <- NULL
z$No.Inf.Loci.x<- NULL
z$FATHER.x<- NULL
z$MOTHER.x<- NULL

z$First.Inf.Order.x<- NULL
z$Last.Inf.Order.x<- NULL

z$RecombCount.y<- NULL
z$No.Inf.Loci.y<- NULL
z$FATHER.y<- NULL
z$MOTHER.y<- NULL

z$First.Inf.Order.y<- NULL
z$Last.Inf.Order.y<- NULL

z <- merge(z,sparrow.xover4c, by=c("ANIMAL","parent","Family", "RRID"))
z$RecombCount<- NULL
z$No.Inf.Loci<- NULL
z$FATHER<- NULL
z$MOTHER<- NULL

z$First.Inf.Order<- NULL
z$Last.Inf.Order<- NULL

names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c")


z$data <- paste(z$data.a, z$data.b, z$data.c, sep = "")
nchar(z$data[1])
z$data.a <- NULL
z$data.b<- NULL
z$data.c<- NULL
z$data.d<- NULL
z$data.e<- NULL
z$analysisID.a <-NULL
z$analysisID.b <-NULL
z$analysisID.c <-NULL
z$analysisID.d <-NULL
z$analysisID.e <-NULL
z$analysisID <- "Chr4"
z$UniqueID<- substr(z$UniqueID.a, 4, nchar(z$UniqueID.a))
z$UniqueID.a <-NULL
z$UniqueID.b <-NULL
z$UniqueID.c <-NULL
z$UniqueID.d <-NULL
z$UniqueID.e <-NULL


sparrow.doubles <- check_double_crossovers(parsed.xovers = z)

physmap <- data.frame(SNP.Name = snpnames(sparrowgen)[chromosome(sparrowgen) == 4], 
                      Position = map(sparrowgen)[chromosome(sparrowgen) == 4], 
                      Order = 1:length(which(chromosome(sparrowgen) == 4)), 
                      analysisID = "Chr4")
sparrow.doubles <- check_double_crossovers(parsed.xovers = z, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean4 <- revise_double_crossovers(parsed.xovers = z, removesections = sparrow.remove)


write.table(sparrow.doubles,file = "sparrow.doubles.Chr4.txt", row.names = FALSE, col.names = TRUE)
write.table(sparrow.xovers.clean4, file = "sparrow.xovers.clean.Chr4.txt", row.names = FALSE, col.names = TRUE)
rm(sparrow.doubles)
rm(sparrow.xover4a)
rm(sparrow.xover4b)
rm(sparrow.xover4c)
rm(sparrow.xovers.clean)
rm(sparrow.remove)








sparrow.xover5a$data <- substr(sparrow.xover5a$data, 1, nchar(sparrow.xover5a$data)-101)
sparrow.xover5b$data <- substr(sparrow.xover5b$data, 1, nchar(sparrow.xover5b$data)-101)

z <- merge(sparrow.xover5a, sparrow.xover5b, by=c("ANIMAL","parent","Family", "RRID"))

z$RecombCount.x <- NULL
z$No.Inf.Loci.x<- NULL
z$FATHER.x<- NULL
z$MOTHER.x<- NULL

z$First.Inf.Order.x<- NULL
z$Last.Inf.Order.x<- NULL

z$RecombCount.y<- NULL
z$No.Inf.Loci.y<- NULL
z$FATHER.y<- NULL
z$MOTHER.y<- NULL

z$First.Inf.Order.y<- NULL
z$Last.Inf.Order.y<- NULL

z <- merge(z,sparrow.xover5c, by=c("ANIMAL","parent","Family", "RRID"))
z$RecombCount<- NULL
z$No.Inf.Loci<- NULL
z$FATHER<- NULL
z$MOTHER<- NULL

z$First.Inf.Order<- NULL
z$Last.Inf.Order<- NULL

names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c")


z$data <- paste(z$data.a, z$data.b, z$data.c, sep = "")
nchar(z$data[1])
z$data.a <- NULL
z$data.b<- NULL
z$data.c<- NULL
z$data.d<- NULL
z$data.e<- NULL
z$analysisID.a <-NULL
z$analysisID.b <-NULL
z$analysisID.c <-NULL
z$analysisID.d <-NULL
z$analysisID.e <-NULL
z$analysisID <- "Chr5"
z$UniqueID<- substr(z$UniqueID.a, 4, nchar(z$UniqueID.a))
z$UniqueID.a <-NULL
z$UniqueID.b <-NULL
z$UniqueID.c <-NULL
z$UniqueID.d <-NULL
z$UniqueID.e <-NULL


sparrow.doubles <- check_double_crossovers(parsed.xovers = z)

physmap <- data.frame(SNP.Name = snpnames(sparrowgen)[chromosome(sparrowgen) == 5], 
                      Position = map(sparrowgen)[chromosome(sparrowgen) == 5], 
                      Order = 1:length(which(chromosome(sparrowgen) == 5)), 
                      analysisID = "Chr5")
sparrow.doubles <- check_double_crossovers(parsed.xovers = z, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean5 <- revise_double_crossovers(parsed.xovers = z, removesections = sparrow.remove)


write.table(sparrow.doubles,file = "sparrow.doubles.Chr5.txt", row.names = FALSE, col.names = TRUE)
write.table(sparrow.xovers.clean5, file = "sparrow.xovers.clean.Chr5.txt", row.names = FALSE, col.names = TRUE)
#sparrow.doubles.chr5 <- read.table("sparrow.doubles.chr5.txt", header = TRUE)
#sparrow.xovers.clean5 <- read.table("sparrow.xovers.clean.Chr5.txt", header = TRUE)

rm(sparrow.doubles)
rm(sparrow.xover5a)
rm(sparrow.xover5b)
rm(sparrow.xover5c)
rm(sparrow.xovers.clean)
rm(sparrow.remove)






sparrow.xover6a$data <- substr(sparrow.xover6a$data, 1, nchar(sparrow.xover6a$data)-101)


z <- merge(sparrow.xover6a, sparrow.xover6b, by=c("ANIMAL","parent","Family", "RRID"))

z$RecombCount.x <- NULL
z$No.Inf.Loci.x<- NULL
z$FATHER.x<- NULL
z$MOTHER.x<- NULL

z$First.Inf.Order.x<- NULL
z$Last.Inf.Order.x<- NULL

z$RecombCount.y<- NULL
z$No.Inf.Loci.y<- NULL
z$FATHER.y<- NULL
z$MOTHER.y<- NULL

z$First.Inf.Order.y<- NULL
z$Last.Inf.Order.y<- NULL

names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b")


z$data <- paste(z$data.a, z$data.b, sep = "")
nchar(z$data[1])
z$data.a <- NULL
z$data.b<- NULL
z$data.c<- NULL
z$data.d<- NULL
z$data.e<- NULL
z$analysisID.a <-NULL
z$analysisID.b <-NULL
z$analysisID <- "Chr6"
z$UniqueID<- substr(z$UniqueID.a, 4, nchar(z$UniqueID.a))
z$UniqueID.a <-NULL
z$UniqueID.b <-NULL
z$UniqueID.c <-NULL
z$UniqueID.d <-NULL
z$UniqueID.e <-NULL


sparrow.doubles <- check_double_crossovers(parsed.xovers = z)

physmap <- data.frame(SNP.Name = snpnames(sparrowgen)[chromosome(sparrowgen) == 6], 
                      Position = map(sparrowgen)[chromosome(sparrowgen) == 6], 
                      Order = 1:length(which(chromosome(sparrowgen) == 6)), 
                      analysisID = "Chr6")
sparrow.doubles <- check_double_crossovers(parsed.xovers = z, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean6 <- revise_double_crossovers(parsed.xovers = z, removesections = sparrow.remove)


write.table(sparrow.doubles,file = "sparrow.doubles.Chr6.txt", row.names = FALSE, col.names = TRUE)
write.table(sparrow.xovers.clean6, file = "sparrow.xovers.clean.Chr6.txt", row.names = FALSE, col.names = TRUE)
rm(sparrow.doubles)
rm(sparrow.xover6a)
rm(sparrow.xover6b)
rm(sparrow.xovers.clean)
rm(sparrow.remove)









sparrow.xover7a$data <- substr(sparrow.xover7a$data, 1, nchar(sparrow.xover7a$data)-101)


z <- merge(sparrow.xover7a, sparrow.xover7b, by=c("ANIMAL","parent","Family", "RRID"))

z$RecombCount.x <- NULL
z$No.Inf.Loci.x<- NULL
z$FATHER.x<- NULL
z$MOTHER.x<- NULL

z$First.Inf.Order.x<- NULL
z$Last.Inf.Order.x<- NULL

z$RecombCount.y<- NULL
z$No.Inf.Loci.y<- NULL
z$FATHER.y<- NULL
z$MOTHER.y<- NULL

z$First.Inf.Order.y<- NULL
z$Last.Inf.Order.y<- NULL

names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b")


z$data <- paste(z$data.a, z$data.b, sep = "")
nchar(z$data[1])
z$data.a <- NULL
z$data.b<- NULL
z$data.c<- NULL
z$data.d<- NULL
z$data.e<- NULL
z$analysisID.a <-NULL
z$analysisID.b <-NULL
z$analysisID <- "Chr7"
z$UniqueID<- substr(z$UniqueID.a, 4, nchar(z$UniqueID.a))
z$UniqueID.a <-NULL
z$UniqueID.b <-NULL
z$UniqueID.c <-NULL
z$UniqueID.d <-NULL
z$UniqueID.e <-NULL


sparrow.doubles <- check_double_crossovers(parsed.xovers = z)

physmap <- data.frame(SNP.Name = snpnames(sparrowgen)[chromosome(sparrowgen) == 7], 
                      Position = map(sparrowgen)[chromosome(sparrowgen) == 7], 
                      Order = 1:length(which(chromosome(sparrowgen) == 7)), 
                      analysisID = "Chr7")
sparrow.doubles <- check_double_crossovers(parsed.xovers = z, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean7 <- revise_double_crossovers(parsed.xovers = z, removesections = sparrow.remove)


write.table(sparrow.doubles,file = "sparrow.doubles.Chr7.txt", row.names = FALSE, col.names = TRUE)
write.table(sparrow.xovers.clean7, file = "sparrow.xovers.clean.Chr7.txt", row.names = FALSE, col.names = TRUE)
rm(sparrow.doubles)
rm(sparrow.xover7a)
rm(sparrow.xover7b)
rm(sparrow.xovers.clean)
rm(sparrow.remove)







sparrow.xover8a$data <- substr(sparrow.xover8a$data, 1, nchar(sparrow.xover8a$data)-101)
#sparrow.xover8b$data <- substr(sparrow.xover8b$data, 1, nchar(sparrow.xover8b$data)-101)

z <- merge(sparrow.xover8a, sparrow.xover8b, by=c("ANIMAL","parent","Family", "RRID"))

z$RecombCount.x <- NULL
z$No.Inf.Loci.x<- NULL
z$FATHER.x<- NULL
z$MOTHER.x<- NULL

z$First.Inf.Order.x<- NULL
z$Last.Inf.Order.x<- NULL

z$RecombCount.y<- NULL
z$No.Inf.Loci.y<- NULL
z$FATHER.y<- NULL
z$MOTHER.y<- NULL

z$First.Inf.Order.y<- NULL
z$Last.Inf.Order.y<- NULL
names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b")

#z <- merge(z,sparrow.xover8c, by=c("ANIMAL","parent","Family", "RRID"))
#z$RecombCount<- NULL
#z$No.Inf.Loci<- NULL
#z$FATHER<- NULL
#z$MOTHER<- NULL

#z$First.Inf.Order<- NULL
#z$Last.Inf.Order<- NULL

#names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c")

z$data <- paste(z$data.a, z$data.b, sep = "")
#z$data <- paste(z$data.a, z$data.b, z$data.c, sep = "")
nchar(z$data[1])
z$data.a <- NULL
z$data.b<- NULL
z$data.c<- NULL
z$data.d<- NULL
z$data.e<- NULL
z$analysisID.a <-NULL
z$analysisID.b <-NULL
z$analysisID.c <-NULL
z$analysisID.d <-NULL
z$analysisID.e <-NULL
z$analysisID <- "Chr8"
z$UniqueID<- substr(z$UniqueID.a, 4, nchar(z$UniqueID.a))
z$UniqueID.a <-NULL
z$UniqueID.b <-NULL
z$UniqueID.c <-NULL
z$UniqueID.d <-NULL
z$UniqueID.e <-NULL


sparrow.doubles <- check_double_crossovers(parsed.xovers = z)

physmap <- data.frame(SNP.Name = snpnames(sparrowgen)[chromosome(sparrowgen) == 8], 
                      Position = map(sparrowgen)[chromosome(sparrowgen) == 8], 
                      Order = 1:length(which(chromosome(sparrowgen) == 8)), 
                      analysisID = "Chr8")
sparrow.doubles <- check_double_crossovers(parsed.xovers = z, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean8 <- revise_double_crossovers(parsed.xovers = z, removesections = sparrow.remove)


write.table(sparrow.doubles,file = "sparrow.doubles.Chr8.txt", row.names = FALSE, col.names = TRUE)
write.table(sparrow.xovers.clean8, file = "sparrow.xovers.clean.Chr8.txt", row.names = FALSE, col.names = TRUE)
rm(sparrow.doubles)
rm(sparrow.xover8a)
rm(sparrow.xover8b)
rm(sparrow.xover8c)
rm(sparrow.xovers.clean)
rm(sparrow.remove)








sparrow.xover29a$data <- substr(sparrow.xover29a$data, 1, nchar(sparrow.xover29a$data)-101)
sparrow.xover29b$data <- substr(sparrow.xover29b$data, 1, nchar(sparrow.xover29b$data)-101)

z <- merge(sparrow.xover29a, sparrow.xover29b, by=c("ANIMAL","parent","Family", "RRID"))

z$RecombCount.x <- NULL
z$No.Inf.Loci.x<- NULL
z$FATHER.x<- NULL
z$MOTHER.x<- NULL

z$First.Inf.Order.x<- NULL
z$Last.Inf.Order.x<- NULL

z$RecombCount.y<- NULL
z$No.Inf.Loci.y<- NULL
z$FATHER.y<- NULL
z$MOTHER.y<- NULL

z$First.Inf.Order.y<- NULL
z$Last.Inf.Order.y<- NULL

z <- merge(z,sparrow.xover29c, by=c("ANIMAL","parent","Family", "RRID"))
z$RecombCount<- NULL
z$No.Inf.Loci<- NULL
z$FATHER<- NULL
z$MOTHER<- NULL

z$First.Inf.Order<- NULL
z$Last.Inf.Order<- NULL

names(z) <- c("ANIMAL", "parent","Family","RRID","data.a","analysisID.a","UniqueID.a","data.b","analysisID.b","UniqueID.b","data.c","analysisID.c","UniqueID.c")


z$data <- paste(z$data.a, z$data.b, z$data.c, sep = "")
nchar(z$data[1])
z$data.a <- NULL
z$data.b<- NULL
z$data.c<- NULL
z$data.d<- NULL
z$data.e<- NULL
z$analysisID.a <-NULL
z$analysisID.b <-NULL
z$analysisID.c <-NULL
z$analysisID.d <-NULL
z$analysisID.e <-NULL
z$analysisID <- "Chr29"
z$UniqueID<- substr(z$UniqueID.a, 4, nchar(z$UniqueID.a))
z$UniqueID.a <-NULL
z$UniqueID.b <-NULL
z$UniqueID.c <-NULL
z$UniqueID.d <-NULL
z$UniqueID.e <-NULL


sparrow.doubles <- check_double_crossovers(parsed.xovers = z)

physmap <- data.frame(SNP.Name = snpnames(sparrowgen)[chromosome(sparrowgen) == 29], 
                      Position = map(sparrowgen)[chromosome(sparrowgen) == 29], 
                      Order = 1:length(which(chromosome(sparrowgen) == 29)), 
                      analysisID = "Chr29")
sparrow.doubles <- check_double_crossovers(parsed.xovers = z, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean29 <- revise_double_crossovers(parsed.xovers = z, removesections = sparrow.remove)


write.table(sparrow.doubles,file = "sparrow.doubles.Chr29.txt", row.names = FALSE, col.names = TRUE)
write.table(sparrow.xovers.clean29, file = "sparrow.xovers.clean.Chr29.txt", row.names = FALSE, col.names = TRUE)
rm(sparrow.doubles)
rm(sparrow.xover29a)
rm(sparrow.xover29b)
rm(sparrow.xover29c)
rm(sparrow.xovers.clean)
rm(sparrow.remove)
rm(z)







sparrow.xovers.clean9$No.Inf.Loci <- NULL
sparrow.xovers.clean9$First.Inf.Order <- NULL
sparrow.xovers.clean9$Last.Inf.Order<- NULL
sparrow.xovers.clean9$FATHER<- NULL
sparrow.xovers.clean9$MOTHER<- NULL
sparrow.xovers.clean9$UniqueID <- substr(sparrow.xovers.clean9$UniqueID, 4, nchar(sparrow.xovers.clean9$UniqueID))
sparrow.xovers.clean9$analysisID <- "Chr9"

sparrow.xovers.clean10$No.Inf.Loci <- NULL
sparrow.xovers.clean10$First.Inf.Order <- NULL
sparrow.xovers.clean10$Last.Inf.Order<- NULL
sparrow.xovers.clean10$FATHER<- NULL
sparrow.xovers.clean10$MOTHER<- NULL
sparrow.xovers.clean10$UniqueID <- substr(sparrow.xovers.clean10$UniqueID, 5, nchar(sparrow.xovers.clean10$UniqueID))
sparrow.xovers.clean10$analysisID <- "Chr10"

sparrow.xovers.clean11$No.Inf.Loci <- NULL
sparrow.xovers.clean11$First.Inf.Order <- NULL
sparrow.xovers.clean11$Last.Inf.Order<- NULL
sparrow.xovers.clean11$FATHER<- NULL
sparrow.xovers.clean11$MOTHER<- NULL
sparrow.xovers.clean11$UniqueID <- substr(sparrow.xovers.clean11$UniqueID, 5, nchar(sparrow.xovers.clean11$UniqueID))
sparrow.xovers.clean11$analysisID <- "Chr11"

sparrow.xovers.clean12$No.Inf.Loci <- NULL
sparrow.xovers.clean12$First.Inf.Order <- NULL
sparrow.xovers.clean12$Last.Inf.Order<- NULL
sparrow.xovers.clean12$FATHER<- NULL
sparrow.xovers.clean12$MOTHER<- NULL
sparrow.xovers.clean12$UniqueID <- substr(sparrow.xovers.clean12$UniqueID, 5, nchar(sparrow.xovers.clean12$UniqueID))
sparrow.xovers.clean12$analysisID <- "Chr12"

sparrow.xovers.clean13$No.Inf.Loci <- NULL
sparrow.xovers.clean13$First.Inf.Order <- NULL
sparrow.xovers.clean13$Last.Inf.Order<- NULL
sparrow.xovers.clean13$FATHER<- NULL
sparrow.xovers.clean13$MOTHER<- NULL
sparrow.xovers.clean13$UniqueID <- substr(sparrow.xovers.clean13$UniqueID, 5, nchar(sparrow.xovers.clean13$UniqueID))
sparrow.xovers.clean13$analysisID <- "Chr13"

sparrow.xovers.clean14$No.Inf.Loci <- NULL
sparrow.xovers.clean14$First.Inf.Order <- NULL
sparrow.xovers.clean14$Last.Inf.Order<- NULL
sparrow.xovers.clean14$FATHER<- NULL
sparrow.xovers.clean14$MOTHER<- NULL
sparrow.xovers.clean14$UniqueID <- substr(sparrow.xovers.clean14$UniqueID, 5, nchar(sparrow.xovers.clean14$UniqueID))
sparrow.xovers.clean14$analysisID <- "Chr14"

sparrow.xovers.clean15$No.Inf.Loci <- NULL
sparrow.xovers.clean15$First.Inf.Order <- NULL
sparrow.xovers.clean15$Last.Inf.Order<- NULL
sparrow.xovers.clean15$FATHER<- NULL
sparrow.xovers.clean15$MOTHER<- NULL
sparrow.xovers.clean15$UniqueID <- substr(sparrow.xovers.clean15$UniqueID, 5, nchar(sparrow.xovers.clean15$UniqueID))
sparrow.xovers.clean15$analysisID <- "Chr15"

sparrow.xovers.clean17$No.Inf.Loci <- NULL
sparrow.xovers.clean17$First.Inf.Order <- NULL
sparrow.xovers.clean17$Last.Inf.Order<- NULL
sparrow.xovers.clean17$FATHER<- NULL
sparrow.xovers.clean17$MOTHER<- NULL
sparrow.xovers.clean17$UniqueID <- substr(sparrow.xovers.clean17$UniqueID, 5, nchar(sparrow.xovers.clean17$UniqueID))
sparrow.xovers.clean17$analysisID <- "Chr17"

sparrow.xovers.clean18$No.Inf.Loci <- NULL
sparrow.xovers.clean18$First.Inf.Order <- NULL
sparrow.xovers.clean18$Last.Inf.Order<- NULL
sparrow.xovers.clean18$FATHER<- NULL
sparrow.xovers.clean18$MOTHER<- NULL
sparrow.xovers.clean18$UniqueID <- substr(sparrow.xovers.clean18$UniqueID, 5, nchar(sparrow.xovers.clean18$UniqueID))
sparrow.xovers.clean18$analysisID <- "Chr18"

sparrow.xovers.clean19$No.Inf.Loci <- NULL
sparrow.xovers.clean19$First.Inf.Order <- NULL
sparrow.xovers.clean19$Last.Inf.Order<- NULL
sparrow.xovers.clean19$FATHER<- NULL
sparrow.xovers.clean19$MOTHER<- NULL
sparrow.xovers.clean19$UniqueID <- substr(sparrow.xovers.clean19$UniqueID, 5, nchar(sparrow.xovers.clean19$UniqueID))
sparrow.xovers.clean19$analysisID <- "Chr19"

sparrow.xovers.clean20$No.Inf.Loci <- NULL
sparrow.xovers.clean20$First.Inf.Order <- NULL
sparrow.xovers.clean20$Last.Inf.Order<- NULL
sparrow.xovers.clean20$FATHER<- NULL
sparrow.xovers.clean20$MOTHER<- NULL
sparrow.xovers.clean20$UniqueID <- substr(sparrow.xovers.clean20$UniqueID, 5, nchar(sparrow.xovers.clean20$UniqueID))
sparrow.xovers.clean20$analysisID <- "Chr20"

sparrow.xovers.clean21$No.Inf.Loci <- NULL
sparrow.xovers.clean21$First.Inf.Order <- NULL
sparrow.xovers.clean21$Last.Inf.Order<- NULL
sparrow.xovers.clean21$FATHER<- NULL
sparrow.xovers.clean21$MOTHER<- NULL
sparrow.xovers.clean21$UniqueID <- substr(sparrow.xovers.clean21$UniqueID, 5, nchar(sparrow.xovers.clean21$UniqueID))
sparrow.xovers.clean21$analysisID <- "Chr21"

sparrow.xovers.clean22$No.Inf.Loci <- NULL
sparrow.xovers.clean22$First.Inf.Order <- NULL
sparrow.xovers.clean22$Last.Inf.Order<- NULL
sparrow.xovers.clean22$FATHER<- NULL
sparrow.xovers.clean22$MOTHER<- NULL
sparrow.xovers.clean22$UniqueID <- substr(sparrow.xovers.clean22$UniqueID, 5, nchar(sparrow.xovers.clean22$UniqueID))
sparrow.xovers.clean22$analysisID <- "Chr22"

sparrow.xovers.clean23$No.Inf.Loci <- NULL
sparrow.xovers.clean23$First.Inf.Order <- NULL
sparrow.xovers.clean23$Last.Inf.Order<- NULL
sparrow.xovers.clean23$FATHER<- NULL
sparrow.xovers.clean23$MOTHER<- NULL
sparrow.xovers.clean23$UniqueID <- substr(sparrow.xovers.clean23$UniqueID, 5, nchar(sparrow.xovers.clean23$UniqueID))
sparrow.xovers.clean23$analysisID <- "Chr23"

sparrow.xovers.clean24$No.Inf.Loci <- NULL
sparrow.xovers.clean24$First.Inf.Order <- NULL
sparrow.xovers.clean24$Last.Inf.Order<- NULL
sparrow.xovers.clean24$FATHER<- NULL
sparrow.xovers.clean24$MOTHER<- NULL
sparrow.xovers.clean24$UniqueID <- substr(sparrow.xovers.clean24$UniqueID, 5, nchar(sparrow.xovers.clean24$UniqueID))
sparrow.xovers.clean24$analysisID <- "Chr24"

sparrow.xovers.clean25$No.Inf.Loci <- NULL
sparrow.xovers.clean25$First.Inf.Order <- NULL
sparrow.xovers.clean25$Last.Inf.Order<- NULL
sparrow.xovers.clean25$FATHER<- NULL
sparrow.xovers.clean25$MOTHER<- NULL
sparrow.xovers.clean25$UniqueID <- substr(sparrow.xovers.clean25$UniqueID, 5, nchar(sparrow.xovers.clean25$UniqueID))
sparrow.xovers.clean25$analysisID <- "Chr25"

sparrow.xovers.clean26$No.Inf.Loci <- NULL
sparrow.xovers.clean26$First.Inf.Order <- NULL
sparrow.xovers.clean26$Last.Inf.Order<- NULL
sparrow.xovers.clean26$FATHER<- NULL
sparrow.xovers.clean26$MOTHER<- NULL
sparrow.xovers.clean26$UniqueID <- substr(sparrow.xovers.clean26$UniqueID, 5, nchar(sparrow.xovers.clean26$UniqueID))
sparrow.xovers.clean26$analysisID <- "Chr26"

sparrow.xovers.clean27$No.Inf.Loci <- NULL
sparrow.xovers.clean27$First.Inf.Order <- NULL
sparrow.xovers.clean27$Last.Inf.Order<- NULL
sparrow.xovers.clean27$FATHER<- NULL
sparrow.xovers.clean27$MOTHER<- NULL
sparrow.xovers.clean27$UniqueID <- substr(sparrow.xovers.clean27$UniqueID, 5, nchar(sparrow.xovers.clean27$UniqueID))
sparrow.xovers.clean27$analysisID <- "Chr27"

sparrow.xovers.clean28$No.Inf.Loci <- NULL
sparrow.xovers.clean28$First.Inf.Order <- NULL
sparrow.xovers.clean28$Last.Inf.Order<- NULL
sparrow.xovers.clean28$FATHER<- NULL
sparrow.xovers.clean28$MOTHER<- NULL
sparrow.xovers.clean28$UniqueID <- substr(sparrow.xovers.clean28$UniqueID, 5, nchar(sparrow.xovers.clean28$UniqueID))
sparrow.xovers.clean28$analysisID <- "Chr28"












z <- merge(sparrow.xovers.clean1, sparrow.xovers.clean2, by=c("ANIMAL","parent","Family", "RRID"))
z$data.x <- NULL
z$analysisID.x <- NULL
z$UniqueID.x <- NULL
z$data.y <- NULL
z$analysisID.y <- NULL
z$UniqueID.y <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean3, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean4, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean5, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean6, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean7, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean8, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean9, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean10, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean11, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean12, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean13, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean14, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean15, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean17, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean18, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean19, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean20, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean21, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clea22, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean23, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean24, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean25, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean26, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean27, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean28, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL
z <- merge(z, sparrow.xovers.clean29, by=c("ANIMAL","parent","Family", "RRID"))
z$data <- NULL
z$analysisID <- NULL
z$UniqueID <- NULL
z$RecombCount <- z$RecombCount.x + z$RecombCount.y
z$RecombCount.x <- NULL
z$RecombCount.y <- NULL

xovers.plot.all <- z %>%
  ggplot(aes(x = RecombCount), xlim = 5000)+
  geom_histogram()+
  ggtitle("Xovers.Clean.AllChrs")
xovers.plot.all

badids_HighRecomb <- unique(z$RRID[which(z$RecombCount > 200)])
write.table(badids_HighRecomb, file = "badids_HiRecomb.txt", row.names = FALSE, col.names = FALSE)
