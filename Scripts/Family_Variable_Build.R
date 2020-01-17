rm(list = ls())
#setwd("C:/Users/johnb/Downloads/crimaptools-master/crimaptools-master/data")
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/")
load("sparrowgen_Helgeland_01_2018.RData")

library(crimaptools)
#load("deer.RData")
library(reshape2)

op.pair <- sparrowgen.Helgeland@phdata
row.names(op.pair) <- NULL
op.pair$sex <- NULL
op.pair$id.og <- op.pair$id
op.pair$Father.og <- op.pair$Father
op.pair$Mother.og <- op.pair$Mother

 
op.pair$id <- gsub("M", "1", op.pair$id.og)
op.pair$Father <- gsub("M", "1", op.pair$Father.og)
op.pair$Mother <- gsub("M", "1", op.pair$Mother.og)

op.pair$id <- gsub("N", "2", op.pair$id)
op.pair$Father <- gsub("N", "2", op.pair$Father)
op.pair$Mother <- gsub("N", "2", op.pair$Mother)

op.pair$id <- gsub("L", "3", op.pair$id)
op.pair$Father <- gsub("L", "3", op.pair$Father)
op.pair$Mother <- gsub("L", "3", op.pair$Mother)

op.pair$id <- gsub("F", "4", op.pair$id)
op.pair$Father <- gsub("F", "4", op.pair$Father)
op.pair$Mother <- gsub("F", "4", op.pair$Mother)

op.pair$id.og <- NULL
op.pair$Father.og <- NULL
op.pair$Mother.og <- NULL

names(op.pair) <- c("ANIMAL", "FATHER", "MOTHER")
op.pair$ANIMAL <- as.numeric(as.character(op.pair$ANIMAL))
op.pair$FATHER <- as.numeric(as.character(op.pair$FATHER))
op.pair$MOTHER <- as.numeric(as.character(op.pair$MOTHER))

op.pair.zeros <- melt(op.pair, id = "ANIMAL")
row.names(op.pair.zeros) <- NULL
GPa <- op.pair.zeros[1:3116,]
GMa <- op.pair.zeros[3117:6232,]


#-----Maternal Grandfather Setup----
PGrand <- list()
counter <- 0
for (i in 1:length(op.pair$ANIMAL)) {
  counter <- counter +1
  if(op.pair$MOTHER[i] != 0){
    PGrand[[counter]] <- op.pair$FATHER[which(op.pair$ANIMAL == op.pair$MOTHER[i])]
  }else{PGrand[[counter]] <- 0} #Ind not a mother therefore the maternal Grandfather is unknown
}

idx <- !(sapply(PGrand, length))
PGrand[idx] <- 0
PGrand <- do.call(rbind, PGrand)
GMa$PGrand <- PGrand[,1]



#-----GrandMother Setup----
MGrand <- list()
counter <- 0
for (i in 1:length(op.pair$ANIMAL)) {
  counter <- counter + 1
  MGrand[[counter]] <- GMa$value[which(GMa$ANIMAL == GMa$value[i])]
}
idx <- !(sapply(MGrand, length))
MGrand[idx] <- 0
MGrand <- do.call(rbind, MGrand)
GMa$MGrand <- MGrand[,1]


#-----Rebuild .ped-----
op.pair$MGrandM <- MGrand[,1]
op.pair$MGrandF <- PGrand[,1]
df <- op.pair
row_sub <- apply(op.pair, 1, function(row) all(row !=0 ))
df <- df[row_sub,]
row.names(df) <- NULL

OffspringVec <- df$ANIMAL
faths <- df$FATHER
mother <- df$MOTHER
Mgrandf <- df$MGrandF
MgrandM <- df$MGrandM
pedvec <- as.data.frame(c(OffspringVec,faths,mother, Mgrandf,MgrandM))
names(pedvec) <- "ANIMAL"

FATHER <- vector()
MOTHER <- vector()
counter <- 0
for(i in 1:length(pedvec$ANIMAL)){
  counter <- counter + 1
  if(any(pedvec$ANIMAL[i] == op.pair$ANIMAL)){
  FATHER[counter] <- op.pair$FATHER[which(pedvec$ANIMAL[i] == op.pair$ANIMAL)]
  MOTHER[counter] <- op.pair$MOTHER[which(pedvec$ANIMAL[i] == op.pair$ANIMAL)]}
  else{
    FATHER[counter] <- 0
    MOTHER[counter] <- 0
  }
  }
pedvec$FATHER <- FATHER
pedvec$MOTHER <- MOTHER


FAM <- vector()
counter <- 0
for(i in 1:length(df$ANIMAL)){
  counter <- counter + 1
  FAM[counter] <- paste("Offspring_Mum", as.character(df$ANIMAL[i]), sep = "_")
}
for(i in 1:length(df$FATHER)){
  counter <- counter + 1
  FAM[counter] <- paste("Offspring_Mum", as.character(df$ANIMAL[i]), sep = "_")
  }
for(i in 1:length(df$MOTHER)){
  counter <- counter + 1
FAM[counter] <- paste("Offspring_Mum", as.character(df$ANIMAL[i]), sep = "_")
}
for(i in 1:length(df$MGrandM)){
  counter <- counter + 1
  FAM[counter] <- paste("Offspring_Mum", as.character(df$ANIMAL[i]), sep = "_")
}
for(i in 1:length(df$MGrandF)){
  counter <- counter + 1
FAM[counter] <- paste("Offspring_Mum", as.character(df$ANIMAL[i]), sep = "_")
}
pedvec$Family <- FAM
length(unique(pedvec$Family))

sparrow.famped <- pedvec
