rm(list = ls())
setwd("C:/Users/johnb/Downloads/crimaptools-master/crimaptools-master/data")
library(crimaptools)
load("deer.RData")
library(reshape2)

op.pair <- deer.ped
op.pair.zeros <- melt(deer.ped, id = "ANIMAL")
row.names(op.pair.zeros) <- NULL
GPa <- op.pair.zeros[1:14,]
GMa <- op.pair.zeros[15:28,]

#-----Maternal Grandfather Setup----
PGrand <- list()
counter <- 0
for (i in 1:14) {
  counter <- counter +1
  if(op.pair$MOTHER[i] != 0){
    PGrand[[counter]] <- op.pair$FATHER[which(op.pair$ANIMAL == op.pair$MOTHER[i])]
  }else{PGrand[[counter]] <- 0} #Ind not a mother therefore the maternal Grandfather is unknown
}

idx <- !(sapply(PGrand, length))
PGrand[idx] <- 0
PGrand <- do.call(rbind, PGrand)
GMa$PGrand <- PGrand



#-----GrandMother Setup----
MGrand <- list()
counter <- 0
for (i in 1:14) {
  counter <- counter + 1
  MGrand[[counter]] <- GMa$value[which(GMa$ANIMAL == GMa$value[i])]
}
idx <- !(sapply(MGrand, length))
MGrand[idx] <- 0
MGrand <- do.call(rbind, MGrand)
GMa$MGrand <- MGrand

op.pair.zeros <- rbind(GPa, GMa)




#-----Rebuild .ped-----
deer.ped$MGrandM <- MGrand
deer.ped$MGrandF <- PGrand
df <- deer.ped
row_sub <- apply(deer.ped, 1, function(row) all(row !=0 ))
df <- df[row_sub,]

OffspringVec <- df$ANIMAL
faths <- df$FATHER
mother <- df$MOTHER
Mgrandf <- df$MGrandF
MgrandM <- df$MGrandM
pedvec <- as.data.frame(c(OffspringVec, df$FATHER, df$MOTHER, df$MGrandM, df$MGrandF))
names(pedvec) <- "ANIMAL"

FATHER <- vector()
MOTHER <- vector()
counter <- 0
for(i in 1:length(pedvec$ANIMAL)){
  counter <- counter + 1
  FATHER[counter] <- op.pair$FATHER[which(pedvec$ANIMAL[i] == op.pair$ANIMAL)]
  MOTHER[counter] <- op.pair$MOTHER[which(pedvec$ANIMAL[i] == op.pair$ANIMAL)]
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


#Long Format----
op.pair.zeros <- melt(deer.ped, id = "ANIMAL")
op.pair.zeros <- op.pair.zeros[op.pair.zeros[,3] != 0,]
row.names(op.pair.zeros) <- NULL

#Wide Format----
op.pair <- dcast(op.pair.zeros, ANIMAL ~ variable, value.var = "value")

op.pair$Family.Interaction <- interaction(op.pair$MOTHER, op.pair$FATHER)

grep(op.pair$ANIMAL[12], op.pair$MOTHER)











