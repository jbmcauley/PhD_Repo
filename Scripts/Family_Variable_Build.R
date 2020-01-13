setwd("C:/Users/johnb/Downloads/crimaptools-master/crimaptools-master/data")
library(crimaptools)
load("deer.RData")
library(reshape2)
op.pair.zeros <- melt(deer.ped, id = "ANIMAL")

op.pair.zeros <- op.pair.zeros[op.pair.zeros[,3] != 0,]
row.names(op.pair.zeros) <- NULL

GPa <- op.pair.zeros[1:14,]
GMa <- op.pair.zeros[15:28,]


#-----Grandfather Setup----
PGrand <- list()
counter <- 0
for (i in 1:14) {
  counter <- counter +1
  PGrand[[counter]] <- GPa$value[which(GPa$ANIMAL ==GPa$value[i])]
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
deer.ped$MGrand <- MGrand
deer.ped$PGrand <- PGrand


#Long Format
op.pair.zeros <- melt(deer.ped, id = "ANIMAL")
op.pair.zeros <- op.pair.zeros[op.pair.zeros[,3] != 0,]
row.names(op.pair.zeros) <- NULL

#Wide Format
op.pair <- dcast(op.pair.zeros, ANIMAL ~ variable, value.var = "value")

op.pair$Family.Interaction <- interaction(op.pair$MOTHER, op.pair$FATHER)

grep(op.pair$ANIMAL[12], op.pair$MOTHER)


#-----Build Fam Group: Off-Parent----
looplist <- list()
counter <- 0
for(i in 1:14){
  #If ind i is a parent
  if(any(op.pair$FATHER == op.pair$ANIMAL[i]) | any(op.pair$MOTHER == op.pair$ANIMAL[i])){
    #Dads
    if(any(op.pair$FATHER == op.pair$ANIMAL[i])){
      
    }
    #Moms (Exactly the same as for Dads but alternate column)
    else{
      Offspringvec <- op.pair$ANIMAL[which(op.pair$MOTHER == op.pair$ANIMAL[i])]
      if(any(Offspringvec))
    }
  }
  #Not a parent: Simply assign fam var once
  else{
    counter <- counter + 1
    looplist[[counter]] <- paste("Offspring_Mum", as.character(op.pair$ANIMAL[i]), sep = "_")
  }

#Close "for" loop  
}


