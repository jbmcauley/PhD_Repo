setwd("C:/Users/johnb/Downloads/crimaptools-master/crimaptools-master/data")
library(crimaptools)
load("deer.RData")
library(reshape2)
op.pair.zeros <- melt(deer.ped, id = "ANIMAL")

op.pair.zeros <- op.pair.zeros[op.pair.zeros[,3] != 0,]
row.names(op.pair.zeros) <- NULL

GPa <- op.pair.zeros[1:14,]
GMa <- op.pair.zeros[15:28,]
rm(v)
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
df <- dcast(op.pair.zeros, ANIMAL ~ variable, value.var = "value")
df<- df[complete.cases(df),]

OffspringVec <- df$ANIMAL
faths <- df$FATHER
mother <- df$MOTHER
Mgrandf <- df$MGrandF
MgrandM <- df$MGrandM
pedvec <- as.data.frame(c(OffspringVec, df$FATHER, df$MOTHER, df$MGrandM, df$MGrandF))
names(pedvec) <- "ANIMAL"

FATHER <- vector()
MOTHER <- vector()
FAM <- vector()
counter <- 0
for(i in 1:length(pedvec$ANIMAL)){
  counter <- counter + 1
  FATHER[counter] <- op.pair$FATHER[which(pedvec$ANIMAL[i] == op.pair$ANIMAL)]
  MOTHER[counter] <- op.pair$MOTHER[which(pedvec$ANIMAL[i] == op.pair$ANIMAL)]
}

FAM[counter] <- paste("Offspring_Mum", as.character(op.pair$ANIMAL[i]), sep = "_")

pedvec$MOTHER
pedvec$Family
paste("Offspring_Mum", as.character(op.pair$ANIMAL[i]), sep = "_")


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
      
      #If ind is a GPa
      if(any(op.pair$ANIMAL[op.pair$PGrand == op.pair$ANIMAL[i]])){
       OffspringVector  <- op.pair$ANIMAL[op.pair$PGrand == op.pair$ANIMAL[i]] #vector of offspring *Offspring may be parents!

       #Create row for each offspring/GParent Pairing ***This script does not check if offspring of Gparent are a parent themselves.     
       for(j in 1:length(OffspringVector)){
         counter <- counter + 1
         tempdf <- data.frame(op.pair$ANIMAL[i], op.pair$FATHER[i], op.pair$MOTHER[i], paste("Offspring_Mum", as.character(OffspringVector[j]), sep = "_"))
         names(tempdf) <- c(names(op.pair)[c(1,2,3)],"Family")
         looplist[[counter]] <- tempdf
      }
      }
    
      #Non-Grandparental - Dads
      else{
        OffspringVector  <- op.pair$ANIMAL[op.pair$FATHER == op.pair$ANIMAL[i]] #vector of offspring *Offspring may be parents!
        if(anyOffspringVector == op.pair$MOTHER)
      }
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

#Close of Starting "for" loop  
}








