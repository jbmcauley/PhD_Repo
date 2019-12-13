rm(list=ls())



setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/Ped_&_map_file_formatting")
dir()


Helg_Ped <- read.table("SNP_pedigree_Helgeland_05122017.txt", header = TRUE)

Helg_Ped <- Helg_Ped[complete.cases(Helg_Ped),]
row.names(Helg_Ped) <- NULL

setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/")
sexings_table <-read.table("All_sexings_Helgeland_200k_29112017.txt", header = TRUE)
#sexings_table <- sexings_table[complete.cases(sexings_table[,2]),]
#sexings_table$GenABEL <- ifelse(sexings_table$Adult_phenotypic == "m", 1, 2)



sexings_table <- sexings_table[, c(1,2)]
names(sexings_table) <- c("id", "sex")
Ped <- merge(Helg_Ped, sexings_table, by = "id")
Ped$sex <- as.character(Ped$sex)
Ped$sex[Ped$sex == "m"] <- 1
Ped$sex[Ped$sex == "f"] <- 2
names(Ped) <- c("id", "mother", "father", "sex")
#Ped <- Ped[complete.cases(Ped[,4]),]
Ped$sex <- as.numeric(Ped$sex)


#========================= Grouping Siblings:====

fam = apply(Ped[2:3], 1, function(x) paste0(sort(x), collapse=" "))

temp_1 = split(Ped, fam)
names(temp_1)[1] = "Dbl NAs"

temp_2 = lapply(1:length(temp_1),
                function(x) temp_1[[1]][which(temp_1[[1]]$Person %in% 
                                                unique(unlist(temp_1[[x]][,c(2,3)], use.names=FALSE))),])

OUT = lapply(1:length(temp_1), function(x) rbind(temp_2[[x]], temp_1[[x]]))
names(OUT) = names(temp_1)

OUT = do.call("rbind", 
              lapply(1:length(OUT), 
                     function(x) cbind(OUT[[x]], fam = names(OUT[x]))))

OUT <- OUT[,c(1,5)]

Ped_wFam <- merge(Ped, OUT, by = "id")


#========================= Building a Family with Susan's data:=====

setwd("C:/Users/s1945757/PhD_Repo/Cri_Map/crimaptools-master/crimaptools-master/data")
load("deer.RData")

#For each ind in x: search the father column, if match (if no match search mother), then grab ANIMAL value and 
#search parental coloumns, repeat until no match, once no match assign current animal value to family vector
#

x<-unique(cbind(c(unique(deer.ped$ANIMAL), unique(deer.ped$FATHER), unique(deer.ped$MOTHER))))
x <- x[x!= 0]
x <- as.data.frame(x)
row.names(deer.ped) <- NULL
FAMILY <- list()
counter <- 0

if(any((deer.ped$ANIMAL %in% deer.ped$MOTHER)|(deer.ped$ANIMAL %in% deer.ped$FATHER))){
  for(i in 1:14){
    if(deer.ped$MOTHER[i] != 0){                                                                           #condition 1; ind has a mother
      if(deer.ped$MOTHER[deer.ped$ANIMAL == deer.ped$MOTHER[i]] != 0){                                     #condition 2; ind's mother has a mother
       if(deer.ped$MOTHER[deer.ped$ANIMAL ==deer.ped$MOTHER[deer.ped$ANIMAL == deer.ped$MOTHER[i]]] != 0){ #condition 3; ind's grandmother has a mother
      print("Another Line Needed")
         break } # <- insert code for If all 3 condtions pass
        counter <- counter + 1
        focalind <- deer.ped[i,1:3]
        FAMILY[[counter]] <- c(focalind,  } # <- code for if 1&2 pass, but 3 fails; assigns grandchild as family code
            if(any(deer.ped$ANIMAL[i] %in% deer.ped$MOTHER)){
              deer.ped$FAMILY[i]
            }} #<- if 1 passes, but 2&3 fail
    else {
        }
          } #close loop
            } #Close original if
  
  

#========================== Father/Mother Pairs
19:22
25:28
op.pair.zeros <- melt(deer.ped, id = "ANIMAL")
#op.pair.zeros <- unique(op.pair.zeros)
#op.pair <- op.pair.zeros[op.pair.zeros$value != 0,]
#op.pair.zeros <- op.pair.zeros[!duplicated(t(apply(op.pair.zeros[c("ANIMAL","value")],1,sort))),]
row.names(op.pair.zeros) <- NULL
names(op.pair.zeros) <- c("FID", "Parent", "PID")

MOTHERS <- vector()
FATHERS <- vector()
for(i in 1:28){
  MOTHERS[i] <- deer.ped$MOTHER[deer.ped$ANIMAL == op.pair.zeros[i,1]]
  FATHERS[i] <- deer.ped$FATHER[deer.ped$ANIMAL == op.pair.zeros[i,1]]
}
op.pair.zeros$FATHER <- FATHERS
op.pair.zeros$MOTHER <- MOTHERS
op.pair.zeros$Parent <- NULL
op.pair.zeros$PID <- NULL

table(deer.famped$ANIMAL)
table(op.pair.zeros$ANIMAL)

GRANDMOTHERS <- vector()
GRANDFATHERS <- vector()
for(i in 1:28){
  if(deer.ped$MOTHER[op.pair.zeros$FID == op.pair.zeros[i,1]][1] > 0){
  GRANDMOTHERS[i] <- deer.ped$MOTHER[deer.ped$ANIMAL == deer.ped$MOTHER[deer.ped$ANIMAL == op.pair.zeros[i,1]]]
  GRANDFATHERS[i] <- deer.ped$FATHER[deer.ped$ANIMAL == deer.ped$FATHER[deer.ped$ANIMAL == op.pair.zeros[i,1]]]
  } else {
  GRANDMOTHERS[i] <- 0
  GRANDFATHERS[i] <- 0}
}
op.pair.zeros$GRANDFATHER <- GRANDFATHERS
op.pair.zeros$GRANDMOTHER <- GRANDMOTHERS

