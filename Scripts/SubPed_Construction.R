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

y <- list()
counter <- 0
if(any((deer.ped$ANIMAL %in% deer.ped$MOTHER)|(deer.ped$ANIMAL %in% deer.ped$FATHER))){
  mothers <- which(deer.ped$ANIMAL %in% deer.ped$MOTHER)
  fathers <- which(deer.ped$ANIMAL %in% deer.ped$FATHER)
  idfathers <- deer.ped$ANIMAL[fathers]
  idmothers <- deer.ped$ANIMAL[mothers]
  
  
  deer.ped$FAMILY <-
}

deer.ped$ANIMAL[deer.ped$FATHER %in% idfathers]


deer.ped$ANIMAL[nonfathers]
deer.ped$FAMILY <- NULL

y <- do.call(rbind, y)

#==========================
ped <- deer.ped
if(any(!ped$MOTHER %in% ped$ANIMAL)){
  ped <- rbind(data.frame(ANIMAL = ped$MOTHER[which(!ped$MOTHER %in% ped$ANIMAL)],
                          MOTHER = 0, FATHER = 0,
                          stringsAsFactors = F),
               ped)
}

if(any(!ped$FATHER %in% ped$ANIMAL)){
  ped <- rbind(data.frame(ANIMAL = ped$FATHER[which(!ped$FATHER %in% ped$ANIMAL)],
                          MOTHER = 0, FATHER = 0,
                          stringsAsFactors = F),
               ped)
}


z<-rbind(deer.ped,
data.frame(ANIMAL = ped$MOTHER[which(ped$MOTHER %in% ped$ANIMAL)],
            MOTHER = 0, FATHER = 0,
            stringsAsFactors = F),

data.frame(ANIMAL = ped$FATHER[which(ped$FATHER %in% ped$ANIMAL)],
           MOTHER = 0, FATHER = 0,
           stringsAsFactors = F)
)


