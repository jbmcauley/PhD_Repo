rm(list=ls())

#Fix Sex in GWAA.data file ----
library(GenABEL)
library(crimaptools)
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/")
setwd("C:/Users/johnb/Dropbox/PLINK-files 200k SNP-data/")
dir()

load("sparrowgen_Helgeland_01_2018.RData")
sexings_table <-read.table("All_sexings_Helgeland_200k_29112017.txt", header = TRUE)

#In genABEL the heterogametic sex is coded as 1 (Female WZ in birds...); Homogametic is 0 (Male, ZZ in birds...)
#Adding proper sex allocation to genABEL data

sexings_table <- sexings_table[, c(1,2,4)]
names(sexings_table) <- c("id", "sex", "plinksex")
sex.add <- merge(sparrowgen.Helgeland@phdata, sexings_table, by = "id")
sex.add$sex <- sex.add$sex.y

sex.add$sex <- as.character(sex.add$sex)
sex.add$sex[sex.add$sex == "m"] <- 1   #May need to recode incase swapping was unneccesarry 
sex.add$sex[sex.add$sex == "f"] <- 0   #May need recode if swap unneccessary
sex.add$plinksex <- as.character(sex.add$plinksex)
sex.add$plinksex[sex.add$plinksex == "m"] <- 1   #May need to recode incase swapping was unneccesarry 
sex.add$plinksex[sex.add$plinksex == "f"] <- 2   #May need recode if swap unneccessary
sex.add$plinksex[sex.add$plinksex == "u"] <- 0   #May need recode if swap unneccessary

sex.add$sex.x <- NULL
sex.add$sex.y <- NULL

names(sex.add) <- c("id", "father", "mother", "plinksex","sex")

#Can use this to remove ind with missing sex
#sex.add <- sex.add[complete.cases(sex.add[,4]),]
sex.add$sex <- as.numeric(sex.add$sex)
sex.add$sex <- as.integer(sex.add$sex)
sex.add$plinksex <- as.numeric(sex.add$plinksex)
sex.add$plinksex <- as.integer(sex.add$plinksex)

sparrowgen.Helgeland@phdata[,4] <- sex.add$sex
sparrowgen.Helgeland@phdata[,5] <- sex.add$plinksex
names(sparrowgen.Helgeland@phdata) <- c("id", "father", "mother", "sex","plinksex")

rm(sex.add)
rm(sexings_table)
sparrowgen.Helgeland@gtdata@male <- sparrowgen.Helgeland@phdata$sex

#setwd("C:/Users/johnb/Downloads/crimaptools-master/crimaptools-master/data")----
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/")


library(crimaptools)
#load("deer.RData")
library(reshape2)

op.pair <- read.table("updateParents.txt")
row.names(op.pair) <- NULL
op.pair <- as.data.frame(op.pair)
names(op.pair) <- c("Fam","id","Father","Mother")
op.pair$Fam <- NULL

#-------------------------------------------------------------------------------

names(op.pair) <- c("ANIMAL", "FATHER", "MOTHER")
op.pair$ANIMAL <- as.numeric(as.character(op.pair$ANIMAL))
op.pair$FATHER <- as.numeric(as.character(op.pair$FATHER))
op.pair$MOTHER <- as.numeric(as.character(op.pair$MOTHER))
library(reshape2)

#Fix the following code after QC's have been applied to the sparrow.abel file. This needs to be done 
#because Crimap will get upset if the family pedigree contains IDs not in the gwaa.data object.
load('sparrowABEL_QC.RData')
sparrow.abel <- data1
rm(data1)
id.sub <- sparrow.abel@phdata$id
id.sub <- as.numeric(id.sub)
which(op.pair$ANIMAL %in% id.sub)
op.pair <- op.pair[which(op.pair$ANIMAL %in% id.sub),]
row.names(op.pair) <- NULL
####
####
####

op.pair.zeros <- melt(op.pair, id = "ANIMAL")
row.names(op.pair.zeros) <- NULL
GPa <- op.pair.zeros[1:length(op.pair$ANIMAL),]
GMa <- op.pair.zeros[(length(op.pair$ANIMAL)+1):length(op.pair.zeros$ANIMAL),]



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



#-----Maternal GrandMother Setup----
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


#-----Paternal GrandFather Setup----

Grand <- list()
counter <- 0
for (i in 1:length(op.pair$ANIMAL)) {
  counter <- counter +1
  if(op.pair$FATHER[i] != 0){
    Grand[[counter]] <- op.pair$FATHER[which(op.pair$ANIMAL == op.pair$FATHER[i])]
  }else{Grand[[counter]] <- 0} #Ind not a mother therefore the maternal Grandfather is unknown
}

idx <- !(sapply(Grand, length))
Grand[idx] <- 0
Grand <- do.call(rbind, Grand)
GPa$PGrandF <- Grand[,1]

#Paternal GrandMother Setup-
PGrandM <- list()
counter <- 0
for (i in 1:length(op.pair$ANIMAL)) {
  counter <- counter + 1
  PGrandM[[counter]] <- GMa$value[which(GPa$ANIMAL == GPa$value[i])]
}
idx <- !(sapply(PGrandM, length))
PGrandM[idx] <- 0
PGrandM <- do.call(rbind, PGrandM)
GPa$PGrandM <- PGrandM[,1]






#-----Rebuild .ped-----
op.pair$MGrandM <- MGrand[,1]
op.pair$MGrandF <- PGrand[,1]


df <- op.pair
df$PGrandM <- NULL
df$PGrandF <- NULL

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
sparrow.ped <- df
sparrow.famped <- sparrow.famped[order(sparrow.famped$Family),]
row.names(sparrow.famped) <- NULL

#----Making Father and Grandparent-Parents 0s----
counter <- 2
for (i in 1:3415) {
  
  sparrow.famped$FATHER[counter] <- 0
  counter <- counter + 5
}
counter <- 2
for (i in 1:1474) {
  
  sparrow.famped$MOTHER[counter] <- 0
  counter <- counter + 5
}
counter <- 4
for (i in 1:1474) {
  
  sparrow.famped$FATHER[counter] <- 0
  counter <- counter + 5
}
counter <- 4
for (i in 1:1474) {
  
  sparrow.famped$MOTHER[counter] <- 0
  counter <- counter + 5
}
counter <- 5
for (i in 1:1474) {
  
  sparrow.famped$MOTHER[counter] <- 0
  counter <- counter + 5
}
counter <- 5
for (i in 1:1474) {
  
  sparrow.famped$FATHER[counter] <- 0
  counter <- counter + 5
}
sparrow.famped.ALL_Mothers <- sparrow.famped


##############
##############
##############
#Making a fam.ped for the Dads. 

op.pair$MGrandM <- PGrandM[,1]
op.pair$MGrandF <- Grand[,1]
op.pair$PGrandF <- NULL
op.pair$PGrandM <- NULL

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
  FAM[counter] <- paste("Offspring_Dad", as.character(df$ANIMAL[i]), sep = "_")
}
for(i in 1:length(df$FATHER)){
  counter <- counter + 1
  FAM[counter] <- paste("Offspring_Dad", as.character(df$ANIMAL[i]), sep = "_")
}
for(i in 1:length(df$MOTHER)){
  counter <- counter + 1
  FAM[counter] <- paste("Offspring_Dad", as.character(df$ANIMAL[i]), sep = "_")
}
for(i in 1:length(df$MGrandM)){
  counter <- counter + 1
  FAM[counter] <- paste("Offspring_Dad", as.character(df$ANIMAL[i]), sep = "_")
}
for(i in 1:length(df$MGrandF)){
  counter <- counter + 1
  FAM[counter] <- paste("Offspring_Dad", as.character(df$ANIMAL[i]), sep = "_")
}
pedvec$Family <- FAM
length(unique(pedvec$Family))

sparrow.famped.Fathers <- pedvec
sparrow.ped <- df
sparrow.famped.Fathers <- sparrow.famped.Fathers[order(sparrow.famped.Fathers$Family),]
row.names(sparrow.famped.Fathers) <- NULL


#----Making Mother and Grandparent-Parents 0s----
counter <- 3
for (i in 1:3415) {
  
  sparrow.famped.Fathers$FATHER[counter] <- 0
  counter <- counter + 5
}
counter <- 3
for (i in 1:1474) {
  
  sparrow.famped.Fathers$MOTHER[counter] <- 0
  counter <- counter + 5
}
counter <- 4
for (i in 1:1474) {
  
  sparrow.famped.Fathers$FATHER[counter] <- 0
  counter <- counter + 5
}
counter <- 4
for (i in 1:1474) {
  
  sparrow.famped.Fathers$MOTHER[counter] <- 0
  counter <- counter + 5
}
counter <- 5
for (i in 1:1474) {
  
  sparrow.famped.Fathers$MOTHER[counter] <- 0
  counter <- counter + 5
}
counter <- 5
for (i in 1:1474) {
  
  sparrow.famped.Fathers$FATHER[counter] <- 0
  counter <- counter + 5
}




sparrow.famped <- rbind(sparrow.famped.ALL_Mothers,sparrow.famped.Fathers)
write.table(sparrow.famped,file = "fam.ped_unabridged.txt", row.names = FALSE)




#Tidying----
rm(df)
rm(GPa)
rm(GMa)
rm(PGrand)
rm(MGrand)
rm(deer.abel)
rm(deer.ped)
rm(deer.famped)
sparrow.ped$MGrandM <- NULL
sparrow.ped$MGrandF <- NULL
rm(list = ls.str(mode = 'numeric'))
rm(list = ls.str(mode = 'logical'))
rm(list = ls.str(mode = 'character'))


####
####
####
####


#Attempting CryMap with previously created Fam Var----
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/")
library(crimaptools)
library(GenABEL)

#ensure pheno file has proper sex assignments.
newfile.pheno <- read.table("updateSex.txt", header = FALSE)
newfile.pheno$V1 <- NULL
names(newfile.pheno) <- c("id", "sex")
newfile.pheno$sex[newfile.pheno$sex == 2] <- 0
write.table(newfile.pheno, file = "newfile.pheno", row.names = FALSE)

#Use pheno file and most recent .gen file to construct genABEL file
sparrow.abel <- load.gwaa.data(phe = "newfile.pheno", gen = "newfile-ME-fixed.gen")
sparrow.abel@phdata$father <- as.character(op.pair$FATHER)
sparrow.abel@phdata$mother <- as.character(op.pair$MOTHER)
sparrow.abel@phdata$id <- sparrow.abel@gtdata@idnames

x <- as.data.frame(sparrow.abel@phdata$id)
x$gtdata <- sparrow.abel@gtdata@idnames

#rm(sparrowgen.Helgeland)
#rm(op.pair)
rm(op.pair.zeros)
getwd()
sparrow.famped_backup <- sparrow.famped
#sparrow.famped <- sparrow.famped_backup
sparrow.ped_backup <- sparrow.ped
#sparrow.ped <- sparrow.ped_backup
vec <- c(0,0,0)
sparrow.ped <- rbind(sparrow.ped,vec)
sparrow.famped <- sparrow.famped[which(!(sparrow.famped$FATHER %in% sparrow.famped$ANIMAL)),]
sparrow.famped <- sparrow.famped[which(sparrow.famped$MOTHER %in% sparrow.ped$ANIMAL),]
sparrow.famped <- sparrow.famped[which(sparrow.famped$ANIMAL %in% sparrow.ped$ANIMAL),]


#Remove incomplete families of less than 5
sparrow.famped$Family <- as.factor(sparrow.famped$Family)
sparrow.famped <- sparrow.famped[!(as.numeric(sparrow.famped$Family) %in% which(table(sparrow.famped$Family)<5)),]
sparrow.famped$Family <- as.character(sparrow.famped$Family)


library(dplyr)
sparrow.famped$ANIMAL <- as.character(sparrow.famped$ANIMAL)
sparrow.famped$FATHER <- as.character(sparrow.famped$FATHER)
sparrow.famped$MOTHER <- as.character(sparrow.famped$MOTHER)
sparrow.famped$Family <- as.character(sparrow.famped$Family)
sparrow.famped <- data.frame(lapply(sparrow.famped, as.character), stringsAsFactors=FALSE)
#sparrow.famped <-z

#Removing following families: 537, 663, 665, 706, 764, 770, 781, 805, 819, 821, 835, 836, 856, 870, 878, 885, 940, 1000, 1119,
                          #   1120, 1121, 1940, 1982, 1202, 1203, 1204, 1512, 1556, 1629, 1704, 1715, 1782, 1940, 1982, 1993, 2030, 2201, 2203,
                          #   2231, 2253, 2262, 2292, 2340, 2350, 2352, 2357, 2366, 2380, 2396, 2423, 2429, 2434, 2440, 2448, 2486,
                          #   2493, 2509, 2615, 2695, 2731, 2732, 2784, 2785, 2789, 2946, 2978, 2979, 3135, 3136, 3149, 3156, 3158,
                          #   3241, 3314, 3316, 3342, 3343, 3344, 3345, 3397, 3424, 3426, 3495)
#consider devopling script to elimnate manual family type?
sparrow.famped <- sparrow.famped[-which(sparrow.famped$Family == "Offspring_Mum_1982"),]


sparrow.famped$ANIMAL <- as.numeric(sparrow.famped$ANIMAL)
sparrow.famped$FATHER <- as.numeric(sparrow.famped$FATHER)
sparrow.famped$MOTHER <- as.numeric(sparrow.famped$MOTHER)

create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "19a", 
                    chr = 19, 
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)



run_crimap_prepare(genfile = "crimap/chr19a.gen", crimap.path = "C:/PathApps/crimap.exe")
#Function has not produced the .loc .par and .dat files suggested in the tutorial
#Using the terminal does work following the crimapinput1 preiously created (n,n,n,n,7,y,y) to generate .loc, .par, and .dat files... mention to Susan a potential issue in crimaptools
#Alternatively the problem is due to miscommunication between R and Path, have had Path issues in past with Plink...

dir("crimap")

parse_mend_err(prefile = "crimap/chr19a.pre", genfile = "crimap/chr19a.gen", familyPedigree = sparrow.famped)
read.table("crimap/chr19a.mnd", header = T)
#No Mendelian errors detected, in tutorial "2" are listed... Data already cleaned?
#Next line only used if mendialian erros are present in order to mask them in the .gen file.
create_crimap_input(gwaa.data = sparrow.abel, 
                     familyPedigree = sparrow.famped, 
                     analysisID = "19a", 
                     chr = 19, 
                     outdir = "crimap", 
                     clear.existing.analysisID = TRUE, 
                     use.mnd = TRUE)


#If mendelian errors detected you will need to rerun the prepare function. Once again, might need to be done in terminal rather than crimap-tools
run_crimap_prepare(genfile = "crimap/chr19a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr19a.pre", genfile = "crimap/chr19a.gen", familyPedigree = sparrow.famped)
dir("crimap")
#Repeat the process of looking for mendelian error until there are none.


#Build a Linkage Map:----

run_crimap_map(genfile = "crimap/chr19a.gen", crimap.path = "C:/PathApps/crimap.exe")
#.map file generated has size "0 B" This is a sign that something is not working in generation of the file. Will need to work with regular functions 
#in terminal.
dir("crimap")
sparrow.map19 <- parse_map(mapfile = "crimap/chr19a.map")
#Produces error, might be an issue with the sex assignments in the sparrow.abel gwaa.data file in some fashion. Needs checking.


#Characterizing recombination events:----

run_crimap_chrompic(genfile = "crimap/chr19a.gen", crimap.path =  "C:/PathApps/crimap.exe")
sparrow.cmpmap <- parse_map_chrompic(chrompicfile = "crimap/chr19a.cmp")
head(sparrow.cmpmap)

sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr19a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers[1:2,]


#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

head(sparrow.doubles)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 9], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 9], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 9)), 
                      analysisID = "9a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)

sparrow.xovers.clean




