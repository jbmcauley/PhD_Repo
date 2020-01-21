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
sparrowgen.Helgeland@gtdata@idnames<- as.character(as.numeric(op.pair$id))
names(sparrowgen.Helgeland@gtdata@male) <- as.character(as.numeric(op.pair$id))


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
sparrow.ped <- df

#Tidying----
rm(df)
rm(pedvec)
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


#Attempting CryMap with previously created Fam Var----
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/")
library(crimaptools)
library(GenABEL)
sparrow.abel <- sparrowgen.Helgeland
sparrow.abel@phdata$Father <- as.character(op.pair$FATHER)
sparrow.abel@phdata$Mother <- as.character(op.pair$MOTHER)


#rm(sparrowgen.Helgeland)
#rm(op.pair)
rm(op.pair.zeros)
getwd()
sparrow.famped_backup <- sparrow.famped
sparrow.ped_backup <- sparrow.ped
vec <- c(0,0,0)
sparrow.ped <- rbind(sparrow.ped,vec)
sparrow.famped <- sparrow.famped[which(sparrow.famped$FATHER %in% sparrow.ped$ANIMAL),]
sparrow.famped <- sparrow.famped[which(sparrow.famped$MOTHER %in% sparrow.ped$ANIMAL),]
sparrow.famped <- sparrow.famped[which(sparrow.famped$ANIMAL %in% sparrow.ped$ANIMAL),]



create_crimap_input(gwaa.data = sparrow.abel, familyPedigree = sparrow.famped, analysisID = "8a", chr = 8, outdir = "crimap", clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr8a.gen", crimap.path = "C:/PathApps/crimap.exe")
#Function has not produced the .loc .par and .dat files suggested in the tutorial
#Using the terminal does work following the crimapinput1 preiously created (n,n,n,n,7,y,y) to generate .loc, .par, and .dat files... mention to Susan a potential issue in crimaptools
#Alternatively the problem is due to miscommunication between R and Path, have had Path issues in past with Plink...

dir("crimap")

parse_mend_err(prefile = "crimap/chr8a.pre", genfile = "crimap/chr8a.gen", familyPedigree = sparrow.famped)
read.table("crimap/chr8a.mnd", header = T)
#No Mendelian errors detected, in tutorial "2" are listed... Data already cleaned?
#Next line only used if mendialian erros are present in order to mask them in the .gen file.
create_crimap_input (gwaa.data = sparrow.abel, 
                     familyPedigree = sparrow.famped, 
                     analysisID = "8a", 
                     chr = 8, 
                     outdir = "crimap", 
                     clear.existing.analysisID = TRUE, 
                     use.mnd = TRUE)


#If mendelian errors detected you will need to rerun the prepare function. Once again, might need to be done in terminal rather than crimap-tools
run_crimap_prepare(genfile = "crimap/chr8a.gen", crimap.path = "C:/PathApps/crimap.exe") 
parse_mend_err(prefile = "crimap/chr8a.pre", genfile = "crimap/chr8a.gen", familyPedigree = sparrow.famped)
dir("crimap")
#Repeat the process of looking for mendelian error until there are none.


#Build a Linkage Map:----

run_crimap_map(genfile = "crimap/chr8a.gen", crimap.path = "C:/PathApps/crimap.exe")
#.map file generated has size "0 B" This is a sign that something is not working in generation of the file. Will need to work with regular functions 
#in terminal.
dir("crimap")
sparrow.map <- parse_map(mapfile = "crimap/chr8a.map")
#Produces error, might be an issue with the sex assignments in the sparrow.abel gwaa.data file in some fashion. Needs checking.


#Characterizing recombination events:----

run_crimap_chrompic(genfile = "crimap/chr8a.gen", crimap.path =  "C:/PathApps/crimap.exe")
sparrow.cmpmap <- parse_map_chrompic(chrompicfile = "crimap/chr8a.cmp")
head(sparrow.cmpmap)

sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr8a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers[1:2,]


#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

head(sparrow.doubles)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 3], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 3], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 3)), 
                      analysisID = "3a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)

sparrow.xovers.clean




