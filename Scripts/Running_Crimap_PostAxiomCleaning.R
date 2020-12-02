
##Run a full map with new GenABEL Object##
library(GenABEL)
library(crimaptools)
library(dplyr)


#Load in new  GenABEL object
setwd("C:/Users/s1945757/PhD_Repo/Pdom200K_Sparrow_NTNU_96/clean_genotype_data")
load("Clean_Map_SNPs.RData")

#Load in Famped file
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/crimap/crimap/crimap")
sparrow.famped <- read.table("FamPed_20200414.txt", header = TRUE)

#Famped has ALL ind included, crimap requires the genABEL and family pedigree have all the same ind. So we remove the extra from the Famped object.
sparrow.famped <- sparrow.famped[which(sparrow.famped$ANIMAL %in% cleanspar@phdata$id),]
sparrow.famped <- sparrow.famped[sparrow.famped$Family %in% names(which(table(sparrow.famped$Family)>=5)),]
sparrow.famped$Family <- as.character(sparrow.famped$Family)
sparrow.famped$Family <- as.factor(sparrow.famped$Family)



##Perform basic GenABEL QC

qc1 <- check.marker(cleanspar, p.level = 0)

data1 <- cleanspar[qc1$idok, qc1$snpok]
# Way too much removal, continue without this step


#sparrow.abel <- data1
#or
sparrow.abel <- cleanspar

##Crimap tools shorthand function
Crymap <- function(zzz) {
  create_crimap_input(gwaa.data = sparrow.abel, 
                      familyPedigree = sparrow.famped, 
                      analysisID = paste(zzz, "a", sep = ""), 
                      chr = zzz, 
                      outdir = "crimap", 
                      clear.existing.analysisID = TRUE)
  
  run_crimap_prepare(genfile = paste("crimap/chr", zzz, "a.gen", sep = ""), crimap.path = "C:/PathApps/crimap.exe")
  
  dir("crimap")
  
  parse_mend_err(prefile = paste("crimap/chr", zzz, "a.pre", sep = ""), genfile = paste("crimap/chr", zzz ,"a.gen", sep=""), familyPedigree = sparrow.famped)
  
  create_crimap_input(gwaa.data = sparrow.abel, 
                      familyPedigree = sparrow.famped, 
                      analysisID = paste(zzz, "a", sep = ""), 
                      chr = zzz, 
                      outdir = "crimap", 
                      clear.existing.analysisID = TRUE, 
                      use.mnd = TRUE)
  
  
  #If mendelian errors detected you will need to rerun the prepare function. Once again, might need to be done in terminal rather than crimap-tools
  run_crimap_prepare(genfile = paste("crimap/chr", zzz, "a.gen", sep = ""), crimap.path = "C:/PathApps/crimap.exe") 
  parse_mend_err(prefile = paste("crimap/chr", zzz, "a.pre", sep = ""), genfile = "crimap/chr", zzz, "a.gen", familyPedigree = sparrow.famped)
  dir("crimap")
  #Repeat the process of looking for mendelian error until there are none.
  
  
  #Build a Linkage Map:----
  
  run_crimap_map(genfile = paste("crimap/chr",zzz,"a.gen", sep = ""), crimap.path = "C:/PathApps/crimap.exe")
  #.map file generated has size "0 B" This is a sign that something is not working in generation of the file. Will need to work with regular functions 
  #in terminal.
  dir("crimap")
  
  #Produces error, might be an issue with the sex assignments in the data1 gwaa.data file in some fashion. Needs checking.
  
  
  #Characterizing recombination events:----
  
  run_crimap_chrompic(genfile = paste("crimap/chr",zzz,"a.gen", sep = ""), crimap.path =  "C:/PathApps/crimap.exe")
  
}

##Choose output directory where crimap function will create a folder crimap/ and place all crimap files
setwd("C:/Users/s1945757/PhD_Repo/Crimap Outputs Post Axiom Cleaning")

Crymap(zzz= 9)
sparrow.map9 <- parse_map(mapfile = "crimap/chr9a.map")
sparrow.cmpmap9 <- parse_map_chrompic(chrompicfile = ("crimap/chr9a.cmp"))

#===
#===

###Further cleaning based on initial .pre file error warnings###

sparrow.famped  <-  sparrow.famped[!(sparrow.famped$Family == "ID_1_MOTHER_1092"),]
##Try again##
##Choose output directory where crimap function will create a folder crimap/ and place all crimap files
setwd("C:/Users/s1945757/PhD_Repo/Crimap Outputs Post Axiom Cleaning")
#Famped fixing again required
sparrow.famped <- sparrow.famped[which(sparrow.famped$ANIMAL %in% cleanspar@phdata$id),]
sparrow.famped <- sparrow.famped[sparrow.famped$Family %in% names(which(table(sparrow.famped$Family)>=5)),]
sparrow.famped$Family <- as.character(sparrow.famped$Family)
sparrow.famped$Family <- as.factor(sparrow.famped$Family)
#Run Crimap
Crymap(zzz= 9)
sparrow.map9 <- parse_map(mapfile = "crimap/chr9a.map")
sparrow.cmpmap9 <- parse_map_chrompic(chrompicfile = ("crimap/chr9a.cmp"))

#===
#===

##New family flagged
sparrow.famped  <-  sparrow.famped[!(sparrow.famped$Family == "ID_1009_FATHER_617"),]
##Try again##
##Choose output directory where crimap function will create a folder crimap/ and place all crimap files
setwd("C:/Users/s1945757/PhD_Repo/Crimap Outputs Post Axiom Cleaning")
#Famped fixing again required
sparrow.famped <- sparrow.famped[which(sparrow.famped$ANIMAL %in% cleanspar@phdata$id),]
sparrow.famped <- sparrow.famped[sparrow.famped$Family %in% names(which(table(sparrow.famped$Family)>=5)),]
sparrow.famped$Family <- as.character(sparrow.famped$Family)
sparrow.famped$Family <- as.factor(sparrow.famped$Family)
#Run Crimap
Crymap(zzz= 9)
sparrow.map9 <- parse_map(mapfile = "crimap/chr9a.map")
sparrow.cmpmap9 <- parse_map_chrompic(chrompicfile = ("crimap/chr9a.cmp"))

#===
#===

##Same Parents flagged, remove Parents
sparrow.famped  <-  sparrow.famped[!(sparrow.famped$FATHER == "617"),]
sparrow.famped  <-  sparrow.famped[!(sparrow.famped$MOTHER == "549"),]
##Try again##
##Choose output directory where crimap function will create a folder crimap/ and place all crimap files
setwd("C:/Users/s1945757/PhD_Repo/Crimap Outputs Post Axiom Cleaning")
#Famped fixing again required
sparrow.famped <- sparrow.famped[which(sparrow.famped$ANIMAL %in% cleanspar@phdata$id),]
sparrow.famped <- sparrow.famped[sparrow.famped$Family %in% names(which(table(sparrow.famped$Family)>=5)),]
sparrow.famped$Family <- as.character(sparrow.famped$Family)
sparrow.famped$Family <- as.factor(sparrow.famped$Family)
#Run Crimap
Crymap(zzz= 9)
sparrow.map9 <- parse_map(mapfile = "crimap/chr9a.map")
sparrow.cmpmap9 <- parse_map_chrompic(chrompicfile = ("crimap/chr9a.cmp"))


##New families & Parents constantly being flagged... Also ind 1 family A/A, A/B, B/A, B/B lots of warnings, not sure why...
#Discuss with Susie
