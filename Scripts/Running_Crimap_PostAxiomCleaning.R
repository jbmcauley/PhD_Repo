
##Run a full map with new GenABEL Object##
library(GenABEL)
library(crimaptools)
library(dplyr)


#Load in new  GenABEL object
setwd("C:/Users/s1945757/PhD_Repo/Pdom200K_Sparrow_NTNU_96/clean_genotype_data")
load("Clean_Map_SNPs.RData")
   #Alternative if .RData file not loading properly:
        #setwd("C:/Users/johnb/Documents/PhD/Pdom200K_Sparrow_NTNU_96/clean_genotype_data")
        #cleanspar <- load.gwaa.data(phe = "Clean_Map_SNPs.phe", gen = "Clean_Map_SNPs.gen")


#Load in Famped file
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/crimap/crimap/crimap")
    #setwd("C:/Users/johnb/Documents/PhD/PhD_Repo/data")
sparrow.famped <- read.table("FamPed_20200414.txt", header = TRUE)

#Famped has ALL ind included, crimap requires the genABEL and family pedigree have all the same ind. So we remove the extra from the Famped object.
sparrow.famped <- sparrow.famped[which(sparrow.famped$ANIMAL %in% cleanspar@phdata$id),]
sparrow.famped <- sparrow.famped[sparrow.famped$Family %in% names(which(table(sparrow.famped$Family)>=5)),]
sparrow.famped$Family <- as.character(sparrow.famped$Family)
sparrow.famped$Family <- as.factor(sparrow.famped$Family)



##Perform basic GenABEL QC

qc0snp <- summary.snp.data(cleanspar@gtdata)
hist(qc0snp$CallRate, breaks = 50)

qc0id <- perid.summary(cleanspar)
hist(qc0id$CallPP, breaks = 50)
#Choose decent thresholds based on above histograms

qc1 <- check.marker(cleanspar, callrate = 0.9, perid.call = 0.8, p.level = 0)

data1 <- cleanspar[qc1$idok, qc1$snpok]
sparrow.abel <- data1


#Once again cleanup sparrow.famped
sparrow.famped <- sparrow.famped[which(sparrow.famped$ANIMAL %in% sparrow.abel@phdata$id),]
sparrow.famped <- sparrow.famped[sparrow.famped$Family %in% names(which(table(sparrow.famped$Family)>=5)),]
sparrow.famped$Family <- as.character(sparrow.famped$Family)
sparrow.famped$Family <- as.factor(sparrow.famped$Family)




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
#Evaluate

#===
#===
#Run Remaining Chr
Crymap(zzz= 10)
Crymap(zzz= 11)
Crymap(zzz= 12)
Crymap(zzz= 13)
Crymap(zzz= 14)
Crymap(zzz= 15)
Crymap(zzz= 17)
Crymap(zzz= 18)
Crymap(zzz= 19)
Crymap(zzz= 20)
Crymap(zzz= 21)
Crymap(zzz= 22)
Crymap(zzz= 23)
Crymap(zzz= 24)
Crymap(zzz= 25)
Crymap(zzz= 26)
Crymap(zzz= 27)
Crymap(zzz= 28)
Crymap(zzz= 30)


#Large Chrs: 1, 2, 3, 4, 5, 6, 7, 8, 29
chr1_length <- length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)])
chr2_length <- length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)])
chr3_length <- length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)])
chr4_length <- length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)])
chr5_length <- length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 5)])
chr6_length <- length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 6)])
chr7_length <- length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 7)])
chr8_length <- length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 8)])
chr29_length <- length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)])

#Run Chrs which are no longer too many snps for single crimap run
Crymap(zzz = 6)
Crymap(zzz = 7)
Crymap(zzz = 8)
Crymap(zzz = 5)







#Chrs longer than 5k: 1,2,3,4,29

length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)]) # 23366
a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][7800:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)])]





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr1a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1a.pre", genfile = "crimap/chr1a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr1a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1a.pre", genfile = "crimap/chr1a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr1a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr1a.gen", crimap.path =  "C:/PathApps/crimap.exe")








create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr1b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1b.pre", genfile = "crimap/chr1b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr1b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1b.pre", genfile = "crimap/chr1b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr1b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr1b.gen", crimap.path =  "C:/PathApps/crimap.exe")




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr1c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1c.pre", genfile = "crimap/chr1c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr1c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1c.pre", genfile = "crimap/chr1c.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr1c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr1c.gen", crimap.path =  "C:/PathApps/crimap.exe")







#chromosomes 1, 2 (32016), 3 (23580), 4 (15319), 5 (12216), 6 (6537), 7 (7122), 8 (8191), 29 (14950)


length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)]) # 23366

a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][7800:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)])]

create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr2a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2a.pre", genfile = "crimap/chr2a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr2a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2a.pre", genfile = "crimap/chr2a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr2a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr2b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2b.pre", genfile = "crimap/chr2b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr2b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2b.pre", genfile = "crimap/chr2b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr2b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2b.gen", crimap.path =  "C:/PathApps/crimap.exe")




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr2c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2c.pre", genfile = "crimap/chr2c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr2c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2c.pre", genfile = "crimap/chr2c.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr2c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2c.gen", crimap.path =  "C:/PathApps/crimap.exe")



##
##
##
##


#chromosomes 1, 2 (32016), 3 (23580), 4 (15319), 5 (12216), 6 (6537), 7 (7122), 8 (8191), 29 (14950)

length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)]) 


a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][7800:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)])]




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr3a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3a.pre", genfile = "crimap/chr3a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr3a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3a.pre", genfile = "crimap/chr3a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr3a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr3a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr3b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3b.pre", genfile = "crimap/chr3b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr3b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3b.pre", genfile = "crimap/chr3b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr3b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr3b.gen", crimap.path =  "C:/PathApps/crimap.exe")




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr3c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3c.pre", genfile = "crimap/chr3c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr3c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3c.pre", genfile = "crimap/chr3c.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr3c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr3c.gen", crimap.path =  "C:/PathApps/crimap.exe")


##
##
##
##

length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)]) # 5218


a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)][3900:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)])]





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr4a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4a.pre", genfile = "crimap/chr4a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr4a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4a.pre", genfile = "crimap/chr4a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr4a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr4a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr4b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4b.pre", genfile = "crimap/chr4b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr4b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4b.pre", genfile = "crimap/chr4b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr4b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr4b.gen", crimap.path =  "C:/PathApps/crimap.exe")






a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)][3900:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)])]


create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr29a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29a.pre", genfile = "crimap/chr29a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr29a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29a.pre", genfile = "crimap/chr29a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr29a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr29a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr29b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29b.pre", genfile = "crimap/chr29b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr29b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29b.pre", genfile = "crimap/chr29b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr29b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr29b.gen", crimap.path =  "C:/PathApps/crimap.exe")






###
###Graphs
###
###
###

setwd("C:/Users/s1945757/PhD_Repo/Crimap Outputs Post Axiom Cleaning")
#setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Data/Crimap Runs/Final")

#Set up file with BP
map <- as.data.frame(sparrow.abel@gtdata@map)
names(map) <- "BP"
map$SNP.Name <- row.names(map)
row.names(map) <- NULL

sparrow.map1a <- parse_map(mapfile = "crimap/chr1a.map")
sparrow.map1b <- parse_map(mapfile = "crimap/chr1b.map")
sparrow.map1c <- parse_map(mapfile = "crimap/chr1c.map")
#sparrow.map1d <- parse_map(mapfile = "crimap/chr1d.map")
#Combining Chr chunks

#adding cM pos within chunks
sparrow.map1b$cMPosition <- sparrow.map1b$cMPosition + (sparrow.map1a[which(sparrow.map1b$SNP.Name[1] == sparrow.map1a$SNP.Name),]$cMPosition)
sparrow.map1c$cMPosition <- sparrow.map1c$cMPosition + (sparrow.map1b[which(sparrow.map1c$SNP.Name[1] == sparrow.map1b$SNP.Name),]$cMPosition)
#sparrow.map1d$cMPosition <- sparrow.map1d$cMPosition + (sparrow.map1c[which(sparrow.map1d$SNP.Name[1] == sparrow.map1c$SNP.Name),]$cMPosition)
#adding M&F cM pos within Chunks
sparrow.map1b$cMPosition.Female <- sparrow.map1b$cMPosition.Female + (sparrow.map1a[which(sparrow.map1b$SNP.Name[1] == sparrow.map1a$SNP.Name),]$cMPosition.Female)
sparrow.map1c$cMPosition.Female <- sparrow.map1c$cMPosition.Female + (sparrow.map1b[which(sparrow.map1c$SNP.Name[1] == sparrow.map1b$SNP.Name),]$cMPosition.Female)
sparrow.map1d$cMPosition.Female <- sparrow.map1d$cMPosition.Female + (sparrow.map1c[which(sparrow.map1d$SNP.Name[1] == sparrow.map1c$SNP.Name),]$cMPosition.Female)
sparrow.map1b$cMPosition.Male <- sparrow.map1b$cMPosition.Male + (sparrow.map1a[which(sparrow.map1b$SNP.Name[1] == sparrow.map1a$SNP.Name),]$cMPosition.Male)
sparrow.map1c$cMPosition.Male <- sparrow.map1c$cMPosition.Male + (sparrow.map1b[which(sparrow.map1c$SNP.Name[1] == sparrow.map1b$SNP.Name),]$cMPosition.Male)
#sparrow.map1d$cMPosition.Male <- sparrow.map1d$cMPosition.Male + (sparrow.map1c[which(sparrow.map1d$SNP.Name[1] == sparrow.map1c$SNP.Name),]$cMPosition.Male)
#Adding order within chunks
sparrow.map1b$Order <- sparrow.map1b$Order + (sparrow.map1a[which(sparrow.map1b$SNP.Name[1] == sparrow.map1a$SNP.Name),]$Order)
sparrow.map1c$Order <- sparrow.map1c$Order + (sparrow.map1b[which(sparrow.map1c$SNP.Name[1] == sparrow.map1b$SNP.Name),]$Order)
#sparrow.map1d$Order <- sparrow.map1d$Order + (sparrow.map1c[which(sparrow.map1d$SNP.Name[1] == sparrow.map1c$SNP.Name),]$Order)

#Remove overlap
sparrow.map1a <- sparrow.map1a[-which(sparrow.map1a$SNP.Name %in% sparrow.map1b$SNP.Name),]
sparrow.map1b <- sparrow.map1b[-which(sparrow.map1b$SNP.Name %in% sparrow.map1c$SNP.Name),]
#sparrow.map1c <- sparrow.map1c[-which(sparrow.map1c$SNP.Name %in% sparrow.map1d$SNP.Name),]

#Merge
sparrow.map1a <- merge(sparrow.map1a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map1b <- merge(sparrow.map1b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map1c <- merge(sparrow.map1c, map, by.x = "SNP.Name",by.y = "SNP.Name")
#sparrow.map1d <- merge(sparrow.map1d, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map1 <- rbind(sparrow.map1a,sparrow.map1b,sparrow.map1c)
sparrow.map1[,c(5,6,7,8,10,11)] <- NULL

rm(sparrow.map1a)
rm(sparrow.map1b)
rm(sparrow.map1c)
rm(sparrow.map1d)





#Chr2
sparrow.map2a <- parse_map(mapfile = "crimap/chr2a.map")
sparrow.map2b <- parse_map(mapfile = "crimap/chr2b.map")
sparrow.map2c <- parse_map(mapfile = "crimap/chr2c.map")
#sparrow.map2d <- parse_map(mapfile = "crimap/chr2d.map")
#sparrow.map2e <- parse_map(mapfile = "crimap/chr2e.map")
#sparrow.map2f <- parse_map(mapfile = "crimap/chr2f.map")

sparrow.map2b$cMPosition <- sparrow.map2b$cMPosition + (sparrow.map2a[which(sparrow.map2b$SNP.Name[1] == sparrow.map2a$SNP.Name),]$cMPosition)
sparrow.map2c$cMPosition <- sparrow.map2c$cMPosition + (sparrow.map2b[which(sparrow.map2c$SNP.Name[1] == sparrow.map2b$SNP.Name),]$cMPosition)
#sparrow.map2d$cMPosition <- sparrow.map2d$cMPosition + (sparrow.map2c[which(sparrow.map2d$SNP.Name[1] == sparrow.map2c$SNP.Name),]$cMPosition)
#sparrow.map2e$cMPosition <- sparrow.map2e$cMPosition + (sparrow.map2d[which(sparrow.map2e$SNP.Name[1] == sparrow.map2d$SNP.Name),]$cMPosition)
#sparrow.map2f$cMPosition <- sparrow.map2f$cMPosition + (sparrow.map2e[which(sparrow.map2f$SNP.Name[1] == sparrow.map2e$SNP.Name),]$cMPosition)

sparrow.map2b$cMPosition.Female <- sparrow.map2b$cMPosition.Female + (sparrow.map2a[which(sparrow.map2b$SNP.Name[1] == sparrow.map2a$SNP.Name),]$cMPosition.Female)
sparrow.map2c$cMPosition.Female <- sparrow.map2c$cMPosition.Female + (sparrow.map2b[which(sparrow.map2c$SNP.Name[1] == sparrow.map2b$SNP.Name),]$cMPosition.Female)
#sparrow.map2d$cMPosition.Female <- sparrow.map2d$cMPosition.Female + (sparrow.map2c[which(sparrow.map2d$SNP.Name[1] == sparrow.map2c$SNP.Name),]$cMPosition.Female)
#sparrow.map2e$cMPosition.Female <- sparrow.map2e$cMPosition.Female + (sparrow.map2d[which(sparrow.map2e$SNP.Name[1] == sparrow.map2d$SNP.Name),]$cMPosition.Female)
#sparrow.map2f$cMPosition.Female <- sparrow.map2f$cMPosition.Female + (sparrow.map2e[which(sparrow.map2f$SNP.Name[1] == sparrow.map2e$SNP.Name),]$cMPosition.Female)
sparrow.map2b$cMPosition.Male <- sparrow.map2b$cMPosition.Male + (sparrow.map2a[which(sparrow.map2b$SNP.Name[1] == sparrow.map2a$SNP.Name),]$cMPosition.Male)
sparrow.map2c$cMPosition.Male <- sparrow.map2c$cMPosition.Male + (sparrow.map2b[which(sparrow.map2c$SNP.Name[1] == sparrow.map2b$SNP.Name),]$cMPosition.Male)
#sparrow.map2d$cMPosition.Male <- sparrow.map2d$cMPosition.Male + (sparrow.map2c[which(sparrow.map2d$SNP.Name[1] == sparrow.map2c$SNP.Name),]$cMPosition.Male)
#sparrow.map2e$cMPosition.Male <- sparrow.map2e$cMPosition.Male + (sparrow.map2d[which(sparrow.map2e$SNP.Name[1] == sparrow.map2d$SNP.Name),]$cMPosition.Male)
#sparrow.map2f$cMPosition.Male <- sparrow.map2f$cMPosition.Male + (sparrow.map2e[which(sparrow.map2f$SNP.Name[1] == sparrow.map2e$SNP.Name),]$cMPosition.Male)

#Adding order within chunks
sparrow.map2b$Order <- sparrow.map2b$Order + (sparrow.map2a[which(sparrow.map2b$SNP.Name[1] == sparrow.map2a$SNP.Name),]$Order)
sparrow.map2c$Order <- sparrow.map2c$Order + (sparrow.map2b[which(sparrow.map2c$SNP.Name[1] == sparrow.map2b$SNP.Name),]$Order)
#sparrow.map2d$Order <- sparrow.map2d$Order + (sparrow.map2c[which(sparrow.map2d$SNP.Name[1] == sparrow.map2c$SNP.Name),]$Order)
#sparrow.map2e$Order <- sparrow.map2e$Order + (sparrow.map2d[which(sparrow.map2e$SNP.Name[1] == sparrow.map2d$SNP.Name),]$Order)
#sparrow.map2f$Order <- sparrow.map2f$Order + (sparrow.map2e[which(sparrow.map2f$SNP.Name[1] == sparrow.map2e$SNP.Name),]$Order)

#Remove overlap
sparrow.map2a <- sparrow.map2a[-which(sparrow.map2a$SNP.Name %in% sparrow.map2b$SNP.Name),]
sparrow.map2b <- sparrow.map2b[-which(sparrow.map2b$SNP.Name %in% sparrow.map2c$SNP.Name),]
#sparrow.map2c <- sparrow.map2c[-which(sparrow.map2c$SNP.Name %in% sparrow.map2d$SNP.Name),]
#sparrow.map2d <- sparrow.map2d[-which(sparrow.map2d$SNP.Name %in% sparrow.map2e$SNP.Name),]
#sparrow.map2e <- sparrow.map2e[-which(sparrow.map2e$SNP.Name %in% sparrow.map2f$SNP.Name),]

#Merge
sparrow.map2a <- merge(sparrow.map2a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2b <- merge(sparrow.map2b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2c <- merge(sparrow.map2c, map, by.x = "SNP.Name",by.y = "SNP.Name")
#sparrow.map2d <- merge(sparrow.map2d, map, by.x = "SNP.Name",by.y = "SNP.Name")
#sparrow.map2e <- merge(sparrow.map2e, map, by.x = "SNP.Name",by.y = "SNP.Name")
#sparrow.map2f <- merge(sparrow.map2f, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map2 <- rbind(sparrow.map2a,sparrow.map2b,sparrow.map2c)
sparrow.map2[,c(5,6,7,8,10,11)] <- NULL

rm(sparrow.map2a)
rm(sparrow.map2b)
rm(sparrow.map2c)
rm(sparrow.map2d)
rm(sparrow.map2e)
rm(sparrow.map2f)




#Chr3
sparrow.map3a <- parse_map(mapfile = "crimap/chr3a.map")
sparrow.map3b <- parse_map(mapfile = "crimap/chr3b.map")
sparrow.map3c <- parse_map(mapfile = "crimap/chr3c.map")
sparrow.map3d <- parse_map(mapfile = "crimap/chr3d.map")
#sparrow.map3e <- parse_map(mapfile = "crimap/chr3e.map")


sparrow.map3b$cMPosition <- sparrow.map3b$cMPosition + (sparrow.map3a[which(sparrow.map3b$SNP.Name[1] == sparrow.map3a$SNP.Name),]$cMPosition)
sparrow.map3c$cMPosition <- sparrow.map3c$cMPosition + (sparrow.map3b[which(sparrow.map3c$SNP.Name[1] == sparrow.map3b$SNP.Name),]$cMPosition)
#sparrow.map3d$cMPosition <- sparrow.map3d$cMPosition + (sparrow.map3c[which(sparrow.map3d$SNP.Name[1] == sparrow.map3c$SNP.Name),]$cMPosition)
#sparrow.map3e$cMPosition <- sparrow.map3e$cMPosition + (sparrow.map3d[which(sparrow.map3e$SNP.Name[1] == sparrow.map3d$SNP.Name),]$cMPosition)

sparrow.map3b$cMPosition.Female <- sparrow.map3b$cMPosition.Female + (sparrow.map3a[which(sparrow.map3b$SNP.Name[1] == sparrow.map3a$SNP.Name),]$cMPosition.Female)
sparrow.map3c$cMPosition.Female <- sparrow.map3c$cMPosition.Female + (sparrow.map3b[which(sparrow.map3c$SNP.Name[1] == sparrow.map3b$SNP.Name),]$cMPosition.Female)
#sparrow.map3d$cMPosition.Female <- sparrow.map3d$cMPosition.Female + (sparrow.map3c[which(sparrow.map3d$SNP.Name[1] == sparrow.map3c$SNP.Name),]$cMPosition.Female)
#sparrow.map3e$cMPosition.Female <- sparrow.map3e$cMPosition.Female + (sparrow.map3d[which(sparrow.map3e$SNP.Name[1] == sparrow.map3d$SNP.Name),]$cMPosition.Female)

sparrow.map3b$cMPosition.Male <- sparrow.map3b$cMPosition.Male + (sparrow.map3a[which(sparrow.map3b$SNP.Name[1] == sparrow.map3a$SNP.Name),]$cMPosition.Male)
sparrow.map3c$cMPosition.Male <- sparrow.map3c$cMPosition.Male + (sparrow.map3b[which(sparrow.map3c$SNP.Name[1] == sparrow.map3b$SNP.Name),]$cMPosition.Male)
#sparrow.map3d$cMPosition.Male <- sparrow.map3d$cMPosition.Male + (sparrow.map3c[which(sparrow.map3d$SNP.Name[1] == sparrow.map3c$SNP.Name),]$cMPosition.Male)
#sparrow.map3e$cMPosition.Male <- sparrow.map3e$cMPosition.Male + (sparrow.map3d[which(sparrow.map3e$SNP.Name[1] == sparrow.map3d$SNP.Name),]$cMPosition.Male)

#Adding order within chunks
sparrow.map3b$Order <- sparrow.map3b$Order + (sparrow.map3a[which(sparrow.map3b$SNP.Name[1] == sparrow.map3a$SNP.Name),]$Order)
sparrow.map3c$Order <- sparrow.map3c$Order + (sparrow.map3b[which(sparrow.map3c$SNP.Name[1] == sparrow.map3b$SNP.Name),]$Order)
#sparrow.map3d$Order <- sparrow.map3d$Order + (sparrow.map3c[which(sparrow.map3d$SNP.Name[1] == sparrow.map3c$SNP.Name),]$Order)
#sparrow.map3e$Order <- sparrow.map3e$Order + (sparrow.map3d[which(sparrow.map3e$SNP.Name[1] == sparrow.map3d$SNP.Name),]$Order)

#Remove overlap
sparrow.map3a <- sparrow.map3a[-which(sparrow.map3a$SNP.Name %in% sparrow.map3b$SNP.Name),]
sparrow.map3b <- sparrow.map3b[-which(sparrow.map3b$SNP.Name %in% sparrow.map3c$SNP.Name),]
#sparrow.map3c <- sparrow.map3c[-which(sparrow.map3c$SNP.Name %in% sparrow.map3d$SNP.Name),]
#sparrow.map3d <- sparrow.map3d[-which(sparrow.map3d$SNP.Name %in% sparrow.map3e$SNP.Name),]

#Merge
sparrow.map3a <- merge(sparrow.map3a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map3b <- merge(sparrow.map3b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map3c <- merge(sparrow.map3c, map, by.x = "SNP.Name",by.y = "SNP.Name")
#sparrow.map3d <- merge(sparrow.map3d, map, by.x = "SNP.Name",by.y = "SNP.Name")
#sparrow.map3e <- merge(sparrow.map3e, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map3 <- rbind(sparrow.map3a,sparrow.map3b,sparrow.map3c)
sparrow.map3[,c(5,6,7,8,10,11)] <- NULL

rm(sparrow.map3a)
rm(sparrow.map3b)
rm(sparrow.map3c)
rm(sparrow.map3d)
rm(sparrow.map3e)





#Chr 4
sparrow.map4a <- parse_map(mapfile = "crimap/chr4a.map")
sparrow.map4b <- parse_map(mapfile = "crimap/chr4b.map")
sparrow.map4c <- parse_map(mapfile = "crimap/chr4c.map")

#Combining Chr chunks

#adding cM pos within chunks
sparrow.map4b$cMPosition <- sparrow.map4b$cMPosition + (sparrow.map4a[which(sparrow.map4b$SNP.Name[1] == sparrow.map4a$SNP.Name),]$cMPosition)
#sparrow.map4c$cMPosition <- sparrow.map4c$cMPosition + (sparrow.map4b[which(sparrow.map4c$SNP.Name[1] == sparrow.map4b$SNP.Name),]$cMPosition)

#adding M&F cM pos within Chunks
sparrow.map4b$cMPosition.Female <- sparrow.map4b$cMPosition.Female + (sparrow.map4a[which(sparrow.map4b$SNP.Name[1] == sparrow.map4a$SNP.Name),]$cMPosition.Female)
#sparrow.map4c$cMPosition.Female <- sparrow.map4c$cMPosition.Female + (sparrow.map4b[which(sparrow.map4c$SNP.Name[1] == sparrow.map4b$SNP.Name),]$cMPosition.Female)

sparrow.map4b$cMPosition.Male <- sparrow.map4b$cMPosition.Male + (sparrow.map4a[which(sparrow.map4b$SNP.Name[1] == sparrow.map4a$SNP.Name),]$cMPosition.Male)
#sparrow.map4c$cMPosition.Male <- sparrow.map4c$cMPosition.Male + (sparrow.map4b[which(sparrow.map4c$SNP.Name[1] == sparrow.map4b$SNP.Name),]$cMPosition.Male)

#Adding order within chunks
sparrow.map4b$Order <- sparrow.map4b$Order + (sparrow.map4a[which(sparrow.map4b$SNP.Name[1] == sparrow.map4a$SNP.Name),]$Order)
#sparrow.map4c$Order <- sparrow.map4c$Order + (sparrow.map4b[which(sparrow.map4c$SNP.Name[1] == sparrow.map4b$SNP.Name),]$Order)

#Remove overlap
sparrow.map4a <- sparrow.map4a[-which(sparrow.map4a$SNP.Name %in% sparrow.map4b$SNP.Name),]
#sparrow.map4b <- sparrow.map4b[-which(sparrow.map4b$SNP.Name %in% sparrow.map4c$SNP.Name),]

#Merge
sparrow.map4a <- merge(sparrow.map4a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map4b <- merge(sparrow.map4b, map, by.x = "SNP.Name",by.y = "SNP.Name")
#sparrow.map4c <- merge(sparrow.map4c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map4 <- rbind(sparrow.map4a,sparrow.map4b)
sparrow.map4[,c(5,6,7,8,10,11)] <- NULL

rm(sparrow.map4a)
rm(sparrow.map4b)
rm(sparrow.map4c)




#Chr 29
sparrow.map29a <- parse_map(mapfile = "crimap/chr29a.map")
sparrow.map29b <- parse_map(mapfile = "crimap/chr29b.map")
#sparrow.map29c <- parse_map(mapfile = "crimap/chr29c.map")

#Combining Chr chunks
#adding cM pos within chunks
sparrow.map29b$cMPosition <- sparrow.map29b$cMPosition + (sparrow.map29a[which(sparrow.map29b$SNP.Name[1] == sparrow.map29a$SNP.Name),]$cMPosition)
#sparrow.map29c$cMPosition <- sparrow.map29c$cMPosition + (sparrow.map29b[which(sparrow.map29c$SNP.Name[1] == sparrow.map29b$SNP.Name),]$cMPosition)

#adding M&F cM pos within Chunks
sparrow.map29b$cMPosition.Female <- sparrow.map29b$cMPosition.Female + (sparrow.map29a[which(sparrow.map29b$SNP.Name[1] == sparrow.map29a$SNP.Name),]$cMPosition.Female)
#sparrow.map29c$cMPosition.Female <- sparrow.map29c$cMPosition.Female + (sparrow.map29b[which(sparrow.map29c$SNP.Name[1] == sparrow.map29b$SNP.Name),]$cMPosition.Female)

sparrow.map29b$cMPosition.Male <- sparrow.map29b$cMPosition.Male + (sparrow.map29a[which(sparrow.map29b$SNP.Name[1] == sparrow.map29a$SNP.Name),]$cMPosition.Male)
#sparrow.map29c$cMPosition.Male <- sparrow.map29c$cMPosition.Male + (sparrow.map29b[which(sparrow.map29c$SNP.Name[1] == sparrow.map29b$SNP.Name),]$cMPosition.Male)

#Adding order within chunks
sparrow.map29b$Order <- sparrow.map29b$Order + (sparrow.map29a[which(sparrow.map29b$SNP.Name[1] == sparrow.map29a$SNP.Name),]$Order)
#sparrow.map29c$Order <- sparrow.map29c$Order + (sparrow.map29b[which(sparrow.map29c$SNP.Name[1] == sparrow.map29b$SNP.Name),]$Order)

#Remove overlap
sparrow.map29a <- sparrow.map29a[-which(sparrow.map29a$SNP.Name %in% sparrow.map29b$SNP.Name),]
#sparrow.map29b <- sparrow.map29b[-which(sparrow.map29b$SNP.Name %in% sparrow.map29c$SNP.Name),]

#Merge
sparrow.map29a <- merge(sparrow.map29a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map29b <- merge(sparrow.map29b, map, by.x = "SNP.Name",by.y = "SNP.Name")
#sparrow.map29c <- merge(sparrow.map29c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map29 <- rbind(sparrow.map29a,sparrow.map29b)
sparrow.map29[,c(5,6,7,8,10,11)] <- NULL
rm(sparrow.map29a)
rm(sparrow.map29b)
rm(sparrow.map29c)


#Ch 5
sparrow.map5 <- parse_map(mapfile = "crimap/chr5a.map")
sparrow.map5 <- merge(sparrow.map5, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map5[,c(5,6,7,8,10,11)] <- NULL

#Ch 6
sparrow.map6 <- parse_map(mapfile = "crimap/chr6a.map")
sparrow.map6 <- merge(sparrow.map6, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map6[,c(5,6,7,8,10,11)] <- NULL

#Ch 7
sparrow.map7 <- parse_map(mapfile = "crimap/chr7a.map")
sparrow.map7 <- merge(sparrow.map7, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map7[,c(5,6,7,8,10,11)] <- NULL

#Ch 8
sparrow.map8 <- parse_map(mapfile = "crimap/chr8a.map")
sparrow.map8 <- merge(sparrow.map8, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map8[,c(5,6,7,8,10,11)] <- NULL

#Chr 9
sparrow.map9 <- parse_map(mapfile = "crimap/chr9a.map")
sparrow.map9 <- merge(sparrow.map9, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map9[,c(5,6,7,8,10,11)] <- NULL

#Chr 10
sparrow.map10 <- parse_map(mapfile = "crimap/chr10a.map")
sparrow.map10 <- merge(sparrow.map10, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map10[,c(5,6,7,8,10,11)] <- NULL

#Chr 11
sparrow.map11 <- parse_map(mapfile = "crimap/chr11a.map")
sparrow.map11 <- merge(sparrow.map11, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map11[,c(5,6,7,8,10,11)] <- NULL

#Chr 12
sparrow.map12 <- parse_map(mapfile = "crimap/chr12a.map")
sparrow.map12 <- merge(sparrow.map12, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map12[,c(5,6,7,8,10,11)] <- NULL

#Chr 13
sparrow.map13 <- parse_map(mapfile = "crimap/chr13a.map")
sparrow.map13 <- merge(sparrow.map13, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map13[,c(5,6,7,8,10,11)] <- NULL

#Chr 14
sparrow.map14 <- parse_map(mapfile = "crimap/chr14a.map")
sparrow.map14 <- merge(sparrow.map14, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map14[,c(5,6,7,8,10,11)] <- NULL

#Chr 15
sparrow.map15 <- parse_map(mapfile = "crimap/chr15a.map")
sparrow.map15 <- merge(sparrow.map15, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map15[,c(5,6,7,8,10,11)] <- NULL

#Chr 17
sparrow.map17 <- parse_map(mapfile = "crimap/chr17a.map")
sparrow.map17 <- merge(sparrow.map17, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map17[,c(5,6,7,8,10,11)] <- NULL

#Chr 18
sparrow.map18 <- parse_map(mapfile = "crimap/chr18a.map")
sparrow.map18 <- merge(sparrow.map18, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map18[,c(5,6,7,8,10,11)] <- NULL

#Chr 19
sparrow.map19 <- parse_map(mapfile = "crimap/chr19a.map")
sparrow.map19 <- merge(sparrow.map19, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map19[,c(5,6,7,8,10,11)] <- NULL

#Chr 20
sparrow.map20 <- parse_map(mapfile = "crimap/chr20a.map")
sparrow.map20 <- merge(sparrow.map20, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map20[,c(5,6,7,8,10,11)] <- NULL

#Chr 21
sparrow.map21 <- parse_map(mapfile = "crimap/chr21a.map")
sparrow.map21 <- merge(sparrow.map21, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map21[,c(5,6,7,8,10,11)] <- NULL

#Chr 22
sparrow.map22 <- parse_map(mapfile = "crimap/chr22a.map")
sparrow.map22 <- merge(sparrow.map22, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map22[,c(5,6,7,8,10,11)] <- NULL

#Chr 23
sparrow.map23 <- parse_map(mapfile = "crimap/chr23a.map")
sparrow.map23 <- merge(sparrow.map23, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map23[,c(5,6,7,8,10,11)] <- NULL

#Chr 24
sparrow.map24 <- parse_map(mapfile = "crimap/chr24a.map")
sparrow.map24 <- merge(sparrow.map24, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map24[,c(5,6,7,8,10,11)] <- NULL

#Chr 25
sparrow.map25 <- parse_map(mapfile = "crimap/chr25a.map")
sparrow.map25 <- merge(sparrow.map25, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map25[,c(5,6,7,8,10,11)] <- NULL

#Chr 26
sparrow.map26 <- parse_map(mapfile = "crimap/chr26a.map")
sparrow.map26 <- merge(sparrow.map26, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map26[,c(5,6,7,8,10,11)] <- NULL

#Chr 27
sparrow.map27 <- parse_map(mapfile = "crimap/chr27a.map")
sparrow.map27 <- merge(sparrow.map27, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map27[,c(5,6,7,8,10,11)] <- NULL

#Chr 28
sparrow.map28 <- parse_map(mapfile = "crimap/chr28a.map")
sparrow.map28 <- merge(sparrow.map28, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map28[,c(5,6,7,8,10,11)] <- NULL

sparrow.map1$Chr <- "Chr 1"
sparrow.map2$Chr <- "Chr 2"
sparrow.map3$Chr <- "Chr 3"
sparrow.map4$Chr <- "Chr 4"
sparrow.map5$Chr <- "Chr 5"
sparrow.map6$Chr <- "Chr 6"
sparrow.map7$Chr <- "Chr 7"
sparrow.map8$Chr <- "Chr 8"
sparrow.map9$Chr <- "Chr 9"
sparrow.map10$Chr <- "Chr 10"
sparrow.map11$Chr <- "Chr 11"
sparrow.map12$Chr <- "Chr 12"
sparrow.map13$Chr <- "Chr 13"
sparrow.map14$Chr <- "Chr 14"
sparrow.map15$Chr <- "Chr 15"
sparrow.map17$Chr <- "Chr 17"
sparrow.map18$Chr <- "Chr 18"
sparrow.map19$Chr <- "Chr 19"
sparrow.map20$Chr <- "Chr 20"
sparrow.map21$Chr <- "Chr 21"
sparrow.map22$Chr <- "Chr 22"
sparrow.map23$Chr <- "Chr 23" 
sparrow.map24$Chr <- "Chr 24"
sparrow.map25$Chr <- "Chr 25"
sparrow.map26$Chr <- "Chr 26"
sparrow.map27$Chr <- "Chr 27"
sparrow.map28$Chr <- "Chr 28"
sparrow.map29$Chr <- "Chr 29"
sparrow.map1 <- sparrow.map1[order(sparrow.map1$Order),]
sparrow.map2 <- sparrow.map2[order(sparrow.map2$Order),]
sparrow.map3 <- sparrow.map3[order(sparrow.map3$Order),]
sparrow.map4 <- sparrow.map4[order(sparrow.map4$Order),]
sparrow.map5 <- sparrow.map5[order(sparrow.map5$Order),]
sparrow.map6 <- sparrow.map6[order(sparrow.map6$Order),]
sparrow.map7 <- sparrow.map7[order(sparrow.map7$Order),]
sparrow.map8 <- sparrow.map8[order(sparrow.map8$Order),]
sparrow.map9 <- sparrow.map9[order(sparrow.map9$Order),]
sparrow.map10 <- sparrow.map10[order(sparrow.map10$Order),]
sparrow.map11 <- sparrow.map11[order(sparrow.map11$Order),]
sparrow.map12 <- sparrow.map12[order(sparrow.map12$Order),]
sparrow.map13 <- sparrow.map13[order(sparrow.map13$Order),]
sparrow.map14 <- sparrow.map14[order(sparrow.map14$Order),]
sparrow.map15 <- sparrow.map15[order(sparrow.map15$Order),]
sparrow.map17 <- sparrow.map17[order(sparrow.map17$Order),]
sparrow.map18 <- sparrow.map18[order(sparrow.map18$Order),]
sparrow.map19 <- sparrow.map19[order(sparrow.map19$Order),]
sparrow.map20 <- sparrow.map20[order(sparrow.map20$Order),]
sparrow.map21 <- sparrow.map21[order(sparrow.map21$Order),]
sparrow.map22 <- sparrow.map22[order(sparrow.map22$Order),]
sparrow.map23 <- sparrow.map23[order(sparrow.map23$Order),]
sparrow.map24 <- sparrow.map24[order(sparrow.map24$Order),]
sparrow.map25 <- sparrow.map25[order(sparrow.map25$Order),]
sparrow.map26 <- sparrow.map26[order(sparrow.map26$Order),]
sparrow.map27 <- sparrow.map27[order(sparrow.map27$Order),]
sparrow.map28 <- sparrow.map28[order(sparrow.map28$Order),]
sparrow.map29 <- sparrow.map29[order(sparrow.map29$Order),]



sparrow.map.all <- rbind(sparrow.map1,sparrow.map2,sparrow.map3,sparrow.map4,sparrow.map5,sparrow.map6,sparrow.map7,sparrow.map8,
                         sparrow.map9,sparrow.map10,sparrow.map11,sparrow.map12,sparrow.map13,sparrow.map14,sparrow.map15,sparrow.map17,
                         sparrow.map18,sparrow.map19,sparrow.map20,sparrow.map21,sparrow.map22,sparrow.map23,sparrow.map24,sparrow.map25,sparrow.map26,
                         sparrow.map27,sparrow.map28,sparrow.map29)
sparrow.map.all$Chr <- factor(sparrow.map.all$Chr, levels=c("Chr 1","Chr 2","Chr 3","Chr 4","Chr 5","Chr 6","Chr 7","Chr 8","Chr 9",
                                                            "Chr 10","Chr 11","Chr 12","Chr 13","Chr 14","Chr 15","Chr 17","Chr 18",
                                                            "Chr 19","Chr 20","Chr 21","Chr 22","Chr 23","Chr 24","Chr 25","Chr 26","Chr 27",
                                                            "Chr 28","Chr 29"))
library(ggplot2)
p <-ggplot(data = sparrow.map.all, aes(x = BP/1000000, y = cMPosition)) +
  geom_point(size = 1) +
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")
p + facet_wrap(~Chr, scales = "free")

p <-ggplot(data = sparrow.map.all, aes(x = BP/1000000, y = cMPosition.Female)) +
  geom_point(size = 1) +
  labs(y= "Female Linkage Map Length (cM)", x ="Female Chromosome Length (Mb)")
p + facet_wrap(~Chr, scales = "free")

p <-ggplot(data = sparrow.map.all, aes(x = BP/1000000, y = cMPosition.Male)) +
  geom_point(size = 1) +
  labs(y= "Male Linkage Map Length (cM)", x ="Male Chromosome Length (Mb)")
p + facet_wrap(~Chr, scales = "free")

p <-ggplot(data = sparrow.map.all) +
  geom_point(aes(x = BP/1000000, y = cMPosition.Male), color = "blue", size = 1) +
  geom_point(aes(x = BP/1000000, y = cMPosition.Female), color = "red", size = 1) +
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")
p + facet_wrap(~Chr, scales = "free")

#setwd("C:/Users/johnb/Documents/PhD/PhD_Repo/data")
save(sparrow.map.all, file = "sparrow.map.all.RData")
load("sparrow.map.all.RData")

###
###
###
###
###
###Script for removing odd/suspect snps
###
###
###
###


# Goal is to create a loop of some kind which identifies snps worth removing and assigns them to a object which can 
# then be used to remove those snps from the main GenABEL object

# Final remove function to look something like this:
# sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]

sparrow.map.all_list <- split(sparrow.map.all, sparrow.map.all$Chr)



# Generate logical vector for each Chr on SNPs which have cM less than 1 on both adjacent SNPs


chrom_vec <- levels(sparrow.map.all$Chr)
snpremove <- list()
final_snpremove <- list()
counter <- 0
for(j in 1:length(chrom_vec)){

z <- data.frame(sparrow.map.all_list[j])
names(z) <- names(sparrow.map.all)

snpremove <- as.list(NULL)
counter <- 0
for(i in 2:(length(z$SNP.Name)-1)){
  counter <- counter + 1
  snpremove[[counter]] <- (abs(z$cMPosition[i] - z$cMPosition[i-1]) < 1)  & (abs(z$cMPosition[i] - z$cMPosition[i+1]) < 1)
  
}

snpremove <- do.call(rbind, snpremove)
snpremove <- as.vector(snpremove) 
snpremove <- c("TRUE", snpremove)
snpremove <- c(snpremove, "TRUE")
snpremove <- as.logical(snpremove)


final_snpremove[[j]] <- snpremove

}
names(final_snpremove) <- c(chrom_vec)

save(final_snpremove, file = "SNP_Removal_Logical_Vectors.RData")

#Obvious Chrs to apply filtering to: 1,2,3,4,5,6,7,8,10,11,12,14,15

load("SNP_Removal_Logical_Vectors.RData")

final_snpremove$`Chr 9`[which(final_snpremove$`Chr 9` == FALSE)] <- TRUE
final_snpremove$`Chr 13`[which(final_snpremove$`Chr 13` == FALSE)] <- TRUE
final_snpremove$`Chr 17`[which(final_snpremove$`Chr 17` == FALSE)] <- TRUE
final_snpremove$`Chr 18`[which(final_snpremove$`Chr 18` == FALSE)] <- TRUE
final_snpremove$`Chr 19`[which(final_snpremove$`Chr 19` == FALSE)] <- TRUE
final_snpremove$`Chr 20`[which(final_snpremove$`Chr 20` == FALSE)] <- TRUE
final_snpremove$`Chr 21`[which(final_snpremove$`Chr 21` == FALSE)] <- TRUE
final_snpremove$`Chr 22`[which(final_snpremove$`Chr 22` == FALSE)] <- TRUE
final_snpremove$`Chr 23`[which(final_snpremove$`Chr 23` == FALSE)] <- TRUE
final_snpremove$`Chr 24`[which(final_snpremove$`Chr 24` == FALSE)] <- TRUE
final_snpremove$`Chr 26`[which(final_snpremove$`Chr 26` == FALSE)] <- TRUE
final_snpremove$`Chr 27`[which(final_snpremove$`Chr 27` == FALSE)] <- TRUE
final_snpremove$`Chr 28`[which(final_snpremove$`Chr 28` == FALSE)] <- TRUE

test <- c(final_snpremove$`Chr 1`,final_snpremove$`Chr 2`,final_snpremove$`Chr 3`,final_snpremove$`Chr 4`,final_snpremove$`Chr 5`,
          final_snpremove$`Chr 6`,final_snpremove$`Chr 7`,final_snpremove$`Chr 8`,final_snpremove$`Chr 9`,final_snpremove$`Chr 10`,
          final_snpremove$`Chr 11`,final_snpremove$`Chr 12`,final_snpremove$`Chr 13`,final_snpremove$`Chr 14`,final_snpremove$`Chr 15`,
          final_snpremove$`Chr 17`,final_snpremove$`Chr 18`,final_snpremove$`Chr 19`,final_snpremove$`Chr 20`,final_snpremove$`Chr 21`,
          final_snpremove$`Chr 22`,final_snpremove$`Chr 23`,final_snpremove$`Chr 24`,final_snpremove$`Chr 26`,final_snpremove$`Chr 27`,
          final_snpremove$`Chr 28`,final_snpremove$`Chr 29`)

sparrow.abel <- sparrow.abel[,snpnames(sparrow.abel)[test]]



#Above filtering didnt seem to work.... Try again:

chr1logi <- final_snpremove$`Chr 1`
chr2logi <- final_snpremove$`Chr 2`
chr3logi <- final_snpremove$`Chr 3`
chr4logi <- final_snpremove$`Chr 4`
chr5logi <- final_snpremove$`Chr 5`
chr6logi <- final_snpremove$`Chr 6`
chr7logi <- final_snpremove$`Chr 7`
chr8logi <- final_snpremove$`Chr 8`
chr9logi <- final_snpremove$`Chr 9`
chr10logi <- final_snpremove$`Chr 10`
chr11logi <- final_snpremove$`Chr 11`
chr12logi <- final_snpremove$`Chr 12`
chr13logi <- final_snpremove$`Chr 13`
chr14logi <- final_snpremove$`Chr 14`
chr15logi <- final_snpremove$`Chr 15`
#chr16logi <- final_snpremove$`Chr 16`
chr17logi <- final_snpremove$`Chr 17`
chr18logi <- final_snpremove$`Chr 18`
chr19logi <- final_snpremove$`Chr 19`
chr20logi <- final_snpremove$`Chr 20`
chr21logi <- final_snpremove$`Chr 21`
chr22logi <- final_snpremove$`Chr 22`
chr23logi <- final_snpremove$`Chr 23`
chr24logi <- final_snpremove$`Chr 24`
#chr25logi <- final_snpremove$`Chr 25`
chr26logi <- final_snpremove$`Chr 26`
chr27logi <- final_snpremove$`Chr 27`
chr28logi <- final_snpremove$`Chr 28`
chr29logi <- final_snpremove$`Chr 29`

chr1num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "1")])
chr2num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "2")])
chr3num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "3")])
chr4num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "4")])
chr5num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "5")])
chr6num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "6")])
chr7num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "7")])
chr8num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "8")])
chr9num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "9")])
chr10num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "10")])
chr11num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "11")])
chr12num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "12")])
chr13num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "13")])
chr14num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "14")])
chr15num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "15")])
#chr16num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "16")])
chr17num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "17")])
chr18num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "18")])
chr19num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "19")])
chr20num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "20")])
chr21num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "21")])
chr22num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "22")])
chr23num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "23")])
chr24num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "24")])
#chr25num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "25")])
chr26num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "26")])
chr27num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "27")])
chr28num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "28")])
chr29num <- which(sparrow.abel@gtdata@chromosome %in% levels(sparrow.abel@gtdata@chromosome)[which(levels(sparrow.abel@gtdata@chromosome) == "29")])


#1,2,3,4,5,6,7,8,10,11,12,14,15
snpremove <- sparrow.abel@gtdata@snpnames[chr1num[!chr1num %in% chr1num[chr1logi]]]
snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr2num[!chr2num %in% chr2num[chr2logi]]])
snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr3num[!chr3num %in% chr3num[chr3logi]]])
snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr4num[!chr4num %in% chr4num[chr4logi]]])
snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr5num[!chr5num %in% chr5num[chr5logi]]])
snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr6num[!chr6num %in% chr6num[chr6logi]]])
snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr7num[!chr7num %in% chr7num[chr7logi]]])
snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr8num[!chr8num %in% chr8num[chr8logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr9num[!chr9num %in% chr9num[chr9logi]]])
snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr10num[!chr10num %in% chr10num[chr10logi]]])
snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr11num[!chr11num %in% chr11num[chr11logi]]])
snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr12num[!chr12num %in% chr12num[chr12logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr13num[!chr13num %in% chr13num[chr13logi]]])
snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr14num[!chr14num %in% chr14num[chr14logi]]])
snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr15num[!chr15num %in% chr15num[chr15logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr16num[!chr16num %in% chr16num[chr16logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr17num[!chr17num %in% chr17num[chr17logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr18num[!chr18num %in% chr18num[chr18logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr19num[!chr19num %in% chr19num[chr19logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr20num[!chr20num %in% chr20num[chr20logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr21num[!chr21num %in% chr21num[chr21logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr22num[!chr22num %in% chr22num[chr22logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr23num[!chr23num %in% chr23num[chr23logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr24num[!chr24num %in% chr24num[chr24logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr25num[!chr25num %in% chr25num[chr25logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr26num[!chr26num %in% chr26num[chr26logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr27num[!chr27num %in% chr27num[chr27logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr28num[!chr28num %in% chr28num[chr28logi]]])
#snpremove <- c(snpremove, sparrow.abel@gtdata@snpnames[chr29num[!chr29num %in% chr29num[chr29logi]]])



sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snpremove]


#rerun crimap on chrs: 1,2,3,4,5,6,7,8,10,11,12,14,15


Crymap(5)
Crymap(6)
Crymap(7)
Crymap(8)
Crymap(10)
Crymap(11)
Crymap(12)
Crymap(14)
Crymap(15)

# LargeChrs====
length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)]) # 23366
a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][7800:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)])]


create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr1a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1a.pre", genfile = "crimap/chr1a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr1a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1a.pre", genfile = "crimap/chr1a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr1a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr1a.gen", crimap.path =  "C:/PathApps/crimap.exe")








create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr1b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1b.pre", genfile = "crimap/chr1b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr1b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1b.pre", genfile = "crimap/chr1b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr1b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr1b.gen", crimap.path =  "C:/PathApps/crimap.exe")




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr1c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1c.pre", genfile = "crimap/chr1c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr1c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1c.pre", genfile = "crimap/chr1c.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr1c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr1c.gen", crimap.path =  "C:/PathApps/crimap.exe")







#chromosomes 1, 2 (32016), 3 (23580), 4 (15319), 5 (12216), 6 (6537), 7 (7122), 8 (8191), 29 (14950)


length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)]) # 23366

a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][7800:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)])]

create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr2a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2a.pre", genfile = "crimap/chr2a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr2a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2a.pre", genfile = "crimap/chr2a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr2a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr2b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2b.pre", genfile = "crimap/chr2b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr2b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2b.pre", genfile = "crimap/chr2b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr2b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2b.gen", crimap.path =  "C:/PathApps/crimap.exe")


create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr2c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2c.pre", genfile = "crimap/chr2c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr2c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2c.pre", genfile = "crimap/chr2c.gen", familyPedigree = sparrow.famped)
run_crimap_map(genfile = "crimap/chr2c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2c.gen", crimap.path =  "C:/PathApps/crimap.exe")



##
##
##
##


#chromosomes 1, 2 (32016), 3 (23580), 4 (15319), 5 (12216), 6 (6537), 7 (7122), 8 (8191), 29 (14950)

length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)]) 


a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][7800:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)])]


create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr3a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3a.pre", genfile = "crimap/chr3a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr3a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3a.pre", genfile = "crimap/chr3a.gen", familyPedigree = sparrow.famped)
run_crimap_map(genfile = "crimap/chr3a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr3a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr3b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3b.pre", genfile = "crimap/chr3b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr3b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3b.pre", genfile = "crimap/chr3b.gen", familyPedigree = sparrow.famped)
run_crimap_map(genfile = "crimap/chr3b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr3b.gen", crimap.path =  "C:/PathApps/crimap.exe")


create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr3c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3c.pre", genfile = "crimap/chr3c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr3c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3c.pre", genfile = "crimap/chr3c.gen", familyPedigree = sparrow.famped)
run_crimap_map(genfile = "crimap/chr3c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr3c.gen", crimap.path =  "C:/PathApps/crimap.exe")


##
##
##
##

length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)]) # 5218


a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)][3900:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)])]


create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr4a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4a.pre", genfile = "crimap/chr4a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr4a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4a.pre", genfile = "crimap/chr4a.gen", familyPedigree = sparrow.famped)
run_crimap_map(genfile = "crimap/chr4a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr4a.gen", crimap.path =  "C:/PathApps/crimap.exe")


create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr4b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4b.pre", genfile = "crimap/chr4b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr4b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4b.pre", genfile = "crimap/chr4b.gen", familyPedigree = sparrow.famped)
run_crimap_map(genfile = "crimap/chr4b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr4b.gen", crimap.path =  "C:/PathApps/crimap.exe")






a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)][3900:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)])]


create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr29a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29a.pre", genfile = "crimap/chr29a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr29a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29a.pre", genfile = "crimap/chr29a.gen", familyPedigree = sparrow.famped)
run_crimap_map(genfile = "crimap/chr29a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr29a.gen", crimap.path =  "C:/PathApps/crimap.exe")



create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr29b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29b.pre", genfile = "crimap/chr29b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr29b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29b.pre", genfile = "crimap/chr29b.gen", familyPedigree = sparrow.famped)
run_crimap_map(genfile = "crimap/chr29b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr29b.gen", crimap.path =  "C:/PathApps/crimap.exe")



###







#---- 
sparrowgen <- sparrow.abel
setwd("~/PhD/PhD_Repo/data/Crimap Runs/Final")
sparrow.xover1a <- parse_crossovers(chrompicfile = "crimap/chr1a.cmp", familyPedigree = sparrow.famped)
sparrow.xover1b <- parse_crossovers(chrompicfile = "crimap/chr1b.cmp", familyPedigree = sparrow.famped)
sparrow.xover1c <- parse_crossovers(chrompicfile = "crimap/chr1c.cmp", familyPedigree = sparrow.famped)
#sparrow.xover1d <- parse_crossovers(chrompicfile = "crimap/chr1d.cmp", familyPedigree = sparrow.famped)
#sparrow.xover1e <- parse_crossovers(chrompicfile = "crimap/chr1e.cmp", familyPedigree = sparrow.famped)
sparrow.xover2a <- parse_crossovers(chrompicfile = "crimap/chr2a.cmp", familyPedigree = sparrow.famped)
sparrow.xover2b <- parse_crossovers(chrompicfile = "crimap/chr2b.cmp", familyPedigree = sparrow.famped)
sparrow.xover2c <- parse_crossovers(chrompicfile = "crimap/chr2c.cmp", familyPedigree = sparrow.famped)
#sparrow.xover2d <- parse_crossovers(chrompicfile = "crimap/chr2d.cmp", familyPedigree = sparrow.famped)
#sparrow.xover2e <- parse_crossovers(chrompicfile = "crimap/chr2e.cmp", familyPedigree = sparrow.famped)
#sparrow.xover2f <- parse_crossovers(chrompicfile = "crimap/chr2f.cmp", familyPedigree = sparrow.famped)
#sparrow.xover2g <- parse_crossovers(chrompicfile = "crimap/chr2g.cmp", familyPedigree = sparrow.famped)
sparrow.xover3a <- parse_crossovers(chrompicfile = "crimap/chr3a.cmp", familyPedigree = sparrow.famped)
sparrow.xover3b <- parse_crossovers(chrompicfile = "crimap/chr3b.cmp", familyPedigree = sparrow.famped)
sparrow.xover3c <- parse_crossovers(chrompicfile = "crimap/chr3c.cmp", familyPedigree = sparrow.famped)
#sparrow.xover3d <- parse_crossovers(chrompicfile = "crimap/chr3d.cmp", familyPedigree = sparrow.famped)
#sparrow.xover3e <- parse_crossovers(chrompicfile = "crimap/chr3e.cmp", familyPedigree = sparrow.famped)
sparrow.xover4a <- parse_crossovers(chrompicfile = "crimap/chr4a.cmp", familyPedigree = sparrow.famped)
sparrow.xover4b <- parse_crossovers(chrompicfile = "crimap/chr4b.cmp", familyPedigree = sparrow.famped)
#sparrow.xover4c <- parse_crossovers(chrompicfile = "crimap/chr4c.cmp", familyPedigree = sparrow.famped)
sparrow.xover5a <- parse_crossovers(chrompicfile = "crimap/chr5a.cmp", familyPedigree = sparrow.famped)
#sparrow.xover5b <- parse_crossovers(chrompicfile = "crimap/chr5b.cmp", familyPedigree = sparrow.famped)
#sparrow.xover5c <- parse_crossovers(chrompicfile = "crimap/chr5c.cmp", familyPedigree = sparrow.famped)
sparrow.xover6a <- parse_crossovers(chrompicfile = "crimap/chr6a.cmp", familyPedigree = sparrow.famped)
#sparrow.xover6b <- parse_crossovers(chrompicfile = "crimap/chr6b.cmp", familyPedigree = sparrow.famped)
sparrow.xover7a <- parse_crossovers(chrompicfile = "crimap/chr7a.cmp", familyPedigree = sparrow.famped)
#sparrow.xover7b <- parse_crossovers(chrompicfile = "crimap/chr7b.cmp", familyPedigree = sparrow.famped)
sparrow.xover8a <- parse_crossovers(chrompicfile = "crimap/chr8a.cmp", familyPedigree = sparrow.famped)
#sparrow.xover8b <- parse_crossovers(chrompicfile = "crimap/chr8b.cmp", familyPedigree = sparrow.famped)
#sparrow.xover8c <- parse_crossovers(chrompicfile = "crimap/chr8c.cmp", familyPedigree = sparrow.famped)
sparrow.xover29a <- parse_crossovers(chrompicfile = "crimap/chr29a.cmp", familyPedigree = sparrow.famped)
sparrow.xover29b <- parse_crossovers(chrompicfile = "crimap/chr29b.cmp", familyPedigree = sparrow.famped)
#sparrow.xover29c <- parse_crossovers(chrompicfile = "crimap/chr29c.cmp", familyPedigree = sparrow.famped)



sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr5a.cmp", familyPedigree = sparrow.famped)

sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 5], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 5], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 5)), 
                      analysisID = "5a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean5 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)





sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr6a.cmp", familyPedigree = sparrow.famped)

sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 6], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 6], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 6)), 
                      analysisID = "6a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean6 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)





sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr7a.cmp", familyPedigree = sparrow.famped)

sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 7], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 7], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 7)), 
                      analysisID = "7a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean7 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)




sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr8a.cmp", familyPedigree = sparrow.famped)

sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 8], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 8], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 8)), 
                      analysisID = "8a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean8 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)




sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr9a.cmp", familyPedigree = sparrow.famped)

sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 9], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 9], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 9)), 
                      analysisID = "9a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean9 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)




sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr10a.cmp", familyPedigree = sparrow.famped)


sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

head(sparrow.doubles)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 10], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 10], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 10)), 
                      analysisID = "10a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean10 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)





sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr11a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 11], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 11], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 11)), 
                      analysisID = "11a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean11 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr12a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 12], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 12], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 12)), 
                      analysisID = "12a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean12 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)





sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr13a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 13], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 13], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 13)), 
                      analysisID = "13a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean13 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)




sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr14a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 14], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 14], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 14)), 
                      analysisID = "14a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean14 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr15a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 15], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 15], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 15)), 
                      analysisID = "15a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean15 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)





sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr17a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 17], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 17], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 17)), 
                      analysisID = "17a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean17 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)








sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr18a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 18], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 18], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 18)), 
                      analysisID = "18a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean18 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr19a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 19], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 19], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 19)), 
                      analysisID = "19a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean19 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr20a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 20], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 20], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 20)), 
                      analysisID = "20a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean20 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr21a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 21], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 21], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 21)), 
                      analysisID = "21a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean21 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)




sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr22a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 22], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 22], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 22)), 
                      analysisID = "22a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean22 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)







sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr23a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 23], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 23], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 23)), 
                      analysisID = "23a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean23 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)







sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr24a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 24], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 24], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 24)), 
                      analysisID = "24a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean24 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr25a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 25], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 25], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 25)), 
                      analysisID = "25a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean25 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)







sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr26a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 26], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 26], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 26)), 
                      analysisID = "26a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean26 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)







sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr27a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 27], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 27], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 27)), 
                      analysisID = "27a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean27 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)





sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr28a.cmp", familyPedigree = sparrow.famped)

#
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 28], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 28], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 28)), 
                      analysisID = "28a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean28 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)




###Big Chromsomes xover cleaning
###
###
###

sparrow.xover1a$data <- substr(sparrow.xover1a$data, 1, nchar(sparrow.xover1a$data)-101)
sparrow.xover1b$data <- substr(sparrow.xover1b$data, 1, nchar(sparrow.xover1b$data)-101)


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


z$data <- paste(z$data.a, z$data.b, z$data.c, sep = "")
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



z$data <- paste(z$data.a, z$data.b, z$data.c, sep = "")
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



z$data <- paste(z$data.a, z$data.b, z$data.c, sep = "")
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
#sparrow.xover4b$data <- substr(sparrow.xover4b$data, 1, nchar(sparrow.xover4b$data)-101)


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







sparrow.xover29a$data <- substr(sparrow.xover29a$data, 1, nchar(sparrow.xover29a$data)-101)
#sparrow.xover29b$data <- substr(sparrow.xover29b$data, 1, nchar(sparrow.xover29b$data)-101)

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





sparrow.xovers.clean5$No.Inf.Loci <- NULL
sparrow.xovers.clean5$First.Inf.Order <- NULL
sparrow.xovers.clean5$Last.Inf.Order<- NULL
sparrow.xovers.clean5$FATHER<- NULL
sparrow.xovers.clean5$MOTHER<- NULL
sparrow.xovers.clean5$UniqueID <- substr(sparrow.xovers.clean5$UniqueID, 4, nchar(sparrow.xovers.clean5$UniqueID))
sparrow.xovers.clean5$analysisID <- "Chr5"

sparrow.xovers.clean6$No.Inf.Loci <- NULL
sparrow.xovers.clean6$First.Inf.Order <- NULL
sparrow.xovers.clean6$Last.Inf.Order<- NULL
sparrow.xovers.clean6$FATHER<- NULL
sparrow.xovers.clean6$MOTHER<- NULL
sparrow.xovers.clean6$UniqueID <- substr(sparrow.xovers.clean6$UniqueID, 4, nchar(sparrow.xovers.clean6$UniqueID))
sparrow.xovers.clean6$analysisID <- "Chr6"

sparrow.xovers.clean7$No.Inf.Loci <- NULL
sparrow.xovers.clean7$First.Inf.Order <- NULL
sparrow.xovers.clean7$Last.Inf.Order<- NULL
sparrow.xovers.clean7$FATHER<- NULL
sparrow.xovers.clean7$MOTHER<- NULL
sparrow.xovers.clean7$UniqueID <- substr(sparrow.xovers.clean7$UniqueID, 4, nchar(sparrow.xovers.clean7$UniqueID))
sparrow.xovers.clean7$analysisID <- "Chr7"

sparrow.xovers.clean8$No.Inf.Loci <- NULL
sparrow.xovers.clean8$First.Inf.Order <- NULL
sparrow.xovers.clean8$Last.Inf.Order<- NULL
sparrow.xovers.clean8$FATHER<- NULL
sparrow.xovers.clean8$MOTHER<- NULL
sparrow.xovers.clean8$UniqueID <- substr(sparrow.xovers.clean8$UniqueID, 4, nchar(sparrow.xovers.clean8$UniqueID))
sparrow.xovers.clean8$analysisID <- "Chr8"

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
z <- merge(z, sparrow.xovers.clean22, by=c("ANIMAL","parent","Family", "RRID"))
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
  geom_histogram(binwidth = 1)+
  ggtitle("Xovers.Clean.AllChrs")
xovers.plot.all

badids_HighRecomb <- unique(z$RRID[which(z$RecombCount > 200)])
#write.table(badids_HighRecomb, file = "badids_HiRecomb.txt", row.names = FALSE, col.names = FALSE)


sparrow.xovers.clean1<-read.table("sparrow.xovers.clean.Chr1.txt", header = TRUE)
sparrow.xovers.clean2<-read.table("sparrow.xovers.clean.Chr2.txt", header = TRUE)
sparrow.xovers.clean3<-read.table("sparrow.xovers.clean.Chr3.txt", header = TRUE)
sparrow.xovers.clean4<-read.table("sparrow.xovers.clean.Chr4.txt", header = TRUE)
sparrow.xovers.clean5<-read.table("sparrow.xovers.clean.Chr5.txt", header = TRUE)
sparrow.xovers.clean6<-read.table("sparrow.xovers.clean.Chr6.txt", header = TRUE)
sparrow.xovers.clean7<-read.table("sparrow.xovers.clean.Chr7.txt", header = TRUE)
sparrow.xovers.clean8<-read.table("sparrow.xovers.clean.Chr8.txt", header = TRUE)
sparrow.xovers.clean29<-read.table("sparrow.xovers.clean.Chr29.txt", header = TRUE)

x <- sparrow.xovers.clean9
sparrow.xovers.clean9<-x

sparrow.xovers.clean9 <- sparrow.xovers.clean9[,c(1:5,11:13)]
sparrow.xovers.clean10 <- sparrow.xovers.clean10[,c(1:5,11:13)]
sparrow.xovers.clean11 <- sparrow.xovers.clean11[,c(1:5,11:13)]
sparrow.xovers.clean12 <- sparrow.xovers.clean12[,c(1:5,11:13)]
sparrow.xovers.clean13 <- sparrow.xovers.clean13[,c(1:5,11:13)]
sparrow.xovers.clean14 <- sparrow.xovers.clean14[,c(1:5,11:13)]
sparrow.xovers.clean15 <- sparrow.xovers.clean15[,c(1:5,11:13)]
sparrow.xovers.clean16 <- sparrow.xovers.clean16[,c(1:5,11:13)]
sparrow.xovers.clean17 <- sparrow.xovers.clean17[,c(1:5,11:13)]
sparrow.xovers.clean18 <- sparrow.xovers.clean18[,c(1:5,11:13)]
sparrow.xovers.clean19 <- sparrow.xovers.clean19[,c(1:5,11:13)]
sparrow.xovers.clean20 <- sparrow.xovers.clean20[,c(1:5,11:13)]
sparrow.xovers.clean21 <- sparrow.xovers.clean21[,c(1:5,11:13)]
sparrow.xovers.clean22 <- sparrow.xovers.clean22[,c(1:5,11:13)]
sparrow.xovers.clean23 <- sparrow.xovers.clean23[,c(1:5,11:13)]
sparrow.xovers.clean24 <- sparrow.xovers.clean24[,c(1:5,11:13)]
sparrow.xovers.clean25 <- sparrow.xovers.clean25[,c(1:5,11:13)]
sparrow.xovers.clean26 <- sparrow.xovers.clean26[,c(1:5,11:13)]
sparrow.xovers.clean27 <- sparrow.xovers.clean27[,c(1:5,11:13)]
sparrow.xovers.clean28 <- sparrow.xovers.clean28[,c(1:5,11:13)]


all.xovers <- rbind(sparrow.xovers.clean1,sparrow.xovers.clean2,sparrow.xovers.clean3,sparrow.xovers.clean4,
                    sparrow.xovers.clean5,sparrow.xovers.clean6,sparrow.xovers.clean7,sparrow.xovers.clean8,
                    sparrow.xovers.clean9,sparrow.xovers.clean10,sparrow.xovers.clean11,sparrow.xovers.clean12,
                    sparrow.xovers.clean13,sparrow.xovers.clean14,sparrow.xovers.clean15
                    ,sparrow.xovers.clean17,sparrow.xovers.clean18,sparrow.xovers.clean19,sparrow.xovers.clean20,
                    sparrow.xovers.clean21,sparrow.xovers.clean22,sparrow.xovers.clean23,sparrow.xovers.clean24,
                    sparrow.xovers.clean25,sparrow.xovers.clean26,sparrow.xovers.clean27,sparrow.xovers.clean28,
                    sparrow.xovers.clean29)
#save(all.xovers, file = "sparrow.xovers.clean.all.RData")
setwd("~/PhD/PhD_Repo/data/Crimap Runs/Final")
load("sparrow.xovers.clean.all.RData")
Rcount_RRID <- aggregate(all.xovers$RecombCount, by=list(Category=all.xovers$RRID),FUN=sum)
hist(Rcount_RRID$x)
table(Rcount_RRID$x)

Rcount_Chr <- aggregate(all.xovers$RecombCount, by=list(Category=all.xovers$analysisID),FUN=sum)
hist(Rcount_Chr$x)
table(Rcount_Chr$x)





