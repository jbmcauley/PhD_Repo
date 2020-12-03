
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


