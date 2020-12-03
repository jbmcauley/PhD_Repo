
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






###
###Graphs
###
###
###

setwd("C:/Users/s1945757/PhD_Repo/Crimap Outputs Post Axiom Cleaning")


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
sparrow.map1[,c(5,6,7,8,10,11,12)] <- NULL

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
sparrow.map2[,c(5,6,7,8,10,11,12)] <- NULL

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
sparrow.map3[,c(5,6,7,8,10,11,12)] <- NULL

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
sparrow.map4[,c(5,6,7,8,10,11,12)] <- NULL

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
sparrow.map29[,c(5,6,7,8,10,11,12)] <- NULL
rm(sparrow.map29a)
rm(sparrow.map29b)
rm(sparrow.map29c)


#Ch 5
sparrow.map5 <- parse_map(mapfile = "crimap/chr5a.map")
sparrow.map5 <- merge(sparrow.map5, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map5[,c(5,6,7,8,10,11,12)] <- NULL

#Ch 6
sparrow.map6 <- parse_map(mapfile = "crimap/chr6a.map")
sparrow.map6 <- merge(sparrow.map6, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map6[,c(5,6,7,8,10,11,12)] <- NULL

#Ch 7
sparrow.map7 <- parse_map(mapfile = "crimap/chr7a.map")
sparrow.map7 <- merge(sparrow.map7, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map7[,c(5,6,7,8,10,11,12)] <- NULL

#Ch 8
sparrow.map8 <- parse_map(mapfile = "crimap/chr8a.map")
sparrow.map8 <- merge(sparrow.map8, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map8[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 9
sparrow.map9 <- parse_map(mapfile = "crimap/chr9a.map")
sparrow.map9 <- merge(sparrow.map9, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map9[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 10
sparrow.map10 <- parse_map(mapfile = "crimap/chr10a.map")
sparrow.map10 <- merge(sparrow.map10, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map10[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 11
sparrow.map11 <- parse_map(mapfile = "crimap/chr11a.map")
sparrow.map11 <- merge(sparrow.map11, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map11[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 12
sparrow.map12 <- parse_map(mapfile = "crimap/chr12a.map")
sparrow.map12 <- merge(sparrow.map12, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map12[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 13
sparrow.map13 <- parse_map(mapfile = "crimap/chr13a.map")
sparrow.map13 <- merge(sparrow.map13, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map13[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 14
sparrow.map14 <- parse_map(mapfile = "crimap/chr14a.map")
sparrow.map14 <- merge(sparrow.map14, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map14[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 15
sparrow.map15 <- parse_map(mapfile = "crimap/chr15a.map")
sparrow.map15 <- merge(sparrow.map15, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map15[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 17
sparrow.map17 <- parse_map(mapfile = "crimap/chr17a.map")
sparrow.map17 <- merge(sparrow.map17, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map17[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 18
sparrow.map18 <- parse_map(mapfile = "crimap/chr18a.map")
sparrow.map18 <- merge(sparrow.map18, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map18[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 19
sparrow.map19 <- parse_map(mapfile = "crimap/chr19a.map")
sparrow.map19 <- merge(sparrow.map19, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map19[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 20
sparrow.map20 <- parse_map(mapfile = "crimap/chr20a.map")
sparrow.map20 <- merge(sparrow.map20, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map20[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 21
sparrow.map21 <- parse_map(mapfile = "crimap/chr21a.map")
sparrow.map21 <- merge(sparrow.map21, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map21[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 22
sparrow.map22 <- parse_map(mapfile = "crimap/chr22a.map")
sparrow.map22 <- merge(sparrow.map22, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map22[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 23
sparrow.map23 <- parse_map(mapfile = "crimap/chr23a.map")
sparrow.map23 <- merge(sparrow.map23, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map23[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 24
sparrow.map24 <- parse_map(mapfile = "crimap/chr24a.map")
sparrow.map24 <- merge(sparrow.map24, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map24[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 26
sparrow.map26 <- parse_map(mapfile = "crimap/chr26a.map")
sparrow.map26 <- merge(sparrow.map26, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map26[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 27
sparrow.map27 <- parse_map(mapfile = "crimap/chr27a.map")
sparrow.map27 <- merge(sparrow.map27, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map27[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 28
sparrow.map28 <- parse_map(mapfile = "crimap/chr28a.map")
sparrow.map28 <- merge(sparrow.map28, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map28[,c(5,6,7,8,10,11,12)] <- NULL

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
sparrow.map26 <- sparrow.map26[order(sparrow.map26$Order),]
sparrow.map27 <- sparrow.map27[order(sparrow.map27$Order),]
sparrow.map28 <- sparrow.map28[order(sparrow.map28$Order),]
sparrow.map29 <- sparrow.map29[order(sparrow.map29$Order),]



sparrow.map.all <- rbind(sparrow.map1,sparrow.map2,sparrow.map3,sparrow.map4,sparrow.map5,sparrow.map6,sparrow.map7,sparrow.map8,
                         sparrow.map9,sparrow.map10,sparrow.map11,sparrow.map12,sparrow.map13,sparrow.map14,sparrow.map15,sparrow.map17,
                         sparrow.map18,sparrow.map19,sparrow.map20,sparrow.map21,sparrow.map22,sparrow.map23,sparrow.map24,sparrow.map26,
                         sparrow.map27,sparrow.map28,sparrow.map29)
sparrow.map.all$Chr <- factor(sparrow.map.all$Chr, levels=c("Chr 1","Chr 2","Chr 3","Chr 4","Chr 5","Chr 6","Chr 7","Chr 8","Chr 9",
                                                            "Chr 10","Chr 11","Chr 12","Chr 13","Chr 14","Chr 15","Chr 17","Chr 18",
                                                            "Chr 19","Chr 20","Chr 21","Chr 22","Chr 23","Chr 24","Chr 26","Chr 27",
                                                            "Chr 28","Chr 29"))
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

save(final_snpremove, file = "SNP_Logical_Vectors.RData")

#Obvious Chrs to apply filtering to: 1,2,3,4,5,6,7,8,10,11,12,14,15


load("SNP_Logical_Vectors.RData")
