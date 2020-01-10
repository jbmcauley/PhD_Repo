#Path Problems Currently Exist on University PC: Must work with Personal PC for time being. Get FAMILY Variable working in mean time.

rm(list=ls())

library(devtools)
install_github("susjoh/crimaptools", force = TRUE)
setwd("C:/Users/s1945757/PhD_Repo/Cri_Map/crimaptools-master/crimaptools-master/")
getwd()
#Creating Cri-map input files-----
library(crimaptools)
library(GenABEL)

data(deer)
create_crimap_input(gwaa.data = deer.abel, familyPedigree = deer.famped, analysisID = "3a", chr = 3, outdir = "crimap", clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr3a.gen", crimap.path = "C:/PathApps/crimap.exe")
#Function has not produced the .loc .par and .dat files suggested in the tutorial
#Using the terminal does work following the crimapinput1 preiously created (n,n,n,n,7,y,y) to generate .loc, .par, and .dat files... mention to Susan a potential issue in crimaptools
#Alternatively the problem is due to miscommunication between R and Path, have had Path issues in past with Plink...

dir("crimap")

parse_mend_err(prefile = "crimap/chr3a.pre", genfile = "crimap/chr3a.gen", familyPedigree = deer.famped)
read.table("crimap/chr3a.mnd", header = T)
#No Mendelian errors detected, in tutorial "2" are listed... Data already cleaned?
#Next line only used if mendialian erros are present in order to mask them in the .gen file.
create_crimap_input (gwaa.data = deer.abel, 
                     familyPedigree = deer.famped, 
                     analysisID = "3a", 
                     chr = 3, 
                     outdir = "crimap", 
                     clear.existing.analysisID = TRUE, 
                     use.mnd = TRUE)


#If mendelian errors detected you will need to rerun the prepare function. Once again, might need to be done in terminal rather than crimap-tools
run_crimap_prepare(genfile = "crimap/chr3a.gen", crimap.path = "C:/PathApps/crimap.exe") 
parse_mend_err(prefile = "crimap/chr3a.pre", genfile = "crimap/chr3a.gen", familyPedigree = deer.famped)
dir("crimap")
#Repeat the process of looking for mendelian error until there are none.


#Build a Linkage Map:----

run_crimap_map(genfile = "crimap/chr3a.gen", crimap.path = "C:/PathApps/crimap.exe")
#.map file generated has size "0 B" This is a sign that something is not working in generation of the file. Will need to work with regular functions 
#in terminal.
dir("crimap")
deer.map <- parse_map(mapfile = "crimap/chr3a.map")
#Produces error, might be an issue with the sex assignments in the deer.abel gwaa.data file in some fashion. Needs checking.


#Characterizing recombination events:----

run_crimap_chrompic(genfile = "crimap/chr3a.gen", crimap.path =  "C:/PathApps/crimap.exe")
deer.cmpmap <- parse_map_chrompic(chrompicfile = "crimap/chr3a.cmp")
head(deer.cmpmap)

deer.xovers <- parse_crossovers(chrompicfile = "crimap/chr3a.cmp", familyPedigree = deer.famped)
deer.xovers[1:2,]


#Investigating double crossovers----
deer.doubles <- check_double_crossovers(parsed.xovers = deer.xovers)

head(deer.doubles)

physmap <- data.frame(SNP.Name = snpnames(deer.abel)[chromosome(deer.abel) == 3], 
                      Position = map(deer.abel)[chromosome(deer.abel) == 3], 
                      Order = 1:length(which(chromosome(deer.abel) == 3)), 
                      analysisID = "3a")
deer.doubles <- check_double_crossovers(parsed.xovers = deer.xovers, physical.map = physmap)


deer.remove <- subset(deer.doubles, Singleton == "yes")
deer.xovers.clean <- revise_double_crossovers(parsed.xovers = deer.xovers, removesections = deer.remove)

deer.xovers.clean
