rm(list=ls())

library(devtools)
install_github("susjoh/crimaptools")
setwd("C:/Users/s1945757/PhD_Repo/Cri_Map/crimaptools-master/crimaptools-master/")
getwd()


library(crimaptools)
library(GenABEL)

data(deer)
create_crimap_input(gwaa.data = deer.abel, familyPedigree = deer.famped, analysisID = "3a", chr = 3, outdir = "crimap", clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr3a.gen")
#Function has not produced the .loc .par and .dat files suggested in the tutorial
dir("crimap")

parse_mend_err(prefile = "crimap/chr3a.pre", genfile = "crimap/chr3a.gen", familyPedigree = deer.famped)

run_crimap_map(genfile = "crimap/chr3a.gen") 
dir("crimap")

deer.map <- parse_map(mapfile = "crimap/chr3a.map") 
head(deer.map)
