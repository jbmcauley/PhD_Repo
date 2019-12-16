rm(list=ls())
library(GenABEL)
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/")
dir()


famfile <- read.table("Helgeland_01_2018.fam", stringsAsFactors = F)
famfile <- famfile[,2:5]
names(famfile) <- c("id", "Father", "Mother", "sex")


write.table(famfile, "Helgeland_01_2018.phe", row.names = F, sep = "\t", quote = F)

#~~ Create map file

mapfile <- read.table("Helgeland_01_2018.map")
head(mapfile)
write.table(mapfile[,c(1, 2, 4)], "Helgeland_01_2018.genabelmap", row.names = F, col.names = F, quote = F, sep = "\t")
table.check <- read.table("Helgeland_01_2018.genabelmap")
rm(list=ls())
#~~ Make GenAbel files

convert.snp.ped(pedfile = "Helgeland_01_2018.ped", 
                mapfile = "Helgeland_01_2018.genabelmap",
                outfile = "Helgeland_01_2018.genabel",
                strand = "u", bcast = 10000, traits = 1, mapHasHeaderLine = F)

###### PROBLEM HERE - FIX WITH LETTERS????

sparrowgen.Helgeland <- load.gwaa.data(phenofile = "Helgeland_01_2018.phe",
                             genofile  = "Helgeland_01_2018.genabel")



save(sparrowgen.Helgeland, file = "sparrowgen_Helgeland_01_2018.RData")


