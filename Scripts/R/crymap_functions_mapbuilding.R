library(crimaptools)
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

#Produces error, might be an issue with the sex assignments in the sparrow.abel gwaa.data file in some fashion. Needs checking.


#Characterizing recombination events:----

run_crimap_chrompic(genfile = paste("crimap/chr",zzz,"a.gen", sep = ""), crimap.path =  "C:/PathApps/crimap.exe")

}

getwd()
setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Ped_Map_files/crimap/")
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/Dblxoversremoved/GenABEL_QCs/QC1/")

Crymap(zzz= 9)
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
Crymap(zzz= 32)

setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data")
sparrow.map9 <- parse_map(mapfile = "crimap/chr9a.map")
sparrow.map10 <- parse_map(mapfile = "crimap/chr10a.map")
sparrow.map11 <- parse_map(mapfile = "crimap/chr11a.map")
sparrow.map12 <- parse_map(mapfile = "crimap/chr12a.map")
sparrow.map13 <- parse_map(mapfile = "crimap/chr13a.map")
sparrow.map14 <- parse_map(mapfile = "crimap/chr14a.map")
sparrow.map15 <- parse_map(mapfile = "crimap/chr15a.map")
sparrow.map17 <- parse_map(mapfile = "crimap/chr17a.map")
sparrow.map18 <- parse_map(mapfile = "crimap/chr18a.map")
sparrow.map19 <- parse_map(mapfile = "crimap/chr19a.map")
sparrow.map20 <- parse_map(mapfile = "crimap/chr20a.map")
sparrow.map21 <- parse_map(mapfile = "crimap/chr21a.map")
sparrow.map22 <- parse_map(mapfile = "crimap/chr22a.map")
sparrow.map23 <- parse_map(mapfile = "crimap/chr23a.map")
sparrow.map24 <- parse_map(mapfile = "crimap/chr24a.map")
sparrow.map25 <- parse_map(mapfile = "crimap/chr25a.map")
sparrow.map26 <- parse_map(mapfile = "crimap/chr26a.map")
sparrow.map27 <- parse_map(mapfile = "crimap/chr27a.map")
sparrow.map28 <- parse_map(mapfile = "crimap/chr28a.map")
sparrow.map30 <- parse_map(mapfile = "crimap/chr30a.map")
sparrow.map32 <- parse_map(mapfile = "crimap/chr32a.map")

setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data")
sparrow.cmpmap9 <- parse_map_chrompic(chrompicfile = ("crimap/chr9a.cmp"))
sparrow.cmpmap10 <- parse_map_chrompic(chrompicfile = ("crimap/chr10a.cmp"))
sparrow.cmpmap11 <- parse_map_chrompic(chrompicfile = ("crimap/chr11a.cmp"))
sparrow.cmpmap12 <- parse_map_chrompic(chrompicfile = ("crimap/chr12a.cmp"))
sparrow.cmpmap13 <- parse_map_chrompic(chrompicfile = ("crimap/chr13a.cmp"))
sparrow.cmpmap14 <- parse_map_chrompic(chrompicfile = ("crimap/chr14a.cmp"))
sparrow.cmpmap15 <- parse_map_chrompic(chrompicfile = ("crimap/chr15a.cmp"))
sparrow.cmpmap17 <- parse_map_chrompic(chrompicfile = ("crimap/chr17a.cmp"))
sparrow.cmpmap18 <- parse_map_chrompic(chrompicfile = ("crimap/chr18a.cmp"))
sparrow.cmpmap19 <- parse_map_chrompic(chrompicfile = ("crimap/chr19a.cmp"))
sparrow.cmpmap20 <- parse_map_chrompic(chrompicfile = ("crimap/chr20a.cmp"))
sparrow.cmpmap21 <- parse_map_chrompic(chrompicfile = ("crimap/chr21a.cmp"))
sparrow.cmpmap22 <- parse_map_chrompic(chrompicfile = ("crimap/chr22a.cmp"))
sparrow.cmpmap23 <- parse_map_chrompic(chrompicfile = ("crimap/chr23a.cmp"))
sparrow.cmpmap24 <- parse_map_chrompic(chrompicfile = ("crimap/chr24a.cmp"))
sparrow.cmpmap25 <- parse_map_chrompic(chrompicfile = ("crimap/chr25a.cmp"))
sparrow.cmpmap26 <- parse_map_chrompic(chrompicfile = ("crimap/chr26a.cmp"))
sparrow.cmpmap27 <- parse_map_chrompic(chrompicfile = ("crimap/chr27a.cmp"))
sparrow.cmpmap28 <- parse_map_chrompic(chrompicfile = ("crimap/chr28a.cmp"))
sparrow.cmpmap30 <- parse_map_chrompic(chrompicfile = ("crimap/chr30a.cmp"))
sparrow.cmpmap32 <- parse_map_chrompic(chrompicfile = ("crimap/chr32a.cmp"))




#Investigating double crossovers----

sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr10a.cmp", familyPedigree = sparrow.famped)


sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers10)

head(sparrow.doubles)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 10], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 10], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 10)), 
                      analysisID = "10a")
sparrow.doubles10 <- check_double_crossovers(parsed.xovers = sparrow.xovers10, physical.map = physmap10)


sparrow.remove <- subset(sparrow.doubles10, Singleton == "yes")
sparrow.xovers.clean10 <- revise_double_crossovers(parsed.xovers = sparrow.xovers10, removesections = sparrow.remove)








library(dplyr)
install.packages("ggplot2")

dev.new()
par(mfrow=c(3,7))
plot(sparrow.cmpmap9$Order,sparrow.cmpmap9$cMPosition, xlab = "", ylab = "cM SNP Pos", main = "Chr 9" , cex.main = .75)
plot(sparrow.cmpmap10$Order,sparrow.cmpmap10$cMPosition, xlab = "", ylab = "", main = "Chr 10",cex.main = .75)
plot(sparrow.cmpmap11$Order,sparrow.cmpmap11$cMPosition, xlab = "", ylab = "", main = "Chr 11",cex.main = .75)
plot(sparrow.cmpmap12$Order,sparrow.cmpmap12$cMPosition, xlab = "", ylab = "", main = "Chr 12",cex.main = .75)
plot(sparrow.cmpmap13$Order,sparrow.cmpmap13$cMPosition, xlab = "", ylab = "", main = "Chr 13",cex.main = .75)
plot(sparrow.cmpmap14$Order,sparrow.cmpmap14$cMPosition, xlab = "", ylab = "", main = "Chr 14",cex.main = .75)
plot(sparrow.cmpmap15$Order,sparrow.cmpmap15$cMPosition, xlab = "", ylab = "", main = "Chr 15",cex.main = .75)
plot(sparrow.cmpmap17$Order,sparrow.cmpmap17$cMPosition, xlab = "", ylab = "cM SNP Pos", main = "Chr 17",cex.main = .75)
plot(sparrow.cmpmap18$Order,sparrow.cmpmap18$cMPosition, xlab = "", ylab = "", main = "Chr 18",cex.main = .75)
plot(sparrow.cmpmap19$Order,sparrow.cmpmap19$cMPosition, xlab = "", ylab = "", main = "Chr 19",cex.main = .75)
plot(sparrow.cmpmap20$Order,sparrow.cmpmap20$cMPosition, xlab = "", ylab = "", main = "Chr 20",cex.main = .75)
plot(sparrow.cmpmap21$Order,sparrow.cmpmap21$cMPosition, xlab = "", ylab = "", main = "Chr 21",cex.main = .75)
plot(sparrow.cmpmap22$Order,sparrow.cmpmap22$cMPosition, xlab = "", ylab = "", main = "Chr 22",cex.main = .75)
plot(sparrow.cmpmap23$Order,sparrow.cmpmap23$cMPosition, xlab = "", ylab = "", main = "Chr 23",cex.main = .75)
plot(sparrow.cmpmap24$Order,sparrow.cmpmap24$cMPosition, xlab = "SNP Order", ylab = "cM SNP Pos", main = "Chr 24",cex.main = .75)
plot(sparrow.cmpmap25$Order,sparrow.cmpmap25$cMPosition, xlab = "SNP Order", ylab = "", main = "Chr 25",cex.main = .75)
plot(sparrow.cmpmap26$Order,sparrow.cmpmap26$cMPosition, xlab = "SNP Order", ylab = "", main = "Chr 26",cex.main = .75)
plot(sparrow.cmpmap27$Order,sparrow.cmpmap27$cMPosition, xlab = "SNP Order", ylab = "", main = "Chr 27",cex.main = .75)
plot(sparrow.cmpmap28$Order,sparrow.cmpmap28$cMPosition, xlab = "SNP Order", ylab = "", main = "Chr 28",cex.main = .75)
plot(sparrow.cmpmap30$Order,sparrow.cmpmap30$cMPosition, xlab = "SNP Order", ylab = "", main = "Chr 30",cex.main = .75)
plot(sparrow.cmpmap32$Order,sparrow.cmpmap32$cMPosition, xlab = "SNP Order", ylab = "", main = "Chr 32",cex.main = .75)
dev.off()



sparrow.map10 <- parse_map(mapfile = "crimap/chr10a.map")

dev.new()
par(mfrow=c(3,7))
plot(sparrow.map9$Order,sparrow.map9$cMPosition.Female, xlab = "", ylab = "cM SNP Pos", main = "Chr 9" , cex.main = .75)
plot(sparrow.map10$Order,sparrow.map10$cMPosition.Female, xlab = "", ylab = "", main = "Chr 10",cex.main = .75)
plot(sparrow.map11$Order,sparrow.map11$cMPosition.Female, xlab = "", ylab = "", main = "Chr 11",cex.main = .75)
plot(sparrow.map12$Order,sparrow.map12$cMPosition.Female, xlab = "", ylab = "", main = "Chr 12",cex.main = .75)
plot(sparrow.map13$Order,sparrow.map13$cMPosition.Female, xlab = "", ylab = "", main = "Chr 13",cex.main = .75)
plot(sparrow.map14$Order,sparrow.map14$cMPosition.Female, xlab = "", ylab = "", main = "Chr 14",cex.main = .75)
plot(sparrow.map15$Order,sparrow.map15$cMPosition.Female, xlab = "", ylab = "", main = "Chr 15",cex.main = .75)
plot(sparrow.map17$Order,sparrow.map17$cMPosition.Female, xlab = "", ylab = "cM SNP Pos", main = "Chr 17",cex.main = .75)
plot(sparrow.map18$Order,sparrow.map18$cMPosition.Female, xlab = "", ylab = "", main = "Chr 18",cex.main = .75)
plot(sparrow.map19$Order,sparrow.map19$cMPosition.Female, xlab = "", ylab = "", main = "Chr 19",cex.main = .75)
plot(sparrow.map20$Order,sparrow.map20$cMPosition.Female, xlab = "", ylab = "", main = "Chr 20",cex.main = .75)
plot(sparrow.map21$Order,sparrow.map21$cMPosition.Female, xlab = "", ylab = "", main = "Chr 21",cex.main = .75)
plot(sparrow.map22$Order,sparrow.map22$cMPosition.Female, xlab = "", ylab = "", main = "Chr 22",cex.main = .75)
plot(sparrow.map23$Order,sparrow.map23$cMPosition.Female, xlab = "", ylab = "", main = "Chr 23",cex.main = .75)
plot(sparrow.map24$Order,sparrow.map24$cMPosition.Female, xlab = "SNP Order", ylab = "cM SNP Pos", main = "Chr 24",cex.main = .75)
plot(sparrow.map25$Order,sparrow.map25$cMPosition.Female, xlab = "SNP Order", ylab = "", main = "Chr 25",cex.main = .75)
plot(sparrow.map26$Order,sparrow.map26$cMPosition.Female, xlab = "SNP Order", ylab = "", main = "Chr 26",cex.main = .75)
plot(sparrow.map27$Order,sparrow.map27$cMPosition.Female, xlab = "SNP Order", ylab = "", main = "Chr 27",cex.main = .75)
plot(sparrow.map28$Order,sparrow.map28$cMPosition.Female, xlab = "SNP Order", ylab = "", main = "Chr 28",cex.main = .75)
plot(sparrow.map30$Order,sparrow.map30$cMPosition.Female, xlab = "SNP Order", ylab = "", main = "Chr 30",cex.main = .75)
plot(sparrow.map32$Order,sparrow.map32$cMPosition.Female, xlab = "SNP Order", ylab = "", main = "Chr 32",cex.main = .75)
dev.off()

dev.new()
par(mfrow=c(3,7))
plot(sparrow.map9$Order,sparrow.map9$cMPosition.Male, xlab = "", ylab = "cM SNP Pos", main = "Chr 9" , cex.main = .75)
plot(sparrow.map10$Order,sparrow.map10$cMPosition.Male, xlab = "", ylab = "", main = "Chr 10",cex.main = .75)
plot(sparrow.map11$Order,sparrow.map11$cMPosition.Male, xlab = "", ylab = "", main = "Chr 11",cex.main = .75)
plot(sparrow.map12$Order,sparrow.map12$cMPosition.Male, xlab = "", ylab = "", main = "Chr 12",cex.main = .75)
plot(sparrow.map13$Order,sparrow.map13$cMPosition.Male, xlab = "", ylab = "", main = "Chr 13",cex.main = .75)
plot(sparrow.map14$Order,sparrow.map14$cMPosition.Male, xlab = "", ylab = "", main = "Chr 14",cex.main = .75)
plot(sparrow.map15$Order,sparrow.map15$cMPosition.Male, xlab = "", ylab = "", main = "Chr 15",cex.main = .75)
plot(sparrow.map17$Order,sparrow.map17$cMPosition.Male, xlab = "", ylab = "cM SNP Pos", main = "Chr 17",cex.main = .75)
plot(sparrow.map18$Order,sparrow.map18$cMPosition.Male, xlab = "", ylab = "", main = "Chr 18",cex.main = .75)
plot(sparrow.map19$Order,sparrow.map19$cMPosition.Male, xlab = "", ylab = "", main = "Chr 19",cex.main = .75)
plot(sparrow.map20$Order,sparrow.map20$cMPosition.Male, xlab = "", ylab = "", main = "Chr 20",cex.main = .75)
plot(sparrow.map21$Order,sparrow.map21$cMPosition.Male, xlab = "", ylab = "", main = "Chr 21",cex.main = .75)
plot(sparrow.map22$Order,sparrow.map22$cMPosition.Male, xlab = "", ylab = "", main = "Chr 22",cex.main = .75)
plot(sparrow.map23$Order,sparrow.map23$cMPosition.Male, xlab = "", ylab = "", main = "Chr 23",cex.main = .75)
plot(sparrow.map24$Order,sparrow.map24$cMPosition.Male, xlab = "SNP Order", ylab = "cM SNP Pos", main = "Chr 24",cex.main = .75)
plot(sparrow.map25$Order,sparrow.map25$cMPosition.Male, xlab = "SNP Order", ylab = "", main = "Chr 25",cex.main = .75)
plot(sparrow.map26$Order,sparrow.map26$cMPosition.Male, xlab = "SNP Order", ylab = "", main = "Chr 26",cex.main = .75)
plot(sparrow.map27$Order,sparrow.map27$cMPosition.Male, xlab = "SNP Order", ylab = "", main = "Chr 27",cex.main = .75)
plot(sparrow.map28$Order,sparrow.map28$cMPosition.Male, xlab = "SNP Order", ylab = "", main = "Chr 28",cex.main = .75)
plot(sparrow.map30$Order,sparrow.map30$cMPosition.Male, xlab = "SNP Order", ylab = "", main = "Chr 30",cex.main = .75)
plot(sparrow.map32$Order,sparrow.map32$cMPosition.Male, xlab = "SNP Order", ylab = "", main = "Chr 32",cex.main = .75)
dev.off()



#BP vs cM

dev.new()
par(mfrow=c(3,7))
plot(sparrow.map9$BP,sparrow.map9$cMPosition, xlab = "", ylab = "cM SNP Pos", main = "Chr 9" , cex.main = .75)
plot(sparrow.map10$BP,sparrow.map10$cMPosition, xlab = "", ylab = "", main = "Chr 10",cex.main = .75)
plot(sparrow.map11$BP,sparrow.map11$cMPosition, xlab = "", ylab = "", main = "Chr 11",cex.main = .75)
plot(sparrow.map12$BP,sparrow.map12$cMPosition, xlab = "", ylab = "", main = "Chr 12",cex.main = .75)
plot(sparrow.map13$BP,sparrow.map13$cMPosition, xlab = "", ylab = "", main = "Chr 13",cex.main = .75)
plot(sparrow.map14$BP,sparrow.map14$cMPosition, xlab = "", ylab = "", main = "Chr 14",cex.main = .75)
plot(sparrow.map15$BP,sparrow.map15$cMPosition, xlab = "", ylab = "", main = "Chr 15",cex.main = .75)
plot(sparrow.map17$BP,sparrow.map17$cMPosition, xlab = "", ylab = "cM SNP Pos", main = "Chr 17",cex.main = .75)
plot(sparrow.map18$BP,sparrow.map18$cMPosition, xlab = "", ylab = "", main = "Chr 18",cex.main = .75)
plot(sparrow.map19$BP,sparrow.map19$cMPosition, xlab = "", ylab = "", main = "Chr 19",cex.main = .75)
plot(sparrow.map20$BP,sparrow.map20$cMPosition, xlab = "", ylab = "", main = "Chr 20",cex.main = .75)
plot(sparrow.map21$BP,sparrow.map21$cMPosition, xlab = "", ylab = "", main = "Chr 21",cex.main = .75)
plot(sparrow.map22$BP,sparrow.map22$cMPosition, xlab = "", ylab = "", main = "Chr 22",cex.main = .75)
plot(sparrow.map23$BP,sparrow.map23$cMPosition, xlab = "", ylab = "", main = "Chr 23",cex.main = .75)
plot(sparrow.map24$BP,sparrow.map24$cMPosition, xlab = "BP Pos", ylab = "cM SNP Pos", main = "Chr 24",cex.main = .75)
plot(sparrow.map25$BP,sparrow.map25$cMPosition, xlab = "BP Pos", ylab = "", main = "Chr 25",cex.main = .75)
plot(sparrow.map26$BP,sparrow.map26$cMPosition, xlab = "BP Pos", ylab = "", main = "Chr 26",cex.main = .75)
plot(sparrow.map27$BP,sparrow.map27$cMPosition, xlab = "BP Pos", ylab = "", main = "Chr 27",cex.main = .75)
plot(sparrow.map28$BP,sparrow.map28$cMPosition, xlab = "BP Pos", ylab = "", main = "Chr 28",cex.main = .75)
plot(sparrow.map30$BP,sparrow.map30$cMPosition, xlab = "BP Pos", ylab = "", main = "Chr 30",cex.main = .75)
plot(sparrow.map32$BP,sparrow.map32$cMPosition, xlab = "BP Pos", ylab = "", main = "Chr 32",cex.main = .75)
dev.off()




dev.new()
par(mfrow=c(3,7))
plot(sparrow.map9$Order,sparrow.map9$cMPosition, xlab = "", ylab = "cM SNP Pos", main = "Chr 9" , cex.main = .75)
plot(sparrow.map10$Order,sparrow.map10$cMPosition, xlab = "", ylab = "", main = "Chr 10",cex.main = .75)
plot(sparrow.map11$Order,sparrow.map11$cMPosition, xlab = "", ylab = "", main = "Chr 11",cex.main = .75)
plot(sparrow.map12$Order,sparrow.map12$cMPosition, xlab = "", ylab = "", main = "Chr 12",cex.main = .75)
plot(sparrow.map13$Order,sparrow.map13$cMPosition, xlab = "", ylab = "", main = "Chr 13",cex.main = .75)
plot(sparrow.map14$Order,sparrow.map14$cMPosition, xlab = "", ylab = "", main = "Chr 14",cex.main = .75)
plot(sparrow.map15$Order,sparrow.map15$cMPosition, xlab = "", ylab = "", main = "Chr 15",cex.main = .75)
plot(sparrow.map17$Order,sparrow.map17$cMPosition, xlab = "", ylab = "cM SNP Pos", main = "Chr 17",cex.main = .75)
plot(sparrow.map18$Order,sparrow.map18$cMPosition, xlab = "", ylab = "", main = "Chr 18",cex.main = .75)
plot(sparrow.map19$Order,sparrow.map19$cMPosition, xlab = "", ylab = "", main = "Chr 19",cex.main = .75)
plot(sparrow.map20$Order,sparrow.map20$cMPosition, xlab = "", ylab = "", main = "Chr 20",cex.main = .75)
plot(sparrow.map21$Order,sparrow.map21$cMPosition, xlab = "", ylab = "", main = "Chr 21",cex.main = .75)
plot(sparrow.map22$Order,sparrow.map22$cMPosition, xlab = "", ylab = "", main = "Chr 22",cex.main = .75)
plot(sparrow.map23$Order,sparrow.map23$cMPosition, xlab = "", ylab = "", main = "Chr 23",cex.main = .75)
plot(sparrow.map24$Order,sparrow.map24$cMPosition, xlab = "SNP Order", ylab = "cM SNP Pos", main = "Chr 24",cex.main = .75)
plot(sparrow.map25$Order,sparrow.map25$cMPosition, xlab = "SNP Order", ylab = "", main = "Chr 25",cex.main = .75)
plot(sparrow.map26$Order,sparrow.map26$cMPosition, xlab = "SNP Order", ylab = "", main = "Chr 26",cex.main = .75)
plot(sparrow.map27$Order,sparrow.map27$cMPosition, xlab = "SNP Order", ylab = "", main = "Chr 27",cex.main = .75)
plot(sparrow.map28$Order,sparrow.map28$cMPosition, xlab = "SNP Order", ylab = "", main = "Chr 28",cex.main = .75)
plot(sparrow.map30$Order,sparrow.map30$cMPosition, xlab = "SNP Order", ylab = "", main = "Chr 30",cex.main = .75)
plot(sparrow.map32$Order,sparrow.map32$cMPosition, xlab = "SNP Order", ylab = "", main = "Chr 32",cex.main = .75)
dev.off()


sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr6a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 6], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 6], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 6)), 
                      analysisID = "6a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean6 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)




sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr9a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
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

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 11], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 11], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 11)), 
                      analysisID = "11a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean11 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr12a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 12], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 12], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 12)), 
                      analysisID = "12a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean12 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)





sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr13a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 13], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 13], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 13)), 
                      analysisID = "13a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean13 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)




sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr14a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 14], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 14], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 14)), 
                      analysisID = "14a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean14 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr15a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 15], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 15], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 15)), 
                      analysisID = "15a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean15 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)





sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr17a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 17], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 17], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 17)), 
                      analysisID = "17a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean17 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)








sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr18a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 18], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 18], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 18)), 
                      analysisID = "18a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean18 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr19a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 19], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 19], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 19)), 
                      analysisID = "19a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean19 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr20a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 20], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 20], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 20)), 
                      analysisID = "20a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean20 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr21a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 21], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 21], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 21)), 
                      analysisID = "21a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean21 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)




sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr22a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 22], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 22], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 22)), 
                      analysisID = "22a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean22 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)







sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr23a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 23], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 23], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 23)), 
                      analysisID = "23a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean23 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)







sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr24a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 24], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 24], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 24)), 
                      analysisID = "24a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean24 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr25a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 25], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 25], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 25)), 
                      analysisID = "25a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean25 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)







sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr26a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 26], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 26], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 26)), 
                      analysisID = "26a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean26 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)







sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr27a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 27], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 27], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 27)), 
                      analysisID = "27a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean27 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)





sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr28a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 28], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 28], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 28)), 
                      analysisID = "28a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean28 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr30a.cmp", familyPedigree = sparrow.famped)

#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 30], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 30], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 30)), 
                      analysisID = "30a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean30 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)





library(dplyr)

RecombSum_Sm_Chrs <- bind_rows(sparrow.xovers.clean9,sparrow.xovers.clean10.RData,sparrow.xovers.clean11,
          sparrow.xovers.clean12,sparrow.xovers.clean13, sparrow.xovers.clean14,
          sparrow.xovers.clean15,sparrow.xovers.clean17,sparrow.xovers.clean18,
          sparrow.xovers.clean19,sparrow.xovers.clean20,sparrow.xovers.clean21,
          sparrow.xovers.clean22,sparrow.xovers.clean23,sparrow.xovers.clean24,sparrow.xovers.clean25,
          sparrow.xovers.clean26,sparrow.xovers.clean27,sparrow.xovers.clean28) %>%
  group_by(Family) %>%
  summarise_at("RecombCount", sum)


[18:23] JOHNSTON Susan


RecombSum_Sm_Chrs <- bind_rows(sparrow.xovers.clean9,sparrow.xovers.clean10.RData,sparrow.xovers.clean11,
                               sparrow.xovers.clean12,sparrow.xovers.clean13, sparrow.xovers.clean14,
                               sparrow.xovers.clean15,sparrow.xovers.clean17,sparrow.xovers.clean18,
                               sparrow.xovers.clean19,sparrow.xovers.clean20,sparrow.xovers.clean21,
                               sparrow.xovers.clean22,sparrow.xovers.clean23,sparrow.xovers.clean24,sparrow.xovers.clean25,
                               sparrow.xovers.clean26,sparrow.xovers.clean27,sparrow.xovers.clean28) %>%
  group_by(Family) %>%
  summarise(TotalRecombCount = sum(RecombCount))
write.table(RecombSum_Sm_Chrs, file = "RecombSum_Sm_Chrs.txt",row.names = FALSE, col.names = TRUE)          
            

 
sparrow.xovers.sml.chrs <- bind_rows(sparrow.xovers.clean9,sparrow.xovers.clean10.RData,sparrow.xovers.clean11,
          sparrow.xovers.clean12,sparrow.xovers.clean13, sparrow.xovers.clean14,
          sparrow.xovers.clean15,sparrow.xovers.clean17,sparrow.xovers.clean18,
          sparrow.xovers.clean19,sparrow.xovers.clean20,sparrow.xovers.clean21,
          sparrow.xovers.clean22,sparrow.xovers.clean23,sparrow.xovers.clean24,sparrow.xovers.clean25,
          sparrow.xovers.clean26,sparrow.xovers.clean27,sparrow.xovers.clean28)

write.table(sparrow.xovers.sml.chrs, file = "xovers.sml.chrs.txt", row.names = FALSE, col.names = TRUE)
