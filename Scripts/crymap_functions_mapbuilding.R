
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


setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/crimap")

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

sparrow.xovers9 <- parse_crossovers(chrompicfile = "crimap/chr9a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers9[1:2,]


#Investigating double crossovers----
sparrow.doubles9 <- check_double_crossovers(parsed.xovers = sparrow.xovers9)

head(sparrow.doubles)

physmap9 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 9], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 9], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 9)), 
                      analysisID = "9a")
sparrow.doubles9 <- check_double_crossovers(parsed.xovers = sparrow.xovers9, physical.map = physmap9)


sparrow.remove <- subset(sparrow.doubles9, Singleton == "yes")
saveRDS((revise_double_crossovers(parsed.xovers = sparrow.xovers9, removesections = sparrow.remove)), file = "sparrow_xovers_clean9.RData")

sparrow.xovers.clean9






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

