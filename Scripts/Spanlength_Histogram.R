setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/crimap/crimap/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Getting all xovers in same dataframe
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sparrow.xovers.9 <- parse_crossovers(chrompicfile = "crimap/chr9a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.10 <- parse_crossovers(chrompicfile = "crimap/chr10a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.11 <- parse_crossovers(chrompicfile = "crimap/chr11a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.12 <- parse_crossovers(chrompicfile = "crimap/chr12a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.13 <- parse_crossovers(chrompicfile = "crimap/chr13a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.14 <- parse_crossovers(chrompicfile = "crimap/chr14a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.15 <- parse_crossovers(chrompicfile = "crimap/chr15a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.17 <- parse_crossovers(chrompicfile = "crimap/chr17a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.18 <- parse_crossovers(chrompicfile = "crimap/chr18a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.19 <- parse_crossovers(chrompicfile = "crimap/chr19a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.20 <- parse_crossovers(chrompicfile = "crimap/chr20a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.21 <- parse_crossovers(chrompicfile = "crimap/chr21a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.22 <- parse_crossovers(chrompicfile = "crimap/chr22a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.23 <- parse_crossovers(chrompicfile = "crimap/chr23a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.24 <- parse_crossovers(chrompicfile = "crimap/chr24a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.25 <- parse_crossovers(chrompicfile = "crimap/chr25a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.26 <- parse_crossovers(chrompicfile = "crimap/chr26a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.27 <- parse_crossovers(chrompicfile = "crimap/chr27a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.28 <- parse_crossovers(chrompicfile = "crimap/chr28a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers.30 <- parse_crossovers(chrompicfile = "crimap/chr30a.cmp", familyPedigree = sparrow.famped)

library(dplyr)
sparrow.xovers.all <- bind_rows(sparrow.xovers.9,sparrow.xovers.10,sparrow.xovers.11,
                                                           sparrow.xovers.12,sparrow.xovers.13,sparrow.xovers.14,
                                                           sparrow.xovers.15,sparrow.xovers.17,sparrow.xovers.18,
                                                           sparrow.xovers.19,sparrow.xovers.20,sparrow.xovers.21,
                                                           sparrow.xovers.22,sparrow.xovers.23,sparrow.xovers.24,sparrow.xovers.25,
                                                           sparrow.xovers.26,sparrow.xovers.27,sparrow.xovers.28)
#setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Ped_Map_files/")
#write.table(sparrow.xovers.all, file = "xovers_smChrs_all.txt", col.names = TRUE, row.names = FALSE)

  
sparrow.doubles.all <- check_double_crossovers(parsed.xovers = sparrow.xovers.all)


physmap9 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 9], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 9], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 9)), 
                      analysisID = "9a")

physmap10 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 10], 
                        Position = map(sparrow.abel)[chromosome(sparrow.abel) == 10], 
                        Order = 1:length(which(chromosome(sparrow.abel) == 10)), 
                        analysisID = "10a")


physmap11 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 11], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 11], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 11)), 
                      analysisID = "11a")

physmap12 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 12], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 12], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 12)), 
                      analysisID = "12a")

physmap13 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 13], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 13], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 13)), 
                      analysisID = "13a")

physmap14 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 14], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 14], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 14)), 
                      analysisID = "14a")

physmap15 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 15], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 15], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 15)), 
                      analysisID = "15a")

physmap17 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 17], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 17], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 17)), 
                      analysisID = "17a")

physmap18 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 18], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 18], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 18)), 
                      analysisID = "18a")

physmap19 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 19], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 19], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 19)), 
                      analysisID = "19a")

physmap20 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 20], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 20], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 20)), 
                      analysisID = "20a")

physmap21 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 21], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 21], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 21)), 
                      analysisID = "21a")

physmap22 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 22], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 22], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 22)), 
                      analysisID = "22a")

physmap23 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 23], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 23], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 23)), 
                      analysisID = "23a")

physmap24 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 24], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 24], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 24)), 
                      analysisID = "24a")

physmap25 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 25], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 25], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 25)), 
                      analysisID = "25a")

physmap26 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 26], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 26], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 26)), 
                      analysisID = "26a")

physmap27 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 27], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 27], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 27)), 
                      analysisID = "27a")

physmap28 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 28], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 28], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 28)), 
                      analysisID = "28a")

physmap30 <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 30], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 30], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 30)), 
                      analysisID = "30a")

physmap.all <- bind_rows(physmap9,physmap10,physmap11, physmap12, physmap13,physmap14, physmap15, physmap17, physmap18, 
                         physmap19, physmap20, physmap21, physmap22, physmap23, physmap24,physmap25, physmap26, physmap27,
                         physmap28)


sparrow.doubles.all <- check_double_crossovers(parsed.xovers = sparrow.xovers.all, physical.map = physmap.all)


hist <- sparrow.doubles.all %>%
  ggplot(aes(x = SpanLength, color = Singleton))+
  geom_histogram()

hist







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Potentially slower way to get all xovers in same file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr9a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers[1:2,]
#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 9], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 9], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 9)), 
                      analysisID = "9a")
sparrow.doubles9 <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean9 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)


sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr10a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers[1:2,]
#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 10], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 10], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 10)), 
                      analysisID = "10a")
sparrow.doubles9 <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean9 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)





sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr11a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers[1:2,]
#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 11], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 11], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 11)), 
                      analysisID = "11a")
sparrow.doubles11 <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean11 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)






sparrow.xovers <- parse_crossovers(chrompicfile = "crimap/chr12a.cmp", familyPedigree = sparrow.famped)
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
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
sparrow.xovers[1:2,]
#Investigating double crossovers----
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers)

physmap <- data.frame(SNP.Name = snpnames(sparrow.abel)[chromosome(sparrow.abel) == 30], 
                      Position = map(sparrow.abel)[chromosome(sparrow.abel) == 30], 
                      Order = 1:length(which(chromosome(sparrow.abel) == 30)), 
                      analysisID = "30a")
sparrow.doubles <- check_double_crossovers(parsed.xovers = sparrow.xovers, physical.map = physmap)


sparrow.remove <- subset(sparrow.doubles, Singleton == "yes")
sparrow.xovers.clean30 <- revise_double_crossovers(parsed.xovers = sparrow.xovers, removesections = sparrow.remove)
