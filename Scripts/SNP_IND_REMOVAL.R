setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Ped_Map_files")
load("sparrowgen.RData")
library(crimaptools)
sparrow.abel <- sparrowgen
rm(sparrowgen)

setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/crimap/crimap/crimap")
badsnps <- read.table("badsnps.txt")
badsnps <- as.character(badsnps$V1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Smaller Chrs: SNPs which appear to cause probelms
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#CHR 9

snglSNPremove9 <- "SNPa400394"
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]

#CHR 10
snglSNPremove10 <- "SNPa432095"

sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]

#CHR 11

snglSNPremove11 <- c("SNPa95929",
"SNPa95932",
"SNPa95933",
"SNPa95941",
"SNPa95947",
"SNPa95949",
"SNPa95950",
"SNPa95953",	
"SNPa95961",
'SNPa95964',
'SNPa95966',
'SNPa95972',
'SNPa95973',
'SNPa95974',
'SNPa95976',
'SNPi28296',
'SNPa95989',
'SNPa95994',
'SNPa95998',
'SNPa96000',
'SNPa96004',
'SNPa96007',
'SNPa96009',
'SNPa96012',
'SNPa96013',
'SNPa96014',
'SNPa96019',
'SNPa96023',
'SNPa96024',
'SNPa96026',
'SNPa96029',
'SNPa96030',
'SNPa96033',
'SNPa96035',
'SNPa96050',
'SNPa96054',
'SNPa96058',
'SNPa96065',
'SNPa96066',
'SNPa96069')
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]

#Chr 13

snglSNPremove13 <- c('SNPa25768',	
'SNPa25776',
'SNPa25784',
'SNPa25785')
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]
46
#Chr 14
snglSNPremove14 <- c("SNPa472917", "SNPa472918", "SNPi4748", "SNPa472926", "SNPa472930", "SNPa472934", "SNPa472937",
"SNPa493940", "SNPa493937", "SNPi4906",   "SNPa493929", "SNPa493928", "SNPa493926", "SNPa493923",
"SNPa493917", "SNPa493916", "SNPa493912", "SNPa493906", "SNPa493905", "SNPa493902", "SNPi947",   
"SNPa493885", "SNPa493884", "SNPa493880", "SNPa493878", "SNPa493877", "SNPa493876", "SNPi31019", 
"SNPi35135", "SNPa493870", "SNPa493865", "SNPa493864", "SNPa493863", "SNPa493862", "SNPa493855",
"SNPi18511",  "SNPi24827",  "SNPa349987", "SNPa349990", "SNPa349991", "SNPa349995", "SNPa349997",
"SNPa350000", "SNPa350002", "SNPa350016", "SNPa350021", "SNPi1619",   "SNPa350027", "SNPa350029",
"SNPa350030","SNPi7180",   "SNPa350037", "SNPa350039", "SNPi5209",   "SNPa350047", "SNPa350054",
"SNPa350064", "SNPa350070", "SNPa350073", "SNPa350076", "SNPa350083", "SNPa350086", "SNPa350091",
"SNPa350099", "SNPa350100", "SNPa350104", "SNPa350107", "SNPa350109", "SNPa350126", "SNPi25290",
"SNPa350134", "SNPa350135", "SNPa350136", "SNPa350137", "SNPa350144", "SNPa350145", "SNPa350147",
"SNPa350150", "SNPa350161", "SNPa350165", "SNPa350170", "SNPa350178", "SNPa350179", "SNPa350180",
"SNPa350188", "SNPa350209", "SNPa350211", "SNPa350219", "SNPa350220", "SNPi34365",  "SNPa350222")
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]

#Chr 15


#Chr 17

snglSNPremove17 <- c("SNPa168802")
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]
#Chr 23
snglSNPremove23 <- c("SNPa397408")
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]
#Chr 26
snglSNPremove26 <- c("SNPa430919")
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]
#Chr27
snglSNPremove27 <- c("SNPa489894",
"SNPa489893",
"SNPa489892")
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]
#CHr28
snglSNPremove28 <- c("SNPa470708")
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Large Chromosomes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Chr 1
snglSNPremove1 <- "SNPi9889"
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]


#Chr 2
snglSNPremove2 <- c("SNPa328889","SNPa400196")
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]


#Chr 3
snglSNPremove3 <- c("SNPa212936","SNPa264604")
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]


#Chr 4
snglSNPremove4 <- c("SNPa221441")
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]

#Chr 5
snglSNPremove5 <- c("SNPa295169","SNPa304430","SNPa455754", "SNPa445688")
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]

#Chr 6
snglSNPremove6 <- c("SNPa471896")
sparrow.abel <- sparrow.abel[,!snpnames(sparrow.abel) %in% snglSNPremove]

SNP_spotcheck <- c(snglSNPremove1,snglSNPremove2,snglSNPremove3,snglSNPremove4,
                   snglSNPremove5,snglSNPremove6,snglSNPremove9,snglSNPremove10,
                   snglSNPremove11,snglSNPremove13,snglSNPremove14,snglSNPremove17,
                   snglSNPremove23,snglSNPremove26,snglSNPremove27,snglSNPremove28)
write.table(SNP_spotcheck, file = "SNP_spotcheck.txt", row.names = FALSE, col.names = FALSE)
master_snpremove <- c(badsnps,SNP_spotcheck,mnd_all)
write.table(master_snpremove, file = "master_snpRemove.txt", row.names = FALSE, col.names = FALSE)
