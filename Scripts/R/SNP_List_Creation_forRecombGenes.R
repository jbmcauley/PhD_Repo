library(crimaptools)
library(GenABEL)
library(GenABEL.data)
library(dplyr)

setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Ped_Map_files")
load("sparrowgen.RData")

chr1 <- sparrowgen[,sparrowgen@gtdata@chromosome == 1]
chr2 <- sparrowgen[,sparrowgen@gtdata@chromosome == 2]
chr3 <- sparrowgen[,sparrowgen@gtdata@chromosome == 3]
chr4 <- sparrowgen[,sparrowgen@gtdata@chromosome == 4] 
chr5 <- sparrowgen[,sparrowgen@gtdata@chromosome == 5]
chr6 <- sparrowgen[,sparrowgen@gtdata@chromosome ==6]
chr8 <- sparrowgen[,sparrowgen@gtdata@chromosome == 8]
chr10 <- sparrowgen[,sparrowgen@gtdata@chromosome == 10]
chr12 <- sparrowgen[,sparrowgen@gtdata@chromosome == 12]
chr13 <- sparrowgen[,sparrowgen@gtdata@chromosome == 13]
chr14 <- sparrowgen[,sparrowgen@gtdata@chromosome == 14]
chr20 <- sparrowgen[,sparrowgen@gtdata@chromosome == 20]
chr27 <- sparrowgen[,sparrowgen@gtdata@chromosome == 27]
chr29 <- sparrowgen[,sparrowgen@gtdata@chromosome == 29]

#CBLL1; Chr1a: 12380129	12383500
y <- chr29[,between(chr29@gtdata@map, 12380129, 12383500)]
x<-y@gtdata@map
CBLL1 <- c(names(which(chr29@gtdata@map >= 12383500)[1:100]),
           names(which(chr29@gtdata@map <= 12380129)[2884:2984]))


#RECQL; chr1A;	69730714	69742374
y <- chr29[,between(chr29@gtdata@map, 69730714, 69742374)]
x<-y@gtdata@map
RECQL <- c(names(which(chr29@gtdata@map >= 69742374)[1:100]),names(x),
         names(which(chr29@gtdata@map <= 69730714)[16154:16254]))


#MYC chr2		8901363	8903766
y <- chr2[,between(chr2@gtdata@map, 8901363, 8903766)]
MYC <- c(names(which(chr2@gtdata@map >= 8903766)[1:100]),
            names(which(chr2@gtdata@map <= 8901363)[741:841]))



#MAEA; chr4		63332471	63333787
y <- chr4[,between(chr4@gtdata@map, 63332471, 63333787)]
MAEA <- c(names(which(chr4@gtdata@map >= 63333787)[1:100]),
         names(which(chr4@gtdata@map <= 63332471)[15687:15787]))

#FGFRL1; chr4		63927199	63933410
y <- chr4[,between(chr4@gtdata@map, 63927199, 63933410)]
x<-y@gtdata@map
#SNPa454530; 63928755
which(chr4@gtdata@map == 63928755)
FGFRL1 <- c(names(which(chr4@gtdata@map >= 63933410)[1:100]),names(which(chr4@gtdata@map == 63928755)),
           names(which(chr4@gtdata@map <= 63927199)[15822:15922]))

#MRPL35; chr4		64374204	64377025
y <- chr4[,between(chr4@gtdata@map, 64374204, 64377025)]
names(which(chr4@gtdata@map >= 64377025)[1:100])
names(which(chr4@gtdata@map <= 64374204)[15885:15985])
MRPL35 <- c(names(which(chr4@gtdata@map >= 64377025)[1:100]),
            names(which(chr4@gtdata@map <= 64374204)[15885:15985]))

chr4@gtdata@map[which(chr4@gtdata@snpnames == "SNPa454760")]




#GNAI2; chr12		5554811	5560682

y <- chr12[,between(chr12@gtdata@map, 5554811, 5560682)]
names(which(chr12@gtdata@map >= 64377025)[1:100])
names(which(chr12@gtdata@map <= 64374204)[15885:15985])
GNAI2 <- c(names(which(chr12@gtdata@map >= 5560682)[1:100]),
            names(which(chr12@gtdata@map <= 5554811)[798:898]))

#TDG; chr1A		51372033	51378897
y <- chr29[,between(chr29@gtdata@map, 51372033, 51378897)]
x<-y@gtdata@map
TDG <- c(names(which(chr29@gtdata@map >= 51378897)[1:100]),names(x),
          names(which(chr29@gtdata@map <= 51372033)[12242:12342]))




#NFYB; chr1A		51321757	51329193
y <- chr29[,between(chr29@gtdata@map, 51321757, 51329193)]
x<-y@gtdata@map
NFYB <- c(names(which(chr29@gtdata@map >= 51329193)[1:100]),names(x),
          names(which(chr29@gtdata@map <= 51321757)[12225:12325]))





#RPL24; chr1		16311987	16315717

y <- chr1[,between(chr1@gtdata@map, 16311987, 16315717)]

RPL24 <- c(names(which(chr1@gtdata@map >= 16315717)[1:100]),
           names(which(chr1@gtdata@map <= 16311987)[3408:3508]))




#IMPG2; chr1		16165948	16217147
y <- chr1[,between(chr1@gtdata@map, 16165948, 16217147)]
x<-y@gtdata@map
IMPG2 <- c(names(which(chr1@gtdata@map >= 16217147)[1:100]),names(x),
           names(which(chr1@gtdata@map <= 16165948)[3369:3469]))



#RGCC; chr1		54599781	54612288
y <- chr1[,between(chr1@gtdata@map, 54599781, 54612288)]
x<-y@gtdata@map
RGCC <- c(names(which(chr1@gtdata@map >= 54612288)[1:100]),names(x),
           names(which(chr1@gtdata@map <= 54599781)[12185:12285]))


#WBP4; chr1		54711630	54724061
y <- chr1[,between(chr1@gtdata@map, 54711630, 54724061)]
x<-y@gtdata@map
WBP4 <- c(names(which(chr1@gtdata@map >= 54724061)[1:100]), names(x),
           names(which(chr1@gtdata@map <= 54711630)[12214:12314]))



#LECT1; chr1		54851438	54864276
y <- chr1[,between(chr1@gtdata@map, 54851438, 54864276)]
x<-y@gtdata@map
LECT1 <- c(names(which(chr1@gtdata@map >= 54864276)[1:100]),names(x),
           names(which(chr1@gtdata@map <= 54851438)[12241:12341]))


#CCNA1; chr1		48654156	48659220
y <- chr1[,between(chr1@gtdata@map, 48654156, 48659220)]
x<-y@gtdata@map
CCNA1 <- c(names(which(chr1@gtdata@map >= 48659220)[1:100]),names(x),
           names(which(chr1@gtdata@map <= 48654156)[10563:10663]))



#ANKH; chr2		64612378	64617461
y <- chr2[,between(chr2@gtdata@map, 64612378, 64617461)]
x<-y@gtdata@map
ANKH <- c(names(which(chr2@gtdata@map >= 64617461)[1:100]),names(x),
           names(which(chr2@gtdata@map <= 64612378)[15522:15622]))



#PRIM2; chr3		87724070	87749879
y <- chr3[,between(chr3@gtdata@map, 87724070, 87749879)]
x<-y@gtdata@map
PRIM2 <- c(names(which(chr3@gtdata@map >= 87749879)[1:100]),names(x),
          names(which(chr3@gtdata@map <= 87724070)[20358:20458]))



#EPGN; chr4		2633114	2639966
y <- chr4[,between(chr4@gtdata@map, 2633114, 2639966)]
x<-y@gtdata@map
EPGN <- c(names(which(chr4@gtdata@map >= 2639966)[1:100]),
          names(which(chr4@gtdata@map <= 2633114)[273:373]))



#EREG; chr4		2606133	2611444
y <- chr4[,between(chr4@gtdata@map, 2606133, 2611444)]
x<-y@gtdata@map
EREG <- c(names(which(chr4@gtdata@map >= 2611444)[1:100]),names(x),
          names(which(chr4@gtdata@map <= 2606133)[271:371]))



#IGF2; chr5		14489308	14490217
y <- chr5[,between(chr5@gtdata@map, 14489308, 14490217)]
x<-y@gtdata@map
IGF2 <- c(names(which(chr5@gtdata@map >= 14490217)[1:100]),
          names(which(chr5@gtdata@map <= 14489308)[3189:3289]))



#RAG1; chr5		18301496	18304618
y <- chr5[,between(chr5@gtdata@map, 18301496, 18304618)]
x<-y@gtdata@map
RAG1 <- c(names(which(chr5@gtdata@map >= 18304618)[1:100]),
          names(which(chr5@gtdata@map <= 18301496)[3945:4045]))



#RAG2; chr5		18312675	18314282
y <- chr5[,between(chr5@gtdata@map, 18312675, 18314282)]
x<-y@gtdata@map
RAG2 <- c(names(which(chr5@gtdata@map >= 18314282)[1:100]),
          names(which(chr5@gtdata@map <= 18312675)[3945:4045]))



#RPS17L – Missing (alternative RPS17):chr10		3155373	3157987
y <- chr10[,between(chr10@gtdata@map, 3155373, 3157987)]
x<-y@gtdata@map
RPS17 <- c(names(which(chr10@gtdata@map >= 3157987)[1:100]),
          names(which(chr10@gtdata@map <= 3155373)[281:381]))



#YWHAB; chr20		8970492	8976807
y <- chr20[,between(chr20@gtdata@map, 8970492, 8976807)]
x<-y@gtdata@map
YWHAB <- c(names(which(chr20@gtdata@map >= 8976807)[1:100]),names(x),
           names(which(chr20@gtdata@map <= 8970492)[1181:1281]))



#SPO11; chr20		2696827	2709653
y <- chr20[,between(chr20@gtdata@map, 2696827, 2709653)]
x<-y@gtdata@map
SPO11 <- c(names(which(chr20@gtdata@map >= 2709653)[1:100]),names(x),
           names(which(chr20@gtdata@map <= 2696827)[265:365]))




#HORMAD; scaffold00403		114986	125779




#MRE11; chr1		36009953	36027563
y <- chr1[,between(chr1@gtdata@map, 36009953, 36027563)]
x<-y@gtdata@map
MRE11 <- c(names(which(chr1@gtdata@map >= 36027563)[1:100]),names(x),
           names(which(chr1@gtdata@map <= 36009953)[7623:7723]))



#RAD50; chr13		7077388	   7095176
y <- chr13[,between(chr13@gtdata@map, 7077388, 7095176)]
x<-y@gtdata@map
RAD50 <- c(names(which(chr13@gtdata@map >= 7095176)[1:100]),names(x),
           names(which(chr13@gtdata@map <= 7077388)[1060:1160]))



#BRCC3; chr8		34196712	34203152
y <- chr8[,between(chr8@gtdata@map, 34196712, 34203152)]
x<-y@gtdata@map
BRCC3 <- c(names(which(chr8@gtdata@map >= 34203152)[1:100]),
           names(which(chr8@gtdata@map <= 34196712)[5948:6048]))



#MEIOB; chr14		1531265	1537803
y <- chr14[,between(chr14@gtdata@map, 1531265, 1537803)]
x<-y@gtdata@map
MEIOB <- c(names(which(chr14@gtdata@map >= 1537803)[1:100]),
           names(which(chr14@gtdata@map <= 1531265)[1:92]))



#MCMDC2; chr2		32993923	33003490
y <- chr2[,between(chr2@gtdata@map, 32993923, 33003490)]
x<-y@gtdata@map
MCMDC2 <- c(names(which(chr2@gtdata@map >= 33003490)[1:100]),names(x),
          names(which(chr2@gtdata@map <= 32993923)[7557:7657]))



#DMC1; chr1A		47520810	47530717
y <- chr29[,between(chr29@gtdata@map, 47520810, 47530717)]
x<-y@gtdata@map
DMC1 <- c(names(which(chr29@gtdata@map >= 47530717)[1:100]),names(x),
          names(which(chr29@gtdata@map <= 47520810)[6940:7040]))



#RAD51; chr5		8927055	8932713
y <- chr5[,between(chr5@gtdata@map, 8927055, 8932713)]
x<-y@gtdata@map
RAD51 <- c(names(which(chr5@gtdata@map >= 8932713)[1:100]),
          names(which(chr5@gtdata@map <= 8927055)[1986:2086]))



#SYCP1; scaffold00321		59420	68891




#RAD21L – Missing (RAD21 below) chr20		7095064	7096313
y <- chr20[,between(chr20@gtdata@map, 7095064, 7096313)]
x<-y@gtdata@map
RAD21L <- c(names(which(chr20@gtdata@map >= 7096313)[1:100]),names(x),
           names(which(chr20@gtdata@map <= 7095064)[1003:1103]))



#MSH4; chr6		5736	28427
y <- chr6[,between(chr6@gtdata@map, 5736, 28427)]
x<-y@gtdata@map
MSH4 <- c(names(which(chr6@gtdata@map >= 28427)[1:100]),names(x),
            names(which(chr6@gtdata@map <= 5736)[1]))


#MSH5; chr3		58413881	58416316
y <- chr3[,between(chr3@gtdata@map, 58413881, 58416316)]
x<-y@gtdata@map
MSH5 <- c(names(which(chr3@gtdata@map >= 58416316)[1:100]),
          names(which(chr3@gtdata@map <= 58413881)[13515:13615]))



#MLH3; chr5		38276074	38284767
y <- chr5[,between(chr5@gtdata@map, 38276074, 38284767)]
x<-y@gtdata@map
MLH3 <- c(names(which(chr5@gtdata@map >= 38284767)[1:100]),names(x),
          names(which(chr5@gtdata@map <= 38276074)[8228:8328]))



#CNTD1; chr27		2525307	2530033
y <- chr27[,between(chr27@gtdata@map, 2525307, 2530033)]
x<-y@gtdata@map
CNTD1 <- c(names(which(chr27@gtdata@map >= 2530033)[1:100]),
          names(which(chr27@gtdata@map <= 2525307)[52:152]))



#MLH1; chr2		100972180	100990256
y <- chr2[,between(chr2@gtdata@map, 100972180, 100990256)]
x<-y@gtdata@map
MLH1 <- c(names(which(chr2@gtdata@map >= 100990256)[1:100]),names(x),
          names(which(chr2@gtdata@map <= 100972180)[24012:24112]))


#RNF212 chicken blast region 63,770,883 to 63771612
y <- chr4[,between(chr4@gtdata@map, 63770883, 63771612)]
x<-y@gtdata@map
RNF212a <- c(names(which(chr4@gtdata@map >= 63771612)[1:100]),
          names(which(chr4@gtdata@map <= 63770883)[15301:15401]))

#63769455 to 63770249
y <- chr4[,between(chr4@gtdata@map, 63769455, 63770249)]
x<-y@gtdata@map
RNF212b <- c(names(which(chr4@gtdata@map >= 63770249)[1:100]),
             names(which(chr4@gtdata@map <= 63769455)[15301:15401]))



KnownRecombGene_SNPs <- c(ANKH, BRCC3, CBLL1, CCNA1, CNTD1, DMC1, EPGN, EREG,
                          FGFRL1, GNAI2, IGF2, IMPG2, LECT1, MAEA,
                          MCMDC2, MEIOB, MLH1, MLH3, MRE11, MRPL35,
                          MSH4, MSH5, MYC, NFYB, PRIM2, RAD21L, RAD50,
                          RAD51, RAG1, RAG2, RECQL, RGCC, RPL24, RPS17,
                          SPO11, TDG, WBP4, YWHAB, RNF212a, RNF212b)
RecombGene_SNPs_list <- list(ANKH, BRCC3, CBLL1, CCNA1, CNTD1, DMC1, EPGN, EREG,
                          FGFRL1, GNAI2, IGF2, IMPG2, LECT1, MAEA,
                          MCMDC2, MEIOB, MLH1, MLH3, MRE11, MRPL35,
                          MSH4, MSH5, MYC, NFYB, PRIM2, RAD21L, RAD50,
                          RAD51, RAG1, RAG2, RECQL, RGCC, RPL24, RPS17,
                          SPO11, TDG, WBP4, YWHAB, RNF212a, RNF212b)
names(RecombGene_SNPs_list) <- c("ANKH", "BRCC3", "CBLL1", "CCNA1", "CNTD1",
                              "DMC1", "EPGN", "EREG",
                             "FGFRL1", "GNAI2", "IGF2", "IMPG2", "LECT1",
                             "MAEA", "MCMDC2", "MEIOB", "MLH1", "MLH3", "MRE11",
                             "MRPL35", "MSH4", "MSH5", "MYC", "NFYB", "PRIM2",
                             "RAD21L", "RAD50", "RAD51", "RAG1", "RAG2", "RECQL",
                             "RGCC", "RPL24", "RPS17", "SPO11", "TDG", "WBP4", 
                             "YWHAB", "RNF212a", "RNF212b")

RecombGene_BP <- c("ANKH", "BRCC3", "CBLL1", "CCNA1", "CNTD1",
                   "DMC1", "EPGN", "EREG",
                   "FGFRL1", "GNAI2", "IGF2", "IMPG2", "LECT1",
                   "MAEA", "MCMDC2", "MEIOB", "MLH1", "MLH3", "MRE11",
                   "MRPL35", "MSH4", "MSH5", "MYC", "NFYB", "PRIM2",
                   "RAD21L", "RAD50", "RAD51", "RAG1", "RAG2", "RECQL",
                   "RGCC", "RPL24", "RPS17", "SPO11", "TDG", "WBP4", 
                   "YWHAB", "RNF212a", "RNF212b")
RecombGene_BP <- as.data.frame(RecombGene_BP)
names(RecombGene_BP) <- "Name"


RecombGene_BP$Chr <- c("chr 2", "chr 8", "chr 29", "chr 1", "chr 27", "chr 29", "chr 4", "chr 4", "chr 4", "chr 12", "chr 5", "chr 1",
                       "chr 1", "chr 4", "chr 2", "chr 14", "chr 2", "chr 5", "chr 1", "chr 4", "chr 6", "chr 3", "chr 2", "chr 29",
                       "chr 3", "chr 20", "chr 13", "chr 5", "chr 5", "chr 5", "chr 29", "chr 1", "chr 1", "chr 10", "chr 20", "chr 29",
                       "chr 1", "chr 20", "chr 4", "chr 4")
                 
RecombGene_BP$BP_Start <- c(64612378, 34196712, 12380129, 48654156, 2525307, 47520810, 2633114,
                            2606133, 63927199, 5554811, 14489308, 16165948, 54851438, 63332471,
                            32993923, 1531265, 100972180, 38276074, 36009953, 64374204, 5736,
                            58413881, 8901363, 51321757, 87724070, 7095064, 7077388, 8927055,
                            18301496, 18312675, 69730714, 54599781, 16311987, 3155373, 2696827,
                            51372033, 54711630, 8970492, 63770883, 63769455)

RecombGene_BP$BP_Stop <- c(64617461, 34203152, 12383500, 48659220, 2530033, 47530717, 2639966,
                           2611444, 63933410, 5560682, 14490217, 16217147, 54864276, 63333787,
                           33003490, 1537803, 100990256, 38284767, 36027563, 64377025, 28427,
                           58416316, 8903766, 51329193, 87749879, 7096313, 7095176, 8932713,
                           18304618, 18314282, 69742374, 54612288, 16315717, 3157987, 2709653,
                           51378897, 54724061, 8976807, 63771612, 63770249)


KnownRecombGene_SNPs <- unique(KnownRecombGene_SNPs)
KnownRecombGene_SNPs <- KnownRecombGene_SNPs[-873]
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/Dblxoversremoved/")
write.table(KnownRecombGene_SNPs, file = "KnownRecombGene_SNPs.txt", row.names = FALSE, col.names = TRUE)
write.table(RecombGene_BP, file = "RecombGene_Positions.txt", row.names = FALSE, col.names = TRUE)
save(RecombGene_SNPs_list, file = "RecombGene_SNPs_List.RData")

load("RecombGene_SNPs_List.RData")



#Adding MAF for snps; 873


vec <- which(sparrowgen@gtdata@snpnames %in% KnownRecombGene_SNPs)
vec1 <- sparrowgen@gtdata@idnames[vec]
a <- summary(sparrowgen@gtdata[vec1,1:3])
names(a)
afr <- a[, "Q.2"]
maf <- pmin(afr, (1 - afr))
summary(maf)
x <- maf[which(sparrowgen@gtdata@snpnames %in% KnownRecombGene_SNPs)]
KnownRecombGene_SNPs <- as.data.frame(KnownRecombGene_SNPs)
names(KnownRecombGene_SNPs) <- "SNP Name"
KnownRecombGene_SNPs$MAF <- x
