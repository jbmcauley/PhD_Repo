#Inverse Regression
install.packages("MCMCGlmm")
library(lme4)
library(GenABEL)


data1 <- sparrow.abel

#rm(sparrowgen.Helgeland)

#-------------------LEFT OFF HERE 9/12/19 

#Getting the Phenotype in the genABEL object, using mean values as there are multiple measures.
Phenotypes <- read.csv("AdultMorphology-pre4_SJ.csv")
Phenotypes$island<- as.factor(as.character(Phenotypes$island))
Phenotypes.means <- aggregate(Phenotypes[,c(8,10,12:16)],list(Phenotypes$ringnr), mean)
Phenotypes.means$age <- as.factor(Phenotypes.means$age)
names(Phenotypes.means)[1] <- "id"
age <- Phenotypes.means[,c(1,3)]
rm(Phenotypes)
data1 <- add.phdata(data1, Phenotypes.means)
rm(Phenotypes.means)
rm(age)


#Basic QC with genABEL
qc1 <- check.marker(data1, callrate = .95, ibs.threshold = .9, maf = .0001,p.level = 0)


#NEED TO SORT OUT HOW TO PROPERLY CHECK ZZ WZ
data1 <- data1[qc1$idok, qc1$snpok]
data1 <- Xfix(data1)
save(data1, file = "data1_AllSparrows.RData")
load("data1_AllSparrows.RData")
descriptives.marker(data1)[2]

#Construct matrix of genomic kinship
data1.gkin <- ibs(data1[,autosomal(data1)], weight = "freq")
save(data1.gkin, file = "gkin_AllSparrows.RData")
load("gkin_AllSparrows.RData")
#Transform into distance matrix with...
data1.dist <- as.dist(0.5-data1.gkin)
#Perform Classical Multidimensional Scaling..
data1.mds <- cmdscale(data1.dist)
plot(data1.mds)
#Are outliers evident?-------------------------------------------

km <- kmeans(data1.mds, centers = 2, nstart = 1000)
cl1 <- names(which(km$cluster ==1))
cl2 <- names(which(km$cluster==2))
data2 <- data1[cl1,]

qc2 <- check.marker(data2, hweids=(phdata(data2)$billD), fdr=0.2)
#Large portion of markers lost at this step!?



data1.qt <- qtscore(billD, data1)
lambda(data1.qt)


pop <- as.numeric(idnames(data1) %in% cl1)
data1.sa <- qtscore(billD, data= data1, strata = pop)
lambda(data1.sa)

#Adjusting for genetic substructure with Price et al. methods, using PC of gk matrix to adjust both
#phenotypes and genotypes for possible strat.
data1.eg <- egscore(billD~sex+age+age^2,data=data1,kin=data1.gkin)
lambda(data1.eg)
plot(data1.eg)
#Adjustment for strat w/ principal components of genetic variation. 
pcs <- cmdscale(data1.dist, k = 10)
pcs[1:5,]
plot(pcs)
data1.pca <- qtscore(billD~pcs[,1]+pcs[,2]+pcs[,3]+age+age^2+sex,data1)

#Use of full genomic kinship matrix 
h2a <- polygenic(billD~age+age^2+sex,data=data1,kin=data1.gkin)
h2a$esth2 #-------------This estimate of ~0.42 is higher than the papers result of 0.35

#Mixed model approaches
data1.mm <- mmscore(h2a,data1)
lambda(data1.mm)$est
plot(data1.mm)
save(data1.mm, file = "data1mm_Helgeland_01_2018.RData")

par(mfrow=c(3,1))
plot(data1.eg,ylim=c(1,9))  
plot(data1.pca,ylim=c(1,9)) 
plot(data1.mm,ylim=c(1,9))



#Analysis of Family Data Section

erfs <- data1
pkins <-  matrix(rnorm(nids(erfs)^2,sd=0.01),ncol=nids(erfs),nrow=nids(erfs)) 
par(mfrow=c(1,1))
hist(pkins[lower.tri(pkins)])
#gkins <- ibs(erfs[,autosomal(erfs)],weight="freq") 
pkin_gkin_plot <- plot(pkins[lower.tri(pkins)],data1.gkin[lower.tri(data1.gkin)])
cor(pkins[lower.tri(pkins)],data1.gkin[lower.tri(data1.gkin)])
pop <- 1*(pcs[,1]>0) 
table(pop)



#BillL
data1.eg.billL <- egscore(billL~sex+age+age^2,data=data1,kin=data1.gkin)
data1.pca.billL <- qtscore(billL~pcs[,1]+pcs[,2]+pcs[,3]+age+age^2+sex,data1)
h2a.billL <- polygenic(billL~age+age^2+sex,data=data1,kin=data1.gkin)
data1.mm.billL <- mmscore(h2a.billL,data1)

par(mfrow=c(3,1))
plot(data1.eg.billL,ylim=c(1,9))  
plot(data1.pca.billL,ylim=c(1,9)) 
plot(data1.mm.billL,ylim=c(1,9))





lmer(Vp = , data = )

