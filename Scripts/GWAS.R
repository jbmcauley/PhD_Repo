rm(list = ls())
#Load in libraries and data
library(GenABEL)
library(crimaptools)
library(dplyr)
library(asreml)
library(forcats)
library(RepeatABEL)
setwd("C:/Users/s1945757/Dropbox/McAuley PhD - Data/Scripts/GWAS")
source("rGLSadj.R")
load("sparrowgen.RData")
load("sparrow.john.RData") #This is the data after sngltons removed, and InfCount < 5 removed. 
load("sparrow.gkin.RData")

#gkin setup
sparrow.gkin.sym <- sparrow.gkin
sparrow.gkin.sym[upper.tri(sparrow.gkin.sym)] = t(sparrow.gkin.sym)[upper.tri(sparrow.gkin.sym)]
sparrow.gkin.sym <- sparrow.gkin.sym * 2



#RepeatABEL on unaltered sparrowgen object
system.time(gwasACC_prefit2 <- preFitModel(fixed = ACC ~ sex, 
                                            id.name = "id",
                                            genabel.data = sparrowgen,
                                            phenotype.data = sparrow, 
                                            corStruc = list(id = list("GRM")),
                                            GRM = sparrow.gkin.sym))
system.time(gwasACC2 <- rGLSadj(ACC ~ sex,
                                 genabel.data = sparrowgen,
                                 phenotype.data = sparrow, 
                                 id = "id",
                                 V = gwasACC_prefit2$V,
                                 GRM=sparrow.gkin.sym))

# The process_rGLSadj_results function will do the lambda correction for
# population structure and also return the expected distribution of P-values.
gwasACCres2 <- process_rGLSadj_results(gwasACC2, sparrowgen)
head(gwasACCres)
plot_rGLSadj_results(gwasACCres2)
plot_rGLSadj_results(gwasACCres2, PP = TRUE)







#Inspect sparrowgen object
qc0snp <- summary.snp.data(sparrowgen@gtdata)
hist(qc0snp$CallRate, breaks = 50)
qc0id <- perid.summary(sparrowgen)
hist(qc0id$CallPP, breaks = 50) 


#Choose QC thresholds for check.marker() function based on histograms
qc1 <- check.marker(sparrowgen, callrate = 0.9, perid.call = 0.8, p.level = 0)
data1 <- sparrowgen[qc1$idok, qc1$snpok]

#RepeatABEL 
system.time(gwasACC_prefit2 <- preFitModel(fixed = ACC ~ sex, 
                                           id.name = "id",
                                           genabel.data = data1,
                                           phenotype.data = sparrow, 
                                           corStruc = list(id = list("GRM")),
                                           GRM = sparrow.gkin.sym))
system.time(gwasACC2 <- rGLSadj(ACC ~ sex,
                                genabel.data = data1,
                                phenotype.data = sparrow, 
                                id = "id",
                                V = gwasACC_prefit2$V,
                                GRM=sparrow.gkin.sym))

# The process_rGLSadj_results function will do the lambda correction for
# population structure and also return the expected distribution of P-values.

gwasACCres2 <- process_rGLSadj_results(gwasACC2, data1)
head(gwasACCres2)
plot_rGLSadj_results(gwasACCres2)