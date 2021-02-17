#Load in libraries and data
library(GenABEL)
library(crimaptools)
library(dplyr)
library(asreml)
library(forcats)
library(RepeatABEL)

setwd("~/PhD/PhD_Repo/data")
load("sparrowABEL.RData")
load("sparrowgen.RData")
load("sparrow.xovers.clean.all.RData")
Rcount_RRID <- all.xovers %>% group_by(RRID, Family) %>% summarise(ACC = sum(RecombCount))
sparrow.abel <- add.phdata(sparrow.abel, Rcount_RRID, ACC)
setwd("~/PhD/PhD_Repo/data")
names(Rcount_RRID)[1] <- "id"
load("sparrow.gkin.RData")

#Mean_ACC <- aggregate(ACC ~ id, Rcount_RRID, mean)
#sparrow.abel <- add.phdata(sparrow.abel, Mean_ACC)


sparrow.gkin <- ibs(sparrowgen, weight = "freq")
save(sparrow.gkin, file = "sparrow.gkin.RData")
sparrow.dist <- as.dist(0.5 - sparrow.gkin)
#sparrow.mds <- cmdscale(sparrow.dist)
#save(sparrow.mds, file="sparrow.mds.RData")
load("sparrow.mds.RData")

#Mean_ACC <- aggregate(ACC ~ id, Rcount_RRID, mean)
#sparrow.abel <- add.phdata(sparrow.abel, Mean_ACC)

km <- kmeans(sparrow.mds, centers = 2, nstart =1000)
cl1 <- names(which(km$cluster == 1))
cl2 <- names(which(km$cluster == 2))

qt <- qtscore(ACC,sparrow.abel)
qt.kin <- egscore(ACC,sparrow.abel, kin = sparrow.gkin)
qts.e <- qtscore(ACC ~ sex,data =sparrow.abel, times = 200)
h2 <- polygenic(ACC, kin = sparrow.gkin, data = sparrow.abel)
mms <- mmscore(h2, data = sparrow.abel)
grs <- qtscore(h2$pgres, data = sparrow.abel, clambda = FALSE)
grs.e <- qtscore(h2$pgres, data = sparrow.abel, times = 200, clam = FALSE)

plot(an0)



#RepeatABEL
#Redue with original GenABEL object

source("rGLSadj.R")


sparrow.gkin.sym <- sparrow.gkin
sparrow.gkin.sym[upper.tri(sparrow.gkin.sym)] = t(sparrow.gkin.sym)[upper.tri(sparrow.gkin.sym)]
sparrow.gkin.sym <- sparrow.gkin.sym * 2

Rcount_RRID$id <- as.character(Rcount_RRID$id)
Rcount_RRID <- merge(Rcount_RRID, phdata, by = "id")


system.time(gwasACC_prefit2 <- preFitModel(fixed = ACC ~ sex, 
                                            id.name = "id",
                                            genabel.data = sparrowgen,
                                            phenotype.data = Rcount_RRID, 
                                            corStruc = list(id = list("GRM")),
                                            GRM = sparrow.gkin.sym))




system.time(gwasACC2 <- rGLSadj(ACC ~ sex,
                                 genabel.data = sparrowgen,
                                 phenotype.data = Rcount_RRID, 
                                 id = "id",
                                 V = gwasACC_prefit2$V,
                                 GRM=sparrow.gkin.sym))

# The process_rGLSadj_results function will do the lambda correction for
# population structure and also return the expected distribution of P-values.

gwasACCres2 <- process_rGLSadj_results(gwasACC2, sparrowgen)

head(gwasACCres)

plot_rGLSadj_results(gwasACCres2)
plot_rGLSadj_results(gwasACCres2, PP = TRUE)




#gtca alternative option for GWAS?


