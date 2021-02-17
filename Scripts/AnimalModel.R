#Load in libraries and data
library(GenABEL)
library(crimaptools)
library(dplyr)
library(asreml)
library(forcats)

setwd("~/PhD/PhD_Repo/data")
source("ASReml4.EstEffects.R")
load("sparrowABEL.RData")

load("sparrow.xovers.clean.all.RData")
Rcount_RRID <- all.xovers %>% group_by(RRID, Family) %>% summarise(ACC = sum(RecombCount))
names(Rcount_RRID) <- c("id", "ACC")
Rcount_RRID$id <- as.factor(Rcount_RRID$id)

setwd("~/PhD/PhD_Repo/data")


#Set up pedigree using ID-Recode-Key
idkey <- read.table("ID-Recode-Key.txt", header = T, stringsAsFactors = F)[,c(2, 4)]
names(idkey) <- c("ringnr", "id")
ped <- read.table("SNP_pedigree_Helgeland_05122017.txt", header = TRUE, stringsAsFactors = F)

# Susie sanity check - add the dummy IDs to the key (just give them the same ID before and after)

ped$id[which(!ped$id %in% idkey$ringnr)]    # These are the ids in pedigree that are not the key.

x <- data.frame(id = ped$id[which(!ped$id %in% idkey$ringnr)])
x$ringnr <- x$id

idkey <- rbind(idkey, x)
rm(x)

level_vec <- idkey$ringnr
replacement_vec <- idkey$id
ped[] <- lapply(ped, function(x) forcats::lvls_revalue(factor(x, levels = level_vec), replacement_vec))
names(ped) <- c("id", "dam", "sire")
ped <- ped[, c(1,3,2)]
ped <- as.data.frame(ped)




#Create inverse relatedness matrix
ainv <- ainverse(ped[,1:3])

#Data for model
sparrow <- sparrow.abel@phdata
sparrow$id <- as.factor(sparrow$id)
sparrow <- merge(sparrow, Rcount_RRID, by = "id", sort = TRUE)
sparrow$id <- as.factor(sparrow$id)
sparrow$sex[sparrow$sex == 1] <- 2
sparrow$sex[sparrow$sex == 0] <- 1
sparrow$sex <- as.factor(sparrow$sex)


#Run Model
sparrow.asr <- asreml(fixed= ACC ~ 1 + sex, 
                      random= ~ vm(id,ainv),
                      residual= ~ idv(units), 
                      data= sparrow)

vpredict(sparrow.asr, h2 ~  V1/(V1+V2)) # heritability


#GRM
#Make update-ids files 
idkey$OldFamID <- 1
idkey$OldID <- idkey$ringnr
idkey$NewFamID <- 1
idkey$NewID <- idkey$id
idkey$ringnr <- NULL
idkey$id <-NULL
write.table(idkey, file = "updateTDs.txt", col.names = FALSE, row.names = FALSE)

#Make Binary Files
system("plink --file Pdo_200k_n3960_21032017 --autosome-num 31 --maf 0.05 --make-bed --out workfiles")

#Update ids
system("plink --bfile workfiles --update-ids updateIDs.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles")

system("gcta64.exe --bfile workfiles --autosome --autosome-num 31 --make-grm-gz --out workfiles.GRM")
system("gcta64.exe --grm-gz workfiles.GRM --grm-adj 0 --make-grm-gz --out workfiles.GRM.adj")

grm.auto <- read.table("workfiles.GRM.adj.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.auto <- read.table("workfiles.GRM.adj.grm.id")

Rcount_RRID <- as.data.frame(Rcount_RRID)
names(Rcount_RRID) <- c("id","FAMID","ACC")
idvec <- Rcount_RRID %>%
  left_join(idkey, by = "id")

source("makeGRM.R")
grminv <- makeGRM(grm.auto, ids.auto, id.vector = idvec$ringnr) # vector of IDs from the datasset that you use for the asreml model

# You MUST specify this this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE

#~~ Have to make everything that is character in your model into a factor

sparrow <- sparrow.abel@phdata
sparrow$id <- as.factor(sparrow$id)
sparrow <- merge(sparrow, idvec, by = "id", sort = TRUE)
sparrow$ringnr <- as.factor(sparrow$ringnr)
sparrow$sex[sparrow$sex == 1] <- 2
sparrow$sex[sparrow$sex == 0] <- 1
sparrow$sex <- as.factor(sparrow$sex)

#~~ Data can only have pedigreed individuals

sparrow <- subset(sparrow, id %in% ped$id)
sparrow <- droplevels(sparrow) # need to run this command

#~~ Run the model. vm(ID, ainv) is the relatedness matricx, ide(ID) is the
#   remaining variance due to the individual (repeated measures). In your models,
#   you don't need to run ide(ID). You can use grminv if you want to look at GRM.

model1 <- asreml(fixed= ACC ~ 1 + sex, 
                 random= ~ vm(id,ainv) + ide(id),
                 residual= ~ idv(units), 
                 data= sparrow)

model2 <- asreml(fixed = ACC ~ 1 + sex,
                 random = ~ vm(ringnr, grminv) + ide(ringnr),
                 residual= ~ idv(units),
                 data = sparrow)

#~~ Look at the fixed and random effects

summary.asreml(model1, coef = T)$coef.fixed # Fixed effects
summary.asreml(model1, coef = T)$varcomp    # random effects
vpredict(model1, hA ~ V1/(V1+V2+V3))

summary.asreml(model2, coef = T)$coef.fixed # Fixed effects
summary.asreml(model2, coef = T)$varcomp    # random effects
vpredict(model2, hA ~ V1/(V1+V2+V3))







####
####
#### Additional Models
####
####




Rcount_RRID <- as.data.frame(Rcount_RRID)
names(Rcount_RRID) <- c("id","FAMID","ACC")
Rcount_RRID$id <- as.character(Rcount_RRID$id)
ped$id <- as.character(ped$id)
test <- Rcount_RRID %>%
  left_join(ped, by = "id")
test <- test %>%
  left_join(idkey, by = "id")
sparrow <- sparrow.abel@phdata
sparrow$id <- as.factor(sparrow$id)
test <- merge(sparrow, test, by = "id", sort = TRUE)
test$sex[test$sex == 1] <- 2
test$sex[test$sex == 0] <- 1


#Adding in additional variables
Morph <- read.csv("AdultMorphology-pre4_SJ.csv")
Morph$ringnr <- as.character(Morph$ringnr)
Morph <- Morph[!duplicated(Morph$ringnr),]
Morph$sex <- NULL
test <- test %>% 
  left_join(Morph[,1:3], by = "ringnr")
names(test)[8:9] <- c("hatchislandID", "hatchyearID")
names(test)[7] <- "ringnrID"
test$idGAM <- read.table(text = as.character(test$FAMID), sep = "_")$V2
test$idGAM <- as.character(test$idGAM)
test <- test %>% 
  left_join(idkey, by = c("idGAM"="id"))
names(test)[11] <-  "ringnrGAM"
test <- test %>% 
  left_join(Morph[,1:3], by = c("ringnrGAM"="ringnr") )
names(test)[12:13] <- c("hatchislandGAM", "hatchyearGAM")
test <- test %>% 
  left_join(idkey, by = c("dam"="id"))
names(test)[14] <-  "ringnrDAM"

#Make all vars factors
test$id <- as.factor(test$id)
test$idGAM <- as.factor(test$idGAM)
test$dam <- as.factor(test$dam)
test$ringnrID <- as.factor(test$ringnrID)
test$ringnrGAM <- as.factor(test$ringnrGAM)
test$ringnrDAM <- as.factor(test$ringnrDAM)
test$sex <- as.factor(test$sex)
test$hatchislandID <- as.factor(test$hatchislandID)
test$hatchyearID <- as.factor(test$hatchyearID)
test$hatchislandGAM <- as.factor(test$hatchislandGAM)
test$hatchyearGAM <- as.factor(test$hatchyearGAM)
test.female <- test %>% 
  filter(sex == 1)
test.male <- test %>% 
  filter(sex ==2)


#Basic Model with maternal effects
model.ped.mat <- asreml(fixed= ACC ~ 1 + sex, 
                        random= ~ vm(id,ainv) + ide(id) + vm(dam,ainv) + ide(dam),
                        residual= ~ idv(units), 
                        na.action = na.method(x = "omit", y = "omit"),
                        data= test)

model.grm.mat <- asreml(fixed = ACC ~ 1 + sex,
                        random = ~ vm(ringnrID, grminv) + ide(ringnrID) + vm(ringnrDAM, grminv) + ide(ringnrDAM),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = test)
asreml4pin(model.ped.mat) 
asreml4pin(model.grm.mat)
#All Variables in Model: hatch island, hatch year for FID and Gamete

model.ped.all <- asreml(fixed= ACC ~ 1 + sex, 
                 random= ~ vm(id,ainv) + ide(id) + vm(dam,ainv) + ide(dam) + hatchyearID + hatchislandID + hatchislandGAM + hatchyearGAM,
                 residual= ~ idv(units),
                 na.action = na.method(x = "omit", y = "omit"),
                 data= test)
asreml4pin(model.ped.all) #Ha = 4.769 * 10^-7





#Female only Models

#Simple Models
model.ped.F <- asreml(fixed= ACC ~ 1, 
                        random= ~ vm(id,ainv) + ide(id) + vm(dam,ainv) + ide(dam),
                        residual= ~ idv(units), 
                        na.action = na.method(x = "omit", y = "omit"),
                        data= test.female)

model.grm.F <- asreml(fixed = ACC ~ 1,  #WARNING MESSAGE ringnrDAM has levels in dat that are missing in grminv!!! Need to take a look
                        random = ~ vm(ringnrID, grminv) + ide(ringnrID) + vm(ringnrDAM, grminv) + ide(ringnrDAM),
                        residual= ~ idv(units),
                        na.action = na.method(x = "omit", y = "omit"),
                        data = test.female)

asreml4pin(model.ped.F)  # Ha = 7.6 *10^-8
asreml4pin(model.grm.F)  # Ha = 0.209

#All variables inlcuded
model.ped.F.all <- asreml(fixed= ACC ~ 1, 
                        random= ~ vm(id,ainv) + ide(id) + vm(dam,ainv) + ide(dam) + hatchyearID + hatchislandID + hatchislandGAM + hatchyearGAM,
                        residual= ~ idv(units), 
                        na.action = na.method(x = "omit", y = "omit"),
                        data= test.female)
asreml4pin(model.ped.F.all)




#Male only Models

#Simple Models
model.ped.M <- asreml(fixed= ACC ~ 1, 
                      random= ~ vm(id,ainv) + ide(id) + vm(dam,ainv) + ide(dam),
                      residual= ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data= test.male)

model.grm.M <- asreml(fixed = ACC ~ 1,
                      random = ~ vm(ringnrID, grminv) + ide(ringnrID) + vm(ringnrDAM, grminv) + ide(ringnrDAM),
                      residual= ~ idv(units),
                      na.action = na.method(x = "omit", y = "omit"),
                      data = test.male)

asreml4pin(model.ped.M) # Ha = 0.081
asreml4pin(model.grm.M) # Ha = 0.084

#All variables inlcuded
model.ped.M.all <- asreml(fixed= ACC ~ 1, 
                          random= ~ vm(id,ainv) + ide(id) + vm(dam,ainv) + ide(dam) + hatchyearID + hatchislandID + hatchislandGAM + hatchyearGAM,
                          residual= ~ idv(units),
                          na.action = na.method(x = "omit", y = "omit"),
                          data= test.male)
asreml4pin(model.ped.M.all) #8.723 * 10^-7


