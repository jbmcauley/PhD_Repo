setwd("C:/Users/s1945757/Downloads/Sparrow/scripts")
lep_ped <- read.table("lepmap_pedigree_readable.txt", header = TRUE)


names.all <- as.data.frame(row.names(all))

colnames(names.all) <- "OG"
names.all$names.all <- gsub("^.*\\.","", names.all$OG)
names.all %>% 
  group_by(names.all) %>% 
  filter(row_number() <= .5 * n())

names.all.Paternal <- names.all %>% 
  group_by(names.all) %>% 
  filter(row_number() <= .5 * n())

names.all.Paternal <- as.data.frame(names.all.Paternal)
names.all$OG <- as.character(names.all$OG)
names.all.Paternal$OG <- as.character(names.all.Paternal$OG)
names.all.Maternal <- names.all[which(!names.all$OG %in% names.all.Paternal$OG),]

all$gamete <- row.names(all)
maternal <- all[which(!all$gamete %in% names.all.Paternal$OG),]
paternal <- all[which(all$gamete %in% names.all.Paternal$OG),]
paternal <- paternal[1:28]
maternal <- maternal[1:28]
maternal$sum<- rowSums(maternal[1:28])
paternal$sum<- rowSums(paternal[1:28])
#maternal$gamete <- NULL
#paternal$gamete <- NULL
rm(names.all.Maternal)
rm(names.all.Paternal)

hist(maternal$sum, breaks = 50)
hist(paternal$sum, breaks = 50)
mean(paternal$sum)
mean(maternal$sum)


paternal$family <- gsub("^.*\\.","", row.names(paternal))
paternal$family <- as.numeric(paternal$family) -3
offspring.vector <- paternal$family

maternal$family <- gsub("^.*\\.","", row.names(maternal))
maternal$family <- as.numeric(maternal$family) -3



x  <- as.data.frame(table(offspring.vector))
x <- x$Freq
lep_ped$family <- as.character(lep_ped$family)
#lep_ped$family <- as.numeric(lep_ped$family)


lep_ped$ID_by_Group <- ave(lep_ped$id, lep_ped$family, FUN =  seq_along)

lep_ped[lep_ped$ID_by_Group %in% seq(, by = 1), ]



lep_ped %>%
  group_by(family) %>%
  filter(row_number()==n()-(1:5) | row_number()==n())







# loop through families and keep only fathers
ped <- lep_ped
fathers <- data.frame(family = character(),
                      id = numeric(),
                      father = numeric(),
                      mother = numeric(),
                      sex = numeric(),
                      phenotype = numeric())


for (i in unique(ped$family)) {
  fam <- i 
  df <- ped[c(ped$family==i),]
  df2 <- df[6,]
  fathers <- rbind(fathers,df2)
}




# loop through families and keep only mothers

mothers <- data.frame(family = character(),
                      id = numeric(),
                      father = numeric(),
                      mother = numeric(),
                      sex = numeric(),
                      phenotype = numeric())



for (i in unique(ped$family)) {
  fam <- i 
  df <- ped[c(ped$family==i),]
  df2 <- df[5,]
  mothers <- rbind(mothers,df2)
}


paternal$family <- paste("family", paternal$family, sep="_")
maternal$family <- paste("family", maternal$family, sep="_")

pat_comb <- merge(paternal, fathers, by= "family")
mat_comb <- merge(maternal, mothers, by= "family")


pat_comb$id <- as.numeric(as.character(pat_comb$id))
pat_comb <- pat_comb[complete.cases(pat_comb),]
mat_comb$id <- as.numeric(as.character(mat_comb$id))
mat_comb <- mat_comb[complete.cases(mat_comb),]



pat_pheno_with0 <- aggregate(pat_comb[, "sum"], list(pat_comb$id), mean)
mat_pheno_with0 <- aggregate(mat_comb[, "sum"], list(mat_comb$id), mean)

pat_pheno_with0 <- pat_pheno_with0[1:683,]
mat_pheno_with0 <- mat_pheno_with0[1:727,]
#pat_pheno_NO_0s <- aggregate(pat_comb[, "sum"], list(pat_comb$father), mean)

pat_comb_NO.0s <- pat_comb[pat_comb$sum != 0,]
mat_comb_NO.0s <- mat_comb[mat_comb$sum != 0,]

#No high value removal before mean
pat_pheno_NO.0 <- aggregate(pat_comb_NO.0s[, "sum"], list(pat_comb_NO.0s$id), mean)
mat_pheno_NO.0 <- aggregate(mat_comb_NO.0s[, "sum"], list(mat_comb_NO.0s$id), mean)
names(pat_pheno_NO.0) <- c("RRID", "RecombCount")
names(mat_pheno_NO.0) <- c("RRID", "RecombCount")
pat_pheno_NO.0$RRID <- as.character(pat_pheno_NO.0$RRID)
mat_pheno_NO.0$RRID <- as.character(mat_pheno_NO.0$RRID)
Lepmap_Pheno_No.0 <- rbind(pat_pheno_NO.0, mat_pheno_NO.0)

#Removing high values before taking mean
x <- pat_comb_NO.0s[which(pat_comb_NO.0s$sum <= (mean(pat_comb_NO.0s$sum)+2*sd(pat_comb_NO.0s$sum))),] #remove values less than 95% of 2sd
y <- mat_comb_NO.0s[which(mat_comb_NO.0s$sum <= (mean(mat_comb_NO.0s$sum)+2*sd(mat_comb_NO.0s$sum))),] #remove values less than 95% of 2sd
pat_pheno_NO.0 <- aggregate(x[, "sum"], list(x$id), mean)
mat_pheno_NO.0 <- aggregate(y[, "sum"], list(y$id), mean)

both_pheno_NO.0 <- rbind(pat_pheno_NO.0, mat_pheno_NO.0)
names(both_pheno_NO.0) <- c("id","MeanACC")
both_pheno_NO.0$id <- as.character(both_pheno_NO.0$id)
LP_No0s_2sdRemoved <- rbind(x,y)
LP_No0s_2sdRemoved <- LP_No0s_2sdRemoved[c(1,30:34)]




###
###
###
#Animal Models
LP_No0s_2sdRemoved$id <- as.character(LP_No0s_2sdRemoved$id)
names(LP_No0s_2sdRemoved)[2] <- "ACC"
LP_No0s_2sdRemoved <- left_join(LP_No0s_2sdRemoved, idkey, by = "id")
LP_No0s_2sdRemoved$sex <- as.factor(as.character(LP_No0s_2sdRemoved$sex)) 
LP_No0s_2sdRemoved$id <- as.factor(LP_No0s_2sdRemoved$id)
LP_No0s_2sdRemoved$ringnr <- as.factor(LP_No0s_2sdRemoved$ringnr)

grm.auto <- read.table("workfiles.GRM.adj.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
ids.auto <- read.table("workfiles.GRM.adj.grm.id")

source("makeGRM.R")
grminv <- makeGRM(grm.auto, ids.auto, id.vector = LP_No0s_2sdRemoved$ringnr) # vector of IDs from the datasset that you use for the asreml model

# You MUST specify this this is an inverted matrix!
attr(grminv, which = "INVERSE") <- TRUE



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. RUN ANIMAL MODELS                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Quick data explore

ggplot(LP_No0s_2sdRemoved, aes(sex, ACC)) +
  geom_boxplot()

ggplot(LP_No0s_2sdRemoved, aes(age, ACC, col = sex)) +
  geom_point() +
  stat_smooth(method = "lm")

# Pedigree vs GRM

model.ped <- asreml(fixed= ACC ~ 1 + sex, 
                    random= ~ vm(id,ainv) + ide(id),
                    residual= ~ idv(units), 
                    na.action = na.method(x = "omit", y = "omit"),
                    data= LP_No0s_2sdRemoved)

summary.asreml(model.ped, coef = T)$coef.fixed # Fixed effects
asreml4pin(model.ped)

model.grm <- asreml(fixed = ACC ~ 1 + sex,
                    random = ~ vm(ringnr, grminv) + ide(ringnr),
                    residual= ~ idv(units),
                    na.action = na.method(x = "omit", y = "omit"),
                    data = LP_No0s_2sdRemoved)
summary.asreml(model.grm, coef = T)$coef.fixed # Fixed effects
asreml4pin(model.grm)



model.ped.F <- asreml(fixed= ACC ~ 1 + sex, 
                    random= ~ vm(id,ainv) + ide(id),
                    residual= ~ idv(units), 
                    na.action = na.method(x = "omit", y = "omit"),
                    data= LP_No0s_2sdRemoved)

summary.asreml(model.ped, coef = T)$coef.fixed # Fixed effects
asreml4pin(model.ped)

model.grm.F <- asreml(fixed = ACC ~ 1 + sex,
                    random = ~ vm(ringnr, grminv) + ide(ringnr),
                    residual= ~ idv(units),
                    na.action = na.method(x = "omit", y = "omit"),
                    data = LP_No0s_2sdRemoved)
summary.asreml(model.grm, coef = T)$coef.fixed # Fixed effects
asreml4pin(model.grm)

