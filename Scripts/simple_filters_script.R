
#~~ Load libraries

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(GenABEL)

#~~ load the data and functions

load("3_Examine_SNP_Clusters.RData", verbose = T)




z <- fullsnp[which(fullsnp$stdevTheta_AB < .046),]
z <- z[which(z$stdevTheta_AA < .046),]
z <- z[which(z$stdevTheta_BB < .046),]

y <- 1:(length(z$cust_id)/4)
y.1 <- (length(z$cust_id)/4)+1:((length(z$cust_id)/4)*2)
y.2 <- ((length(z$cust_id)/4)*2)+1:((length(z$cust_id)/4)*3)
y.3 <- ((length(z$cust_id)/4)*3)+1:((length(z$cust_id)/4)*4)


x.1 <- snp_extract_func(c(z$cust_id[y]))
x.2 <- snp_extract_func(c(z$cust_id[y.1]))
x.3 <- snp_extract_func(c(z$cust_id[y.2]))
x.4 <- snp_extract_func(c(z$cust_id[y.3]))

new_x1 <- list()
counter <- 0
for (i in 1:length(unique(x.1$cust_id))) {
  x1 <- x.1[which(x.1$cust_id == unique(x.1$cust_id)[i]),] 
  aa <- x1[which(x1$AB_geno == "AA"),]
  ab <- x1[which(x1$AB_geno == "AB"),]
  bb <- x1[which(x1$AB_geno == "BB"),]
  
  counter <- counter + 1
  new_x1[[counter]] <- aa[between(aa$Theta ,mean(aa$Theta)-2*sd(aa$Theta),mean(aa$Theta)+2*sd(aa$Theta)),]
  
  counter <- counter + 1
  new_x1[[counter]] <- ab[between(ab$Theta ,mean(ab$Theta)-2*sd(ab$Theta),mean(ab$Theta)+2*sd(ab$Theta)),]
  
  counter <- counter + 1
  new_x1[[counter]] <- bb[between(bb$Theta ,mean(bb$Theta)-2*sd(bb$Theta),mean(bb$Theta)+2*sd(bb$Theta)),]
  
} 
new_x2 <- list()
counter <- 0
for (i in 1:length(unique(x.2$cust_id))) {
  x1 <- x.2[which(x.2$cust_id == unique(x.2$cust_id)[i]),] 
  aa <- x1[which(x1$AB_geno == "AA"),]
  ab <- x1[which(x1$AB_geno == "AB"),]
  bb <- x1[which(x1$AB_geno == "BB"),]
  
  counter <- counter + 1
  new_x2[[counter]] <- aa[between(aa$Theta ,mean(aa$Theta)-2*sd(aa$Theta),mean(aa$Theta)+2*sd(aa$Theta)),]
  
  counter <- counter + 1
  new_x2[[counter]] <- ab[between(ab$Theta ,mean(ab$Theta)-2*sd(ab$Theta),mean(ab$Theta)+2*sd(ab$Theta)),]
  
  counter <- counter + 1
  new_x2[[counter]] <- bb[between(bb$Theta ,mean(bb$Theta)-2*sd(bb$Theta),mean(bb$Theta)+2*sd(bb$Theta)),]
  
} 
new_x3 <- list()
counter <- 0
for (i in 1:length(unique(x.3$cust_id))) {
  x1 <- x.3[which(x.3$cust_id == unique(x.3$cust_id)[i]),] 
  aa <- x1[which(x1$AB_geno == "AA"),]
  ab <- x1[which(x1$AB_geno == "AB"),]
  bb <- x1[which(x1$AB_geno == "BB"),]
  
  counter <- counter + 1
  new_x3[[counter]] <- aa[between(aa$Theta ,mean(aa$Theta)-2*sd(aa$Theta),mean(aa$Theta)+2*sd(aa$Theta)),]
  
  counter <- counter + 1
  new_x3[[counter]] <- ab[between(ab$Theta ,mean(ab$Theta)-2*sd(ab$Theta),mean(ab$Theta)+2*sd(ab$Theta)),]
  
  counter <- counter + 1
  new_x3[[counter]] <- bb[between(bb$Theta ,mean(bb$Theta)-2*sd(bb$Theta),mean(bb$Theta)+2*sd(bb$Theta)),]
  
} 
new_x4 <- list()
counter <- 0
for (i in 1:length(unique(x.4$cust_id))) {
  x1 <- x.4[which(x.4$cust_id == unique(x.4$cust_id)[i]),] 
  aa <- x1[which(x1$AB_geno == "AA"),]
  ab <- x1[which(x1$AB_geno == "AB"),]
  bb <- x1[which(x1$AB_geno == "BB"),]
  
  counter <- counter + 1
  new_x4[[counter]] <- aa[between(aa$Theta ,mean(aa$Theta)-2*sd(aa$Theta),mean(aa$Theta)+2*sd(aa$Theta)),]
  
  counter <- counter + 1
  new_x4[[counter]] <- ab[between(ab$Theta ,mean(ab$Theta)-2*sd(ab$Theta),mean(ab$Theta)+2*sd(ab$Theta)),]
  
  counter <- counter + 1
  new_x4[[counter]] <- bb[between(bb$Theta ,mean(bb$Theta)-2*sd(bb$Theta),mean(bb$Theta)+2*sd(bb$Theta)),]
  
} 


final_x1 <- do.call(rbind, new_x1)
final_x2 <- do.call(rbind, new_x2)
final_x3 <- do.call(rbind, new_x3)
final_x4 <- do.call(rbind, new_x4)
