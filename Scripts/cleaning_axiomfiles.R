


y <- floor(runif(20, min=1, max=length(fullsnp$cust_id)))
y <- 1:191082
#rand num selection for snps


x <- snp_extract_func(c(fullsnp$cust_id[y]))

snp_extract_plot(x)                     # with AB geno and R ~ Theta
snp_extract_plot(x, theta = F)          # with AB geno and B ~ A
snp_extract_plot(x, AB = F)             # with 0,1,2 geno and R ~ Theta
snp_extract_plot(x, theta = F, AB = F)  # with 0,1,2 geno and B ~ A


ggplot(data = x, aes(Theta, R, col = factor(AB_geno))) +
  geom_point(alpha = 0.5) +
  facet_wrap(~cust_id, scales = "free_y") +
  theme(legend.position = "top")

#Theta Histogram
ggplot(data = x, aes(Theta, col = factor(AB_geno))) +
  geom_histogram(binwidth = .01) +
  facet_wrap(~cust_id, scales = "free_y") +
  theme(legend.position = "top")

#R Histogram
ggplot(data = x, aes(R, col = factor(AB_geno))) +
  geom_histogram(binwidth = 60) +
  facet_wrap(~cust_id, scales = "free_y") +
  theme(legend.position = "top")



x %>% 
  group_by(cust_id, AB_geno) %>% 
  summarise(ThetaMean = mean(Theta), ThetaSD = sd(Theta), Rmean = mean(R), Rsd = sd(R))



z <- x %>% 
  group_by(cust_id, AB_geno) %>% 
  summarise(ThetaMean = mean(Theta), ThetaSD = sd(Theta), ThetaVar = var(Theta),
            Rmean = mean(R), Rsd = sd(R), Rvar = var(Theta))

ggplot(data = z, aes(ThetaSD)) +
  geom_histogram(bins=50) +
  facet_wrap(~cust_id, scales = "free_y") +
  theme(legend.position = "top")














#Testing various removals

#base changes

t <- fullsnp[which(fullsnp$stdevTheta_AB < .046),]
t <- fullsnp[which(fullsnp$madTheta_AB < .041),]


t <- fullsnp[which(fullsnp$maxTheta_AA[] < fullsnp$minTheta_AB[]),]

t <- fullsnp[which(fullsnp$meanTheta_AB >= .4 & fullsnp$meanTheta_AB <= .6),]
t <- fullsnp[which(fullsnp$meanTheta_AA < .14),]
t <- fullsnp[which(fullsnp$meanTheta_BB > .75),]

#additional restrictions
t <- t[which(t$stdevTheta_AB < .05),]
t <- t[which(t$stdevTheta_AA < .046),]
t <- t[which(t$stdevTheta_BB < .046),]

t <- t[which(t$madTheta_AA < .041),]
t <- t[which(t$madTheta_BB < .041),]

t <- t[which(t$maxTheta_AA[] < t$minTheta_AB[]),]
t <- t[which(t$minTheta_BB[] > t$maxTheta_AB[]),]
t <- t[which(t$meanTheta_AB >= .4 & t$meanTheta_AB <= .6),]
t <- t[which(t$meanTheta_AA < .14),]
t <- t[which(t$meanTheta_BB > .75),]


y <- floor(runif(20, min=1, max=length(t$cust_id)))
y <- 1:27966
y.1 <- 27967:41949
y.2 <- 41950:55932
#rand num selection for snps


x <- snp_extract_func(c(t$cust_id[y]))
x.1 <- snp_extract_func(c(t$cust_id[y.1]))
x.2 <- snp_extract_func(c(t$cust_id[y.2]))

z <- x %>% 
  group_by(cust_id, AB_geno) %>% 
  summarise(ThetaMean = mean(Theta), ThetaSD = sd(Theta), ThetaVar = var(Theta),
            Rmean = mean(R), Rsd = sd(R), Rvar = var(Theta))
snp_extract_plot(x)



z <- x %>% 
  group_by(c("cust_id", "AB_geno"))
  
  subset(x, Theta < 2*sd(Theta), select = Probe.Set.ID:AB_geno)
  
z %>% 
  group_by_(.dots=c("cust_id", "AB_geno")) %>% 
  select_()
  summarise(sd(Theta))

  
  
x1 <- x[which(x$cust_id == unique(x$cust_id)[1]),]  
x1aa <- x1[which(x1$AB_geno == "AA"),]
x1ab <- x1[which(x1$AB_geno == "AB"),]
x1bb <- x1[which(x1$AB_geno == "BB"),]

x1aa <- x1aa[which(x1aa$Theta < 2*sd(x1aa$Theta) + mean(x1aa$Theta)) |,]

test <- x1aa %>% filter(Theta %in% ((mean(x1aa$Theta)-2*sd(x1aa$Theta)):(mean(x1aa$Theta)+2*sd(x1aa$Theta))))


test<- x1aa[between(x1aa$Theta ,mean(x1aa$Theta)-2*sd(x1aa$Theta),mean(x1aa$Theta)+2*sd(x1aa$Theta)),]










y <- floor(runif(20, min=1, max=length(t$cust_id)))
y <- 1:13983
y.1 <- 13984:27966
y.2 <- 27967:41949
y.3 <- 41950:55932
#rand num selection for snps
#y <- 1:20
#y.1 <- 21:40
#y.2 <- 41:60
#y.3 <- 61:80


x.1 <- snp_extract_func(c(t$cust_id[y]))
x.2 <- snp_extract_func(c(t$cust_id[y.1]))
x.3 <- snp_extract_func(c(t$cust_id[y.2]))
x.4 <- snp_extract_func(c(t$cust_id[y.3]))

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

final_x12 <- rbind(final_x1, final_x2)

setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/Axiom_Data_Cleaning_Outputs")
#saveRDS(final_x1, file = "final_x1.rds")
#saveRDS(final_x2, file = "final_x2.rds")
#saveRDS(final_x3, file = "final_x3.rds")
#saveRDS(final_x4, file = "final_x4.rds")

setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/Axiom_Data_Cleaning_Outputs")
new_x1 <- readRDS("x1.rds")
new_x2 <- readRDS("x2.rds")
new_x3 <- readRDS("x3.rds")

