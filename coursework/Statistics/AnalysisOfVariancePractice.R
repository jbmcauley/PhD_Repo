#Practice Problems from Statistics: Analysis of Variance
#John MCAuley, Oct 23 2019

rm(list = ls())

#Question 1
control <- c(99,95,120,112,80,106,98,102)
treatment <- c(108,110,105,131,104,96,116,118)
diff <- controls - treatment
mean(diff)
t.test(diff)
t.test(controls, treatment)
lambs <- data.frame(controls, treatment)
t.test(lambs)
summary(aov(control~ treatment))



#Question 2
lv <- c(0,1,2,3,4,5,6,7)
Yield <- c(36,54,39,57,66,60,63,75)
irrig <- data.frame(lv,Yield)
irrig <- transform(irrig, lv = factor(lv))
fit <- aov(Yield ~ lv)
summary(fit)




#Question 3
bristle <- read.csv("Bristles.csv")
colnames(bristle) <- c("Genotype", "Cy", "Me", "Bristles")
fit <- lm(Bristles ~ Genotype + Cy:Genotype + Me:Genotype, data = bristle)
summary(fit)
aov.fit <- aov(Bristles ~ Genotype, data = bristle)
summary(aov.fit)
aov.fit <- aov(Bristles ~ Genotype*Cy, data = bristle)
summary(aov.fit)



#Question 4
bobbins <- read.csv("bobbins.csv")
aov.fit <-  aov(strength ~ bobbin, data = bobbins)
summary(aov.fit, split = list(bobbins$bobbin))
fit <- lm(strength ~ bobbin, data = bobbins)
summary(fit)
model.tables(aov.fit, "mean", se = TRUE)
