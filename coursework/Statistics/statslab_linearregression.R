#----
#Question 1
age <- c(35,45,55,65,75)
mbp <- c(114,124,143,158,166)

fit <- lm(mbp ~ age)
anova(fit)
summary(fit)
confint(fit)
plot(age, mbp)
abline(fit)

#----
#Question 2
Irrig <- c(0:7)
yield <- c(12,18,13,19,22,20,21,25)
N <- 8
x <- 28
y <- 150
x2 <- 140
y2 <- 2948
xy <- 590
Sxx <- x2-(x^2/N)
Sxy <- xy - (x*y/N)
Syy <- y2 - (y^2/N)
b <- Sxy/Sxx  #1.55
mean(yield)
mean(Irrig)
b0 <- 13.32
fit <- lm(yield ~ Irrig)
summary(fit)
new <- data.frame(x = seq(0,8,by =0.5))
preds <- predict(fit,new, interval = 'conf')
matplot(new$x, preds, ann = FALSE, las = 1)

#----
#Question 3
brains <- read.csv("brains.csv")
fit <- lm(log(brains$brain.wt) ~ log(brains$body.wt))
summary(fit)
plot(fit)
plot(y=log(brains$brain.wt), x=log(brains$body.wt)) #b0 = 2.13, b1 = 0.75
abline(fit)
anova(fit)
R2 <- 336.19/(336.19+28.92)
sqrt(R2)
#y = mx + b
#Doubling X -> 2X*b+bo = 2*0.75(X) +2.13
#Doubling X leads to a 1.5 times increase in Y
#Very slight devations due to potential outliers

#----
#Q4

adata <- read.csv("anscombe.csv")

y1lm <- lm(adata$y1 ~adata$x) #Okay
y2lm<- lm(adata$y2 ~ adata$x) #Not okay, concave residual vs fitted values
y3lm<- lm(adata$y3 ~adata$x) #Not okay, linear non-zero residual slope vs fitted values

plot(y3lm)

#----
#Q5
rm(list = ls())
library("MASS")
data(forbes, package = "MASS")
help(forbes)

plot(forbes$bp, forbes$pres)
fit <- lm(forbes$pres ~ forbes$bp)
summary(fit)
plot(fit) #Not great looking: pt12 potential outlier?
forbes2 <- forbes[-12,]
fit2 <- lm(forbes2$pres ~ forbes2$bp)
summary(fit2)
plot(fit2)
#Not really a major improvement to diagnostic graphs..

fitlog <- lm(log(forbes$pres) ~ forbes$bp)
summary(fitlog)
plot(fitlog) #pt 12 more clearly an outlier

#----
#Q6

broth <- c(71,68,66,67,70,71,70,73,72,65,66)
sist<- c(69,64,65,63,65,62,65,64,66,59,62)

cor.test(broth,sist)
cor.test(sist,broth)
plot(broth,sist)
fit <- lm(sist~broth)
summary(fit)
sqrt(.31)
