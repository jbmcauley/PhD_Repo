#Statistics Lab Questions: MSc Quant Gen
#Basic Stats
#Written by John McAuley 2/10/2019 University of Edinburgh

rm(list = ls())
#1-----------------------------------------------------------------------
b <- 1:5
choose(5,3) #10 combinations
mean(b)  #3
bdistn <- combn(b,3)
#bdistn <- data.frame(t(bdistn))
meanbdistn <- apply(bdistn, 2,mean)
mean(meanbdistn) #3 which is equal to the pop mean!
medianbdistn <- apply(bdistn, 2, median)
median(b)
median(medianbdistn) #Also equal to 3!


var(meanbdistn) #
var(medianbdistn) #Median has greater variance


#2-----------------------------------------------------------------------
Carriers <- c(19, 53)
NonCarriers <- c(497, 829)
rbind(Carriers, NonCarriers)
ChiCarriers <- (19-26.58)^2/26.58+(53-45.42)^2/45.42+(497-489.42)^2/489.42+(829-836.58)^2/836.58
#No association
#-----------------------------------------------------------------------

#3----------------------------------------------------------------------
ChiDie <- (19/102-1/6)^2/(1/6)+(21/102-1/6)^2/(1/6)+(17/102-1/6)^2/(1/6)+(10/102-1/6)^2/(1/6)+(14/102-1/6)^2/(1/6)+(21/102-1/6)^2/(1/6)
#Yes all faces occur same prob
#-----------------------------------------------------------------------

#4----------------------

ChiRain <- (139-134)^2/134+(130-134)^2/134+(131-134)^2/134+(141-134)^2/134+(140-134)^2/134+(125-134)^2/134+(132-134)^2/134
Rain <- c(139,130,131,141,140,125,132)
Expected <-c(134,134,134,134,134,134,134)
raindiff <- Rain-Expected
rainfall <- data.frame(Rain,Expected)
chisq.test(rainfall)
sum(Rain)
mean(Rain)  
t.test(raindiff)

#5------------------------------------------------------
medidif <- c(1,3,-2,4,-1,3,6,4,2,0)
mean(medidif)
sqrt(sd(medidif)^2/10)  
2/.78
t.test(medidif)

#6-------------------------------------------------------
Darwin <- read.csv("darwin.csv")
t.test(Darwin$height.cross,Darwin$height.self)

#7----
flower <- c(705,224)
expected <- c(.75*929,.25*929)
flwr <- data.frame(flower,expected)
chisq.test(flwr)
