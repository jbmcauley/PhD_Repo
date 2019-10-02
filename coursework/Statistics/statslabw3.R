#Statistics Lab Questions: MSc Quant Gen
#Basic Stats
#Written by John McAuley 2/10/2019 University of Edinburgh

rm(list = ls())

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


