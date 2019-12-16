#----- Problem 1
choose(8,4)
choose(10,4)
choose(10,4)


TeaTasting <-
  matrix(c(5, 1, 1, 3),
         nrow = 2,
         dimnames = list(Guess = c("Milk", "Tea"),
                         Truth = c("Milk", "Tea")))
fisher.test(TeaTasting, alternative = "greater")
#There is no evidence to reject Null hypothesis

#----- problem2
binom.test(3,10, p=.5)
#p-value = 0.3438



#-----problem 3
x <- c(5,-1,7,13,3,14,-4,10)
wilcox.test(x, alternative = "greater")


#------problem 4
x<- c(61,-84,10,20,7,29,35,51,17,36,70,30,93,75,-60)
x <- x/10
wilcox.test(x, alternative = "two.sided")

#-----problem5 
x <- c(20,23.7,25.1,21,26.4)
y <- c(18.5,22,24.3,19.6,20.5)
z <- y-x
wilcox.test(z, alternative = "less")


#-----problem 6
x <- c(16,22,21,19,15,13,23,17,20,29,18,25)
median(x)

#Check values in R
n <- length(x)
d <- qbinom(0.025, n, 0.5)
x[c(d,n+1-d)]
1-2*pbinom(d-1,n,0.5)


#Values achieved by hand
#rank of lowerlimit
12/2 + sqrt(12)*1.96/2
#approximately ~3 so we will use 16

#rank of upper limit
1+ 12/2 + sqrt(12)*1.96/2
#approximately 10 so will use 23
#95% confidence interval is (16,23)



