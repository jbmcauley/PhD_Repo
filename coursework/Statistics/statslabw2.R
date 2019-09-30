rm(list = ls())
Q1a <- (3/5)*(2/4)
Q1b <- (3/5)*(2/4)+(2/5)*(3/4)
Q2a <- dbinom(1,4,.25)*4
#Pr(A|B) = [Pr(B|A)*Pr(A)]/Pr(B)
#Pr(A) = 2/3; Pr(B|A) = 3/4; Pr(B) = 7/8
Q2b <- ((1)*(2/3))/(7/8)
Q3 <- 'c1 + c2 - 2(c1*c2)'
#c1(1-c2)+c2(1-c1)
#pbinom(x,n,p)
Q4 <- 1 - pbinom(5,10,.5)
Q5a1 <- pnorm(.75)-(1-pnorm(.75))
Q5a2 <- pnorm(1.0)-(1-pnorm(1.0))
Q5a3 <- pnorm(1.5)-(1-pnorm(1.5))
Q5a4 <- pnorm(2.0)-(1-pnorm(2.0))
Q5b_1 <- 0.68
Q5b_2 <- 1.65
Q5b_3 <- 1.96
Q5c_Q1 <- pnorm(-.67) 
Q5c_Q2 <- pnorm(0)
Q5c_Q3 <- pnorm(.67)
#Q6
#m <- 750
#s <- 100
Q6a <- pnorm((850-750)/100)-pnorm((650-750)/100)
Q6b <- pnorm((600-750)/100)

#Q7
x <- 0:16
print(binomprobs <- dbinom(x, 16, 1/4))
m <- weighted.mean(x, binomprobs)
v <- weighted.mean(x^2, binomprobs) - m^2
c(mean = m, variance = v)
barplot(binomprobs, names = 0:16)

#Q8
levels <- c(.25, .5, .75)
quartiles <- qnorm(levels)
curve(dnorm, from = -3, to = 3)
abline(v = quartiles)
#cumulative version
curve(pnorm, from = -3, to=3)
abline(h=levels, v = quartiles)
