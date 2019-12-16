#Stats Problems: Maximum Likelihood and Bayes
#John McAuley

#Problem 1
logl <- function(x)
  -1997*log(2+x)-1810*log(1-x)-32*log(x)

library(stats4)

fit <- mle(logl, start = list(x=0.5), method = 'L', lower = 0.001, upper = .999)
summary(fit)
confint(fit)
print(ML <- coef(fit))
print(CI <- confint(fit))

#Problem 2
logl <- function(x)
  1997*log(2+x)+1810*log(1-x)+32*log(x)
curve(logl, from = 0.01, to = 0.07, las = 1, ann = FALSE)
abline(v=CI)
lines(CI, logl(CI))
LRT <- 2*(logl(ML) - logl(CI))
round(LRT, 2)


#Problem 3
logl <- function(p,q)
  (212+39)*log(p)+(103+39)*log(q)+212*log(p+2*(1-p-q))+103*log(q+2*(1-p-q))+2*148*log(1-p-q)
minuslogl <- function(p,q) {-logl(p,q)}
fit <- mle(minuslogl, start = list(p=1/3, q=1/3), method = 'L', lower = 0.001, upper = .999)
summary(fit)

#Problem 4
curve(dbeta(x,1,19), from = 0, to = 0.35, las =1, ann = FALSE, lty = 2)
curve(dbeta(x, 3, 27), add = TRUE)
legend('topright', legend = c('posterior', 'prior'), lty =1:2, inset = 0.1)

#Probability that the proportion of defective beakers exceeds 10%
pbeta(.1,1,19, lower.tail = FALSE) #Prior 
pbeta(.1,3,27, lower.tail = FALSE) #Posterior

#Mean of a beta with parameters: with parameters a and b: a/(a+b)
#Medians:
qbeta(0.5, 1, 19)
qbeta(.5, 3, 27)
#Variance is: m*(1-m)/(a+b+1), where m is the mean of the distn


#Problem 5 
logl <- function(x)
  1997 * log(2+x) + 1810 * log(1-x) + 32* log(x)
reps <- 3000
proposal <- runif(reps, 0.02, 0.06)
logprobs <- logl(proposal)
probs <- exp(logprobs - max(logprobs))
I <- rbinom(reps, 1, probs)
mean(I)
theta <- proposal[I > 0]
h <- hist(theta, plot = FALSE)
plot(density(theta), las = 1, ann = FALSE)
lines(h, freq = FALSE)
