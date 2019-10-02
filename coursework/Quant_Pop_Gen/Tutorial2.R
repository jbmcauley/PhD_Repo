#Quantitative Genetics MSc Tutorial 2
#Practicing application of important concepts and equations
#Written by John McAuley 30/9/2019 University of Edinburgh

#Question 1-------------------

q2 <- 130/(130+300+180)
pq2 <- 300/(300+130+180)
p2 <- 180/(300+130+180)
sum(q2 + pq2 + p2)
q <- (2*130+300)/(2*180+2*130+2*300)
p <- (2*180+300)/(2*180+2*130+2*300)
p+q

q2.exp <- q^2
pq2.exp <- 2*p*q
p2.exp <- p^2

#p^2 = p^2 +fpq
#p2 = p^2 + f*p*q
#(p2 - p^2)/(pq)=f
F.a1a1 <- (p2-p^2)/(p*q) 
F.a2a2 <- (q2-q^2)/(p*q) 
F.a1a2 <- -(pq2/(2*p*q)-1)

#Question 2-------------------
rm(list = ls())
q <- .05
q2 <- q^2 #Expected Freq of female w/ colorblindness
pq2 <- (2*(1-q)*q) #Expected freq of females
p <- 1-q
p2 <- p^2

#Chance colorblind man has son who is color blind
son <- .05
dter <- .05 

#Question 3---------------------
rm(list = ls())
#   qprime = q^2 + 2pq*(.75)
#          = q2 + 1.5pq

#delta-q = qprime - q
#        = q^2 + 1.5pq - q
#        = q(q + 1.5p - 1)
#        = q(q + p + 0.5p - 1)
#        *p+q=1*
#        = 0.5pq

#Normally w2. = qw22 + pw12
#     w22 & w12 = 1
#             = q + 1.5p


# Delta-q = 0 = p*q(w2. - w1.)/W
#             = pq((q(1-s)+p)-(p+q))/W
#             = (3/32)s/W
#
(1/8) * (3/4)
3/36
3/32
1/1.5
