#Phenotype Visualization
#John McAuley


library(dplyr)
library(ggplot2)

adultmorph <- read.csv("AdultMorphology-pre4_SJ.csv")
phenodata <- read.csv(("SNPtypedind_PhenotypicData.csv"))
adultmorph$sex <- as.factor(adultmorph$sex)
adultmorph$age <- as.factor(adultmorph$age)

summary(adultmorph)

#Histograms for each phenotypic data variable
adultmorph %>% 
ggplot(aes(x = tarsus)) +
  geom_histogram(aes(y=..density..),bins = 51,colour = "darkgrey", fill = "white") +
  geom_density(colour = "blue")+
  theme_bw() +
  ylab("Frequency\n")+
  xlab("\nTarsus Length (mm)")+
 ggtitle("Histogram of Tarsus Length\n")

adultmorph %>% 
  ggplot(aes(wing)) +
  geom_histogram(aes(y=..density..), bins = 51, colour ="darkgrey", fill = "white")+
  geom_density(colour = "blue")+
  theme_bw()+
  ylab("Frequency\n")+
  xlab("\nWing length (mm)")+
  ggtitle("Histogram of Wing Length\n")

adultmorph %>% 
  ggplot(aes(billD)) +
  geom_histogram(aes(y=..density..),bins = 51, colour ="darkgrey", fill = "white")+
  geom_density(colour = "blue")+
  theme_bw()+
  ylab("Frequency\n")+
  xlab("\nBill Depth (mm)")+
  ggtitle("Histogram of Bill Depth")

adultmorph %>% 
  ggplot(aes(billL)) +
  geom_histogram(aes(y=..density..),bins = 51, colour ="darkgrey", fill = "white")+
  geom_density(colour = "blue")+
  theme_bw()+
  ylab("Frequency\n")+
  xlab("\nBill Length (mm)")+
  ggtitle("Histogram of Bill Length")

adultmorph %>% 
  ggplot(aes(mass)) +
  geom_histogram(aes(y=..density..),bins = 51, colour ="darkgrey", fill = "white")+
  geom_density(colour = "blue")+
  theme_bw()+
  ylab("Frequency\n")+
  xlab("\nMass (g)")+
  ggtitle("Histogram of Mass")

#Exploartion of the Mass of Individuals at various ages and differing Sex
(graph <- ggplot(adultmorph, aes(age, mass)) +
    geom_boxplot(aes(fill = sex))+
    theme_bw()+ 
    ylab("Mass of individual\n")+
    xlab("\nAge"))


#Visualization of the distribution of mean values from
#multiple measures of the same individuals for multiple
adultmorph.means <- aggregate(adultmorph[,12:16],list(adultmorph$ringnr), mean)

adultmorph.means %>% 
  ggplot(aes(x = tarsus)) +
  geom_histogram(aes(y=..density..),bins = 51,colour = "darkgrey", fill = "white") +
  geom_density(colour = "blue")+
  theme_bw() +
  ylab("Frequency\n")+
  xlab("\nTarsus Length (mm)")+
  ggtitle("Histogram of Mean Tarsus Length\n")

adultmorph.means %>% 
  ggplot(aes(wing)) +
  geom_histogram(aes(y=..density..), bins = 51, colour ="darkgrey", fill = "white")+
  geom_density(colour = "blue")+
  theme_bw()+
  ylab("Frequency\n")+
  xlab("\nWing length (mm)")+
  ggtitle("Histogram of Mean Wing Length\n")

adultmorph.means %>% 
  ggplot(aes(billD)) +
  geom_histogram(aes(y=..density..),bins = 51, colour ="darkgrey", fill = "white")+
  geom_density(colour = "blue")+
  theme_bw()+
  ylab("Frequency\n")+
  xlab("\nBill Depth (mm)")+
  ggtitle("Histogram of Mean Bill Depth")

adultmorph.means %>% 
  ggplot(aes(billL)) +
  geom_histogram(aes(y=..density..),bins = 51, colour ="darkgrey", fill = "white")+
  geom_density(colour = "blue")+
  theme_bw()+
  ylab("Frequency\n")+
  xlab("\nBill Length (mm)")+
  ggtitle("Histogram of Mean Bill Length")

adultmorph.means %>% 
  ggplot(aes(mass)) +
  geom_histogram(aes(y=..density..),bins = 51, colour ="darkgrey", fill = "white")+
  geom_density(colour = "blue")+
  theme_bw()+
  ylab("Frequency\n")+
  xlab("\nMass (g)")+
  ggtitle("Histogram of Mean Mass")



adultmorph %>% 
  group_by(ringnr) %>% 
  summarise_at(vars(-sex, -hatchisland, -hatchyear, -year, -month, -day, -island, -islandname,-age,-maxadyear), funs(mean(.,na.rm=TRUE)))
  
