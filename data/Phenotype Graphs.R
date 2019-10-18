#Phenotype Visualization
#John McAuley


library(dplyr)
library(ggplot2)

adultmorph <- read.csv("AdultMorphology-pre4_SJ.csv")
phenodata <- read.csv(("SNPtypedind_PhenotypicData.csv"))
adultmorph$sex <- as.factor(adultmorph$sex)
adultmorph$age <- as.factor(adultmorph$age)


(graph <- ggplot(adultmorph, aes(age, mass)) +
    geom_boxplot(aes(fill = age))+
    theme_bw()+ 
    ylab("Mass of individual\n")+
    xlab("\nSex"))
summary(adultmorph)

aggregate(adultmorph[,12:16],list(adultmorph$ringnr), mean)

adultmorph %>% 
  group_by(ringnr) %>% 
  summarise_at(vars(-sex, -hatchisland, -hatchyear, -year, -month, -day, -island, -islandname,-age,-maxadyear), funs(mean(.,na.rm=TRUE)))
  
