#Week 2 Coding Club: Intro to R
#Challenge yourself section
#Written by John McAuley 2/10/2019 University of Edinburgh

#Clear R's Brain
rm(list = ls())

#Load Desired Packages
library('dplyr')
library('ggplot2')


bird_sp <- c("sparrow", 'kingfisher', 'eagle', 'hummingbird')
bird_sp <- c(bird_sp, bird_sp, bird_sp)
bird_sp
wingspan <- c(22,26,195,8,24,23,201,9,21,25,185,9)
x <- as.data.frame(cbind(bird_sp,wingspan))
x$wingspan <- as.numeric(as.character(x$wingspan))
x %>% 
  group_by(bird_sp) %>% 
  summarise(Mean = mean(wingspan))
#Using dplyr to find the mean for each species

x %>%
  group_by(bird_sp) %>%
  summarise(Mean = mean(wingspan)) %>% 
  ggplot(aes(x= bird_sp, y=Mean)) +
  geom_col() +
  theme_bw()
