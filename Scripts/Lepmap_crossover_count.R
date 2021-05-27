setwd("C:/Users/s1945757/Downloads/Sparrow/scripts")
install.packages("tidyr")
library(tidyr)
library(dplyr)

# data with crossovers counts

all = data.frame(matrix(vector(), 4914, 0)) # 23610 is the nr of rows, which is the nr of gamets. you can either try to find out what it should be, or just run it with this nr and you will get an error saying that "all" and "x" dont match and you get the nr of rows in x and can put that here and run it again 


for (a in c(1:15,17:29)) { # 18 is the number of chrosmosomes in my dataset
  # Read offspring haplotype columns from lepmap output My files are called f.eks chr1_evaluatorder.txt
  file <- read.table(paste("chr", a, "_evaluateorder.txt", sep = ""), 
                    skip = 3, colClasses = "character")[,c(1,2,3,seq(7,5181,3))] #6510 is the total nr of columns in the lepmap output. try to run just: file <- read.table(paste("chr1_evaluateorder.txt", sep = ""), skip = 3, colClasses = "character") first and you see the nr of columns in the "file" dataframe. 
  #file <- read.table(paste("chr1_evaluateorder.txt", sep = ""), skip = 3, colClasses = "character")
  # split columns into individual haplotypes. 
  
  haplo <- file[,c(1,2,3)]
  
  for (i in 4:ncol(file)) {
    sep <- as.data.frame((file[,i]))
    colnames(sep) <- "col"
    sep <- separate(sep, col, into = c(paste0("h",1:nchar(file[1,i]),".",i)),sep = 1:nchar(file[1,i]))
    
    haplo <- cbind(haplo,sep)
  }
  
  N_COs <- as.data.frame(ifelse(haplo[,i]!= lag(haplo[,i], n = 1L, default = NA, order_by = NULL),1,0))
  N_COs[is.na(N_COs)] <- 0
  N_COs <- haplo
  
  for (i in 4:ncol(haplo)) {
    # Mb positions for crossovers
    N_COs[,i] <- ifelse(haplo[,i]!= lag(haplo[,i], n = 1L, default = NA, order_by = NULL),1,0)
  }
  
  
  x <- as.data.frame(colSums(N_COs[,c(4:ncol(N_COs))], na.rm = T))
  colnames(x) <- a
  all <- cbind(all,x)
  
  
}


# the final dataframe called "all" has one row per gamet (2rows per individual) and one column per chromosomes with the nr of crossovers.



# split into maternal and paternal gamets. row 1,3,5 etc are the paternal gamets and 2,4,6 etc are the maternal.. so change the numbers to fit your data

setwd("C:/Users/s1945757/Downloads/Sparrow/input")
lep_ped <- read.table("lepmap_pedigree_readable.txt", header = TRUE)


names.all <- as.data.frame(row.names(all))
colnames(names.all) <- "OG"
names.all$names.all <- gsub("^.*\\.","", names.all$OG)
names.all %>% 
  group_by(names.all) %>% 
  filter(row_number() <= .5 * n())

names.all.Paternal <- names.all %>% 
  group_by(names.all) %>% 
  filter(row_number() <= .5 * n())

names.all.Paternal <- as.data.frame(names.all.Paternal)
names.all$OG <- as.character(names.all$OG)
names.all.Paternal$OG <- as.character(names.all.Paternal$OG)
names.all.Maternal <- names.all[which(!names.all$OG %in% names.all.Paternal$OG),]

all$gamete <- row.names(all)
maternal <- all[which(!all$gamete %in% names.all.Paternal$OG),]
paternal <- all[which(all$gamete %in% names.all.Paternal$OG),]
maternal$sum<- rowSums(maternal[1:28])
paternal$sum<- rowSums(paternal[1:28])
#maternal$gamete <- NULL
#paternal$gamete <- NULL
rm(names.all.Maternal)
rm(names.all.Paternal)

hist(c(maternal$sum,paternal$sum), breaks = 65)
hist(paternal$sum, breaks = 50)

# PS: the order of the individuals is the same as in the lepmap pedigree file. so if you exclude granparents and parents from the pedigree you can merge the crossover counts with the correct individuals. 

hist(all$sum, breaks = 50)
hist(c(allpaternal$sum,allmaternal$sum), breaks = 50)
hist(allmaternal$sum, breaks = 50)
hist(allpaternal$sum, breaks = 50)
write.table(allpaternal, file = "Paternal_Gametes.txt", row.names = FALSE,col.names = TRUE)
write.table(allmaternal, file = "Maternal_Gametes.txt", row.names = FALSE,col.names = TRUE)

setwd("C:/Users/s1945757/Downloads/Sparrow/input")
lep_ped <- read.table("lepmap_pedigree_readable.txt", header = TRUE)



