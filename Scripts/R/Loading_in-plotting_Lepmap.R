setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 1
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_1 <- map_Mb_cM
plot(map_Mb_cM_1$marker_n,map_Mb_cM_1$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 2
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_2 <- map_Mb_cM
plot(map_Mb_cM_2$marker_n,map_Mb_cM_2$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)






setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 3
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_3 <- map_Mb_cM
plot(map_Mb_cM_3$marker_n,map_Mb_cM_3$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 4
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_4 <- map_Mb_cM
plot(map_Mb_cM_4$marker_n,map_Mb_cM_4$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 5
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_5 <- map_Mb_cM
plot(map_Mb_cM_5$marker_n,map_Mb_cM_5$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)



setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 6
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_6 <- map_Mb_cM
plot(map_Mb_cM_6$marker_n,map_Mb_cM_6$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 7
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_7 <- map_Mb_cM
plot(map_Mb_cM_7$marker_n,map_Mb_cM_7$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)







setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 8
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_8 <- map_Mb_cM
plot(map_Mb_cM_8$marker_n,map_Mb_cM_8$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 9
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_9 <- map_Mb_cM
plot(map_Mb_cM_9$marker_n,map_Mb_cM_9$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 10
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_10 <- map_Mb_cM
plot(map_Mb_cM_10$marker_n,map_Mb_cM_10$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)







setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 11
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_11 <- map_Mb_cM
plot(map_Mb_cM_11$marker_n,map_Mb_cM_11$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)




setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 12
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_12 <- map_Mb_cM
plot(map_Mb_cM_12$marker_n,map_Mb_cM_12$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 13
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_13 <- map_Mb_cM
plot(map_Mb_cM_13$marker_n,map_Mb_cM_13$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)



setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 14
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_14 <- map_Mb_cM
plot(map_Mb_cM_14$marker_n,map_Mb_cM_14$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 15
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_15 <- map_Mb_cM
plot(map_Mb_cM_15$marker_n,map_Mb_cM_15$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)




setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 16
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_16 <- map_Mb_cM
plot(map_Mb_cM_16$marker_n,map_Mb_cM_16$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 17
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_17 <- map_Mb_cM
plot(map_Mb_cM_17$marker_n,map_Mb_cM_17$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)




setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 18
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_18 <- map_Mb_cM
plot(map_Mb_cM_18$marker_n,map_Mb_cM_18$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)




setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 19
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_19 <- map_Mb_cM
plot(map_Mb_cM_19$marker_n,map_Mb_cM_19$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 20
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_20 <- map_Mb_cM
plot(map_Mb_cM_20$marker_n,map_Mb_cM_20$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)




setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 21
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_21 <- map_Mb_cM
plot(map_Mb_cM_21$marker_n,map_Mb_cM_21$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)




setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 22
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_22 <- map_Mb_cM
plot(map_Mb_cM_22$marker_n,map_Mb_cM_22$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)




setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 23
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_23 <- map_Mb_cM
plot(map_Mb_cM_23$marker_n,map_Mb_cM_23$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)







setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 24
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_24 <- map_Mb_cM
plot(map_Mb_cM_24$marker_n,map_Mb_cM_24$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 25
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_25 <- map_Mb_cM
plot(map_Mb_cM_25$marker_n,map_Mb_cM_25$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 26
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_26 <- map_Mb_cM
plot(map_Mb_cM_26$marker_n,map_Mb_cM_26$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)




setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 27
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_27 <- map_Mb_cM
plot(map_Mb_cM_27$marker_n,map_Mb_cM_27$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)



setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 28
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_28 <- map_Mb_cM
plot(map_Mb_cM_28$marker_n,map_Mb_cM_28$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)





setwd("C:/Users/s1945757/Downloads/Sparrow/input")
chr <- 29
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))

colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)
map_Mb_cM$avg_length <- (lepmap_cM$male_length + lepmap_cM$female_length)/2
points(map_Mb_cM$marker_n,map_Mb_cM$avg_length, col="black", 
       xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)

map_Mb_cM_29 <- map_Mb_cM
plot(map_Mb_cM_29$marker_n,map_Mb_cM_29$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)











#BP
ggplot()+
        geom_point(data = sparrow.map1, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
        geom_point(data = map_Mb_cM_1, aes(y = avg_length, x = bp/1000000))+
        labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
        ggtitle('Chr 1: Lepmap vs Crimap')

ggplot()+
        geom_point(data = sparrow.map2, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
        geom_point(data = map_Mb_cM_2, aes(y = avg_length, x = bp/1000000))+
        labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
        ggtitle('Chr 2: Lepmap vs Crimap')

ggplot()+
        geom_point(data = sparrow.map3, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
        geom_point(data = map_Mb_cM_3, aes(y = avg_length, x = bp/1000000))+
        labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
        ggtitle('Chr 3: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map4, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_4, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 4: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map5, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_5, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 5: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map6, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_6, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 6: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map7, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_7, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 7: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map8, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_8, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 8: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map9, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_9, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 9: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map10, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_10, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 10: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map11, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_11, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 11: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map12, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_12, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 12: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map13, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_13, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 13: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map14, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_14, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 14: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map15, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_15, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 15: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map17, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_17, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 17: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map18, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_18, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 18: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map19, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_19, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 19: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map20, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_20, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 20: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map21, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_21, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 21: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map22, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_22, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 22: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map23, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_23, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 23: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map24, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_24, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 24: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map26, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_26, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 26: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map27, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_27, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 27: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map28, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_28, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 28: Lepmap vs Crimap')

ggplot()+
  geom_point(data = sparrow.map29, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_29, aes(y = avg_length, x = bp/1000000))+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")+
  ggtitle('Chr 29: Lepmap vs Crimap')



sparrow.map.all <- rbind(sparrow.map1,sparrow.map2,sparrow.map3,sparrow.map4,sparrow.map5,sparrow.map6,sparrow.map7,sparrow.map8,
                         sparrow.map9,sparrow.map10,sparrow.map11,sparrow.map12,sparrow.map13,sparrow.map14,sparrow.map15,sparrow.map17,
                         sparrow.map18,sparrow.map19,sparrow.map20,sparrow.map21,sparrow.map22,sparrow.map23,sparrow.map24,sparrow.map26,
                         sparrow.map27,sparrow.map28,sparrow.map29)
map_Mb_cM.all <- rbind(map_Mb_cM_1,map_Mb_cM_2,map_Mb_cM_3,map_Mb_cM_4,map_Mb_cM_5,map_Mb_cM_6,map_Mb_cM_7,map_Mb_cM_8,
                         map_Mb_cM_9,map_Mb_cM_10,map_Mb_cM_11,map_Mb_cM_12,map_Mb_cM_13,map_Mb_cM_14,map_Mb_cM_15,map_Mb_cM_17,
                         map_Mb_cM_18,map_Mb_cM_19,map_Mb_cM_20,map_Mb_cM_21,map_Mb_cM_22,map_Mb_cM_23,map_Mb_cM_24,map_Mb_cM_26,
                         map_Mb_cM_27,map_Mb_cM_28,map_Mb_cM_29)

x<- map_Mb_cM.all
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 1] <- "Chr 1"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 2] <- "Chr 2"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 3] <- "Chr 3"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 4] <- "Chr 4"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 5] <- "Chr 5"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 6] <- "Chr 6"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 7] <- "Chr 7"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 8] <- "Chr 8"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 9] <- "Chr 9"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 10] <- "Chr 10"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 11] <- "Chr 11"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 12] <- "Chr 12"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 13] <- "Chr 13"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 14] <- "Chr 14"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 15] <- "Chr 15"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 16] <- "Chr 16"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 17] <- "Chr 17"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 18] <- "Chr 18"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 19] <- "Chr 19"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 20] <- "Chr 20"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 21] <- "Chr 21"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 22] <- "Chr 22"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 23] <- "Chr 23"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 24] <- "Chr 24"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 25] <- "Chr 25"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 26] <- "Chr 26"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 27] <- "Chr 27"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 28] <- "Chr 28"
map_Mb_cM.all$chr[map_Mb_cM.all$chr == 29] <- "Chr 29"

names(map_Mb_cM.all)[names(map_Mb_cM.all) == "chr"] <- "Chr"


p <-ggplot(data = map_Mb_cM.all, aes(x= bp/1000000, y = avg_length)) +
  geom_point(size = 1) +
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")
p + facet_wrap(~Chr, scales = "free")



p <-ggplot(data = sparrow.map.all, aes(x = BP/1000000, y = cMPosition)) +
  geom_point(size = 1) +
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")

p + facet_wrap(~Chr, scales = "free")


