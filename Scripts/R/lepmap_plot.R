
setwd("C:/Users/s1945757/Downloads/Sparrow/output")

# Results from lepmap

# plot cM to Mb positions

chr <- 13

# Map file with Mb positions 
setwd("C:/Users/s1945757/Downloads/Sparrow/input")

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



chr <- 9

# Map file with Mb positions 
setwd("C:/Users/s1945757/Downloads/Sparrow/input")

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
plot(map_Mb_cM_9$marker_n, map_Mb_cM_9$avg_length, col="black", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)



# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome ", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)



#Overlay Plots between Crimap and Lepmap

#Order
ggplot()+
  geom_point(data = sparrow.map9, aes(x = Order, y = cMPosition), color = "maroon3")+    
  geom_point(data = map_Mb_cM_9, aes(y = avg_length, x = marker_n))
ggplot()+
        geom_point(data = sparrow.map13, aes(x = Order, y = cMPosition), color = "maroon3")+    
        geom_point(data = map_Mb_cM_13, aes(y = avg_length, x = marker_n))

#BP
ggplot()+
        geom_point(data = sparrow.map9, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
        geom_point(data = map_Mb_cM_9, aes(y = avg_length, x = bp/1000000))+
        labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")
        
ggplot()+
        geom_point(data = sparrow.map13, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
        geom_point(data = map_Mb_cM_13, aes(y = avg_length, x = bp/1000000))+
        labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")

p <- ggplot()+
  geom_point(data = map_Mb_cM.all, aes(y = avg_length, x = bp/1000000))+
  geom_point(data = sparrow.map.all, aes(x = BP/1000000, y = cMPosition), color = "maroon3")+    
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")
p + facet_wrap(~Chr, scales = "free")

Lep_Cri_combined <- merge(map_Mb_cM.all, sparrow.map.all, by = "Chr")



p <- ggplot()+
  geom_point(data = map_Mb_cM.all, aes(x = bp/1000000, y = female_length), color = "red")+    
  geom_point(data = map_Mb_cM.all, aes(y = male_length, x = bp/1000000), color = "blue")+
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")
p + facet_wrap(~Chr, scales = "free")

