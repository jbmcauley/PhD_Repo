

# Results from lepmap

# plot cM to Mb positions

chr <- 9

# Map file with Mb positions 
map_Mb <- read.table(paste("chr",chr, "_evaluateorder.map", sep=""))
colnames(map_Mb) <- c("marker_n","chr","marker_name","V4", "bp")

# lepmap results file  
lepmap_cM <- read.table(paste("chr", chr, "_evaluateorder.txt", sep=""), skip = 3)[,c(1,2,3)]
colnames(lepmap_cM) <- c("marker_n","male_length","female_length")

# merge the two
map_Mb_cM <- merge(lepmap_cM,map_Mb)


plot(map_Mb_cM$marker_n,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Marker order", ylab = "cM", main = paste("Chromosome", chr, sep = ""), pch = 20)
points(map_Mb_cM$marker_n,map_Mb_cM$male_length, col="blue", pch = 20)


# Plot cM positions on Mb positions
plot(map_Mb_cM$bp,map_Mb_cM$female_length, col="maroon1", 
     xlab = "Mb", ylab = "cM", main = paste("Chromosome", chr, sep = ""), pch = 20)
points(map_Mb_cM$bp,map_Mb_cM$male_length, col="blue", pch = 20)