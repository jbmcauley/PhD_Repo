setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/Dblxoversremoved/")



#Set up file with BP
map <- as.data.frame(data2@gtdata@map)
names(map) <- "BP"
map$SNP.Name <- row.names(map)
row.names(map) <- NULL

sparrow.map1a <- parse_map(mapfile = "crimap/chr1a.map")
sparrow.map1b <- parse_map(mapfile = "crimap/chr1b.map")
sparrow.map1c <- parse_map(mapfile = "crimap/chr1c.map")
sparrow.map1d <- parse_map(mapfile = "crimap/chr1d.map")
#Combining Chr chunks

#adding cM pos within chunks
sparrow.map1b$cMPosition <- sparrow.map1b$cMPosition + (sparrow.map1a[which(sparrow.map1b$SNP.Name[1] == sparrow.map1a$SNP.Name),]$cMPosition)
sparrow.map1c$cMPosition <- sparrow.map1c$cMPosition + (sparrow.map1b[which(sparrow.map1c$SNP.Name[1] == sparrow.map1b$SNP.Name),]$cMPosition)
sparrow.map1d$cMPosition <- sparrow.map1d$cMPosition + (sparrow.map1c[which(sparrow.map1d$SNP.Name[1] == sparrow.map1c$SNP.Name),]$cMPosition)
#adding M&F cM pos within Chunks
sparrow.map1b$cMPosition.Female <- sparrow.map1b$cMPosition.Female + (sparrow.map1a[which(sparrow.map1b$SNP.Name[1] == sparrow.map1a$SNP.Name),]$cMPosition.Female)
sparrow.map1c$cMPosition.Female <- sparrow.map1c$cMPosition.Female + (sparrow.map1b[which(sparrow.map1c$SNP.Name[1] == sparrow.map1b$SNP.Name),]$cMPosition.Female)
sparrow.map1d$cMPosition.Female <- sparrow.map1d$cMPosition.Female + (sparrow.map1c[which(sparrow.map1d$SNP.Name[1] == sparrow.map1c$SNP.Name),]$cMPosition.Female)
sparrow.map1b$cMPosition.Male <- sparrow.map1b$cMPosition.Male + (sparrow.map1a[which(sparrow.map1b$SNP.Name[1] == sparrow.map1a$SNP.Name),]$cMPosition.Male)
sparrow.map1c$cMPosition.Male <- sparrow.map1c$cMPosition.Male + (sparrow.map1b[which(sparrow.map1c$SNP.Name[1] == sparrow.map1b$SNP.Name),]$cMPosition.Male)
sparrow.map1d$cMPosition.Male <- sparrow.map1d$cMPosition.Male + (sparrow.map1c[which(sparrow.map1d$SNP.Name[1] == sparrow.map1c$SNP.Name),]$cMPosition.Male)
#Adding order within chunks
sparrow.map1b$Order <- sparrow.map1b$Order + (sparrow.map1a[which(sparrow.map1b$SNP.Name[1] == sparrow.map1a$SNP.Name),]$Order)
sparrow.map1c$Order <- sparrow.map1c$Order + (sparrow.map1b[which(sparrow.map1c$SNP.Name[1] == sparrow.map1b$SNP.Name),]$Order)
sparrow.map1d$Order <- sparrow.map1d$Order + (sparrow.map1c[which(sparrow.map1d$SNP.Name[1] == sparrow.map1c$SNP.Name),]$Order)

#Remove overlap
sparrow.map1a <- sparrow.map1a[-which(sparrow.map1a$SNP.Name %in% sparrow.map1b$SNP.Name),]
sparrow.map1b <- sparrow.map1b[-which(sparrow.map1b$SNP.Name %in% sparrow.map1c$SNP.Name),]
sparrow.map1c <- sparrow.map1c[-which(sparrow.map1c$SNP.Name %in% sparrow.map1d$SNP.Name),]

#Merge
sparrow.map1a <- merge(sparrow.map1a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map1b <- merge(sparrow.map1b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map1c <- merge(sparrow.map1c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map1d <- merge(sparrow.map1d, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map1 <- rbind(sparrow.map1a,sparrow.map1b,sparrow.map1c,sparrow.map1d)
sparrow.map1[,c(5,6,7,8,10,11,12)] <- NULL

rm(sparrow.map1a)
rm(sparrow.map1b)
rm(sparrow.map1c)
rm(sparrow.map1d)





#Chr2
sparrow.map2a <- parse_map(mapfile = "crimap/chr2a.map")
sparrow.map2b <- parse_map(mapfile = "crimap/chr2b.map")
sparrow.map2c <- parse_map(mapfile = "crimap/chr2c.map")
sparrow.map2d <- parse_map(mapfile = "crimap/chr2d.map")
sparrow.map2e <- parse_map(mapfile = "crimap/chr2e.map")
#sparrow.map2f <- parse_map(mapfile = "crimap/chr2f.map")

sparrow.map2b$cMPosition <- sparrow.map2b$cMPosition + (sparrow.map2a[which(sparrow.map2b$SNP.Name[1] == sparrow.map2a$SNP.Name),]$cMPosition)
sparrow.map2c$cMPosition <- sparrow.map2c$cMPosition + (sparrow.map2b[which(sparrow.map2c$SNP.Name[1] == sparrow.map2b$SNP.Name),]$cMPosition)
sparrow.map2d$cMPosition <- sparrow.map2d$cMPosition + (sparrow.map2c[which(sparrow.map2d$SNP.Name[1] == sparrow.map2c$SNP.Name),]$cMPosition)
sparrow.map2e$cMPosition <- sparrow.map2e$cMPosition + (sparrow.map2d[which(sparrow.map2e$SNP.Name[1] == sparrow.map2d$SNP.Name),]$cMPosition)
#sparrow.map2f$cMPosition <- sparrow.map2f$cMPosition + (sparrow.map2e[which(sparrow.map2f$SNP.Name[1] == sparrow.map2e$SNP.Name),]$cMPosition)

sparrow.map2b$cMPosition.Female <- sparrow.map2b$cMPosition.Female + (sparrow.map2a[which(sparrow.map2b$SNP.Name[1] == sparrow.map2a$SNP.Name),]$cMPosition.Female)
sparrow.map2c$cMPosition.Female <- sparrow.map2c$cMPosition.Female + (sparrow.map2b[which(sparrow.map2c$SNP.Name[1] == sparrow.map2b$SNP.Name),]$cMPosition.Female)
sparrow.map2d$cMPosition.Female <- sparrow.map2d$cMPosition.Female + (sparrow.map2c[which(sparrow.map2d$SNP.Name[1] == sparrow.map2c$SNP.Name),]$cMPosition.Female)
sparrow.map2e$cMPosition.Female <- sparrow.map2e$cMPosition.Female + (sparrow.map2d[which(sparrow.map2e$SNP.Name[1] == sparrow.map2d$SNP.Name),]$cMPosition.Female)
sparrow.map2f$cMPosition.Female <- sparrow.map2f$cMPosition.Female + (sparrow.map2e[which(sparrow.map2f$SNP.Name[1] == sparrow.map2e$SNP.Name),]$cMPosition.Female)
sparrow.map2b$cMPosition.Male <- sparrow.map2b$cMPosition.Male + (sparrow.map2a[which(sparrow.map2b$SNP.Name[1] == sparrow.map2a$SNP.Name),]$cMPosition.Male)
sparrow.map2c$cMPosition.Male <- sparrow.map2c$cMPosition.Male + (sparrow.map2b[which(sparrow.map2c$SNP.Name[1] == sparrow.map2b$SNP.Name),]$cMPosition.Male)
sparrow.map2d$cMPosition.Male <- sparrow.map2d$cMPosition.Male + (sparrow.map2c[which(sparrow.map2d$SNP.Name[1] == sparrow.map2c$SNP.Name),]$cMPosition.Male)
sparrow.map2e$cMPosition.Male <- sparrow.map2e$cMPosition.Male + (sparrow.map2d[which(sparrow.map2e$SNP.Name[1] == sparrow.map2d$SNP.Name),]$cMPosition.Male)
#sparrow.map2f$cMPosition.Male <- sparrow.map2f$cMPosition.Male + (sparrow.map2e[which(sparrow.map2f$SNP.Name[1] == sparrow.map2e$SNP.Name),]$cMPosition.Male)

#Adding order within chunks
sparrow.map2b$Order <- sparrow.map2b$Order + (sparrow.map2a[which(sparrow.map2b$SNP.Name[1] == sparrow.map2a$SNP.Name),]$Order)
sparrow.map2c$Order <- sparrow.map2c$Order + (sparrow.map2b[which(sparrow.map2c$SNP.Name[1] == sparrow.map2b$SNP.Name),]$Order)
sparrow.map2d$Order <- sparrow.map2d$Order + (sparrow.map2c[which(sparrow.map2d$SNP.Name[1] == sparrow.map2c$SNP.Name),]$Order)
sparrow.map2e$Order <- sparrow.map2e$Order + (sparrow.map2d[which(sparrow.map2e$SNP.Name[1] == sparrow.map2d$SNP.Name),]$Order)
#sparrow.map2f$Order <- sparrow.map2f$Order + (sparrow.map2e[which(sparrow.map2f$SNP.Name[1] == sparrow.map2e$SNP.Name),]$Order)

#Remove overlap
sparrow.map2a <- sparrow.map2a[-which(sparrow.map2a$SNP.Name %in% sparrow.map2b$SNP.Name),]
sparrow.map2b <- sparrow.map2b[-which(sparrow.map2b$SNP.Name %in% sparrow.map2c$SNP.Name),]
sparrow.map2c <- sparrow.map2c[-which(sparrow.map2c$SNP.Name %in% sparrow.map2d$SNP.Name),]
sparrow.map2d <- sparrow.map2d[-which(sparrow.map2d$SNP.Name %in% sparrow.map2e$SNP.Name),]
#sparrow.map2e <- sparrow.map2e[-which(sparrow.map2e$SNP.Name %in% sparrow.map2f$SNP.Name),]

#Merge
sparrow.map2a <- merge(sparrow.map2a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2b <- merge(sparrow.map2b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2c <- merge(sparrow.map2c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2d <- merge(sparrow.map2d, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2e <- merge(sparrow.map2e, map, by.x = "SNP.Name",by.y = "SNP.Name")
#sparrow.map2f <- merge(sparrow.map2f, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map2 <- rbind(sparrow.map2a,sparrow.map2b,sparrow.map2c,sparrow.map2d,sparrow.map2e)
sparrow.map2[,c(5,6,7,8,10,11,12)] <- NULL

rm(sparrow.map2a)
rm(sparrow.map2b)
rm(sparrow.map2c)
rm(sparrow.map2d)
rm(sparrow.map2e)
rm(sparrow.map2f)




#Chr3
sparrow.map3a <- parse_map(mapfile = "crimap/chr3a.map")
sparrow.map3b <- parse_map(mapfile = "crimap/chr3b.map")
sparrow.map3c <- parse_map(mapfile = "crimap/chr3c.map")
sparrow.map3d <- parse_map(mapfile = "crimap/chr3d.map")
#sparrow.map3e <- parse_map(mapfile = "crimap/chr3e.map")


sparrow.map3b$cMPosition <- sparrow.map3b$cMPosition + (sparrow.map3a[which(sparrow.map3b$SNP.Name[1] == sparrow.map3a$SNP.Name),]$cMPosition)
sparrow.map3c$cMPosition <- sparrow.map3c$cMPosition + (sparrow.map3b[which(sparrow.map3c$SNP.Name[1] == sparrow.map3b$SNP.Name),]$cMPosition)
sparrow.map3d$cMPosition <- sparrow.map3d$cMPosition + (sparrow.map3c[which(sparrow.map3d$SNP.Name[1] == sparrow.map3c$SNP.Name),]$cMPosition)
#sparrow.map3e$cMPosition <- sparrow.map3e$cMPosition + (sparrow.map3d[which(sparrow.map3e$SNP.Name[1] == sparrow.map3d$SNP.Name),]$cMPosition)

sparrow.map3b$cMPosition.Female <- sparrow.map3b$cMPosition.Female + (sparrow.map3a[which(sparrow.map3b$SNP.Name[1] == sparrow.map3a$SNP.Name),]$cMPosition.Female)
sparrow.map3c$cMPosition.Female <- sparrow.map3c$cMPosition.Female + (sparrow.map3b[which(sparrow.map3c$SNP.Name[1] == sparrow.map3b$SNP.Name),]$cMPosition.Female)
sparrow.map3d$cMPosition.Female <- sparrow.map3d$cMPosition.Female + (sparrow.map3c[which(sparrow.map3d$SNP.Name[1] == sparrow.map3c$SNP.Name),]$cMPosition.Female)
#sparrow.map3e$cMPosition.Female <- sparrow.map3e$cMPosition.Female + (sparrow.map3d[which(sparrow.map3e$SNP.Name[1] == sparrow.map3d$SNP.Name),]$cMPosition.Female)

sparrow.map3b$cMPosition.Male <- sparrow.map3b$cMPosition.Male + (sparrow.map3a[which(sparrow.map3b$SNP.Name[1] == sparrow.map3a$SNP.Name),]$cMPosition.Male)
sparrow.map3c$cMPosition.Male <- sparrow.map3c$cMPosition.Male + (sparrow.map3b[which(sparrow.map3c$SNP.Name[1] == sparrow.map3b$SNP.Name),]$cMPosition.Male)
sparrow.map3d$cMPosition.Male <- sparrow.map3d$cMPosition.Male + (sparrow.map3c[which(sparrow.map3d$SNP.Name[1] == sparrow.map3c$SNP.Name),]$cMPosition.Male)
#sparrow.map3e$cMPosition.Male <- sparrow.map3e$cMPosition.Male + (sparrow.map3d[which(sparrow.map3e$SNP.Name[1] == sparrow.map3d$SNP.Name),]$cMPosition.Male)

#Adding order within chunks
sparrow.map3b$Order <- sparrow.map3b$Order + (sparrow.map3a[which(sparrow.map3b$SNP.Name[1] == sparrow.map3a$SNP.Name),]$Order)
sparrow.map3c$Order <- sparrow.map3c$Order + (sparrow.map3b[which(sparrow.map3c$SNP.Name[1] == sparrow.map3b$SNP.Name),]$Order)
sparrow.map3d$Order <- sparrow.map3d$Order + (sparrow.map3c[which(sparrow.map3d$SNP.Name[1] == sparrow.map3c$SNP.Name),]$Order)
#sparrow.map3e$Order <- sparrow.map3e$Order + (sparrow.map3d[which(sparrow.map3e$SNP.Name[1] == sparrow.map3d$SNP.Name),]$Order)

#Remove overlap
sparrow.map3a <- sparrow.map3a[-which(sparrow.map3a$SNP.Name %in% sparrow.map3b$SNP.Name),]
sparrow.map3b <- sparrow.map3b[-which(sparrow.map3b$SNP.Name %in% sparrow.map3c$SNP.Name),]
sparrow.map3c <- sparrow.map3c[-which(sparrow.map3c$SNP.Name %in% sparrow.map3d$SNP.Name),]
#sparrow.map3d <- sparrow.map3d[-which(sparrow.map3d$SNP.Name %in% sparrow.map3e$SNP.Name),]

#Merge
sparrow.map3a <- merge(sparrow.map3a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map3b <- merge(sparrow.map3b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map3c <- merge(sparrow.map3c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map3d <- merge(sparrow.map3d, map, by.x = "SNP.Name",by.y = "SNP.Name")
#sparrow.map3e <- merge(sparrow.map3e, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map3 <- rbind(sparrow.map3a,sparrow.map3b,sparrow.map3c,sparrow.map3d)
sparrow.map3[,c(5,6,7,8,10,11,12)] <- NULL

rm(sparrow.map3a)
rm(sparrow.map3b)
rm(sparrow.map3c)
rm(sparrow.map3d)
rm(sparrow.map3e)





#Chr 4
sparrow.map4a <- parse_map(mapfile = "crimap/chr4a.map")
sparrow.map4b <- parse_map(mapfile = "crimap/chr4b.map")
sparrow.map4c <- parse_map(mapfile = "crimap/chr4c.map")

#Combining Chr chunks

#adding cM pos within chunks
sparrow.map4b$cMPosition <- sparrow.map4b$cMPosition + (sparrow.map4a[which(sparrow.map4b$SNP.Name[1] == sparrow.map4a$SNP.Name),]$cMPosition)
sparrow.map4c$cMPosition <- sparrow.map4c$cMPosition + (sparrow.map4b[which(sparrow.map4c$SNP.Name[1] == sparrow.map4b$SNP.Name),]$cMPosition)

#adding M&F cM pos within Chunks
sparrow.map4b$cMPosition.Female <- sparrow.map4b$cMPosition.Female + (sparrow.map4a[which(sparrow.map4b$SNP.Name[1] == sparrow.map4a$SNP.Name),]$cMPosition.Female)
sparrow.map4c$cMPosition.Female <- sparrow.map4c$cMPosition.Female + (sparrow.map4b[which(sparrow.map4c$SNP.Name[1] == sparrow.map4b$SNP.Name),]$cMPosition.Female)
sparrow.map4b$cMPosition.Male <- sparrow.map4b$cMPosition.Male + (sparrow.map4a[which(sparrow.map4b$SNP.Name[1] == sparrow.map4a$SNP.Name),]$cMPosition.Male)
sparrow.map4c$cMPosition.Male <- sparrow.map4c$cMPosition.Male + (sparrow.map4b[which(sparrow.map4c$SNP.Name[1] == sparrow.map4b$SNP.Name),]$cMPosition.Male)
#Adding order within chunks
sparrow.map4b$Order <- sparrow.map4b$Order + (sparrow.map4a[which(sparrow.map4b$SNP.Name[1] == sparrow.map4a$SNP.Name),]$Order)
sparrow.map4c$Order <- sparrow.map4c$Order + (sparrow.map4b[which(sparrow.map4c$SNP.Name[1] == sparrow.map4b$SNP.Name),]$Order)

#Remove overlap
sparrow.map4a <- sparrow.map4a[-which(sparrow.map4a$SNP.Name %in% sparrow.map4b$SNP.Name),]
sparrow.map4b <- sparrow.map4b[-which(sparrow.map4b$SNP.Name %in% sparrow.map4c$SNP.Name),]

#Merge
sparrow.map4a <- merge(sparrow.map4a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map4b <- merge(sparrow.map4b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map4c <- merge(sparrow.map4c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map4 <- rbind(sparrow.map4a,sparrow.map4b,sparrow.map4c)
sparrow.map4[,c(5,6,7,8,10,11,12)] <- NULL

rm(sparrow.map4a)
rm(sparrow.map4b)
rm(sparrow.map4c)


#Chr 5
sparrow.map5a <- parse_map(mapfile = "crimap/chr5a.map")
sparrow.map5b <- parse_map(mapfile = "crimap/chr5b.map")
#sparrow.map5c <- parse_map(mapfile = "crimap/chr5c.map")

#Combining Chr chunks

#adding cM pos within chunks
sparrow.map5b$cMPosition <- sparrow.map5b$cMPosition + (sparrow.map5a[which(sparrow.map5b$SNP.Name[1] == sparrow.map5a$SNP.Name),]$cMPosition)
#sparrow.map5c$cMPosition <- sparrow.map5c$cMPosition + (sparrow.map5b[which(sparrow.map5c$SNP.Name[1] == sparrow.map5b$SNP.Name),]$cMPosition)

#adding M&F cM pos within Chunks
sparrow.map5b$cMPosition.Female <- sparrow.map5b$cMPosition.Female + (sparrow.map5a[which(sparrow.map5b$SNP.Name[1] == sparrow.map5a$SNP.Name),]$cMPosition.Female)
sparrow.map5c$cMPosition.Female <- sparrow.map5c$cMPosition.Female + (sparrow.map5b[which(sparrow.map5c$SNP.Name[1] == sparrow.map5b$SNP.Name),]$cMPosition.Female)
sparrow.map5b$cMPosition.Male <- sparrow.map5b$cMPosition.Male + (sparrow.map5a[which(sparrow.map5b$SNP.Name[1] == sparrow.map5a$SNP.Name),]$cMPosition.Male)
#sparrow.map5c$cMPosition.Male <- sparrow.map5c$cMPosition.Male + (sparrow.map5b[which(sparrow.map5c$SNP.Name[1] == sparrow.map5b$SNP.Name),]$cMPosition.Male)
#Adding order within chunks
sparrow.map5b$Order <- sparrow.map5b$Order + (sparrow.map5a[which(sparrow.map5b$SNP.Name[1] == sparrow.map5a$SNP.Name),]$Order)
#sparrow.map5c$Order <- sparrow.map5c$Order + (sparrow.map5b[which(sparrow.map5c$SNP.Name[1] == sparrow.map5b$SNP.Name),]$Order)

#Remove overlap
sparrow.map5a <- sparrow.map5a[-which(sparrow.map5a$SNP.Name %in% sparrow.map5b$SNP.Name),]
#sparrow.map5b <- sparrow.map5b[-which(sparrow.map5b$SNP.Name %in% sparrow.map5c$SNP.Name),]

#Merge
sparrow.map5a <- merge(sparrow.map5a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map5b <- merge(sparrow.map5b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map5c <- merge(sparrow.map5c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map5 <- rbind(sparrow.map5a,sparrow.map5b)
sparrow.map5[,c(5,6,7,8,10,11,12)] <- NULL

rm(sparrow.map5a)
rm(sparrow.map5b)
#rm(sparrow.map5c)


#chr6
sparrow.map6a <- parse_map(mapfile = "crimap/chr6a.map")
sparrow.map6b <- parse_map(mapfile = "crimap/chr6b.map")
#Combining Chr chunks
#adding cM pos within chunks
sparrow.map6b$cMPosition <- sparrow.map6b$cMPosition + (sparrow.map6a[which(sparrow.map6b$SNP.Name[1] == sparrow.map6a$SNP.Name),]$cMPosition)
#adding M&F cM pos within Chunks
sparrow.map6b$cMPosition.Female <- sparrow.map6b$cMPosition.Female + (sparrow.map6a[which(sparrow.map6b$SNP.Name[1] == sparrow.map6a$SNP.Name),]$cMPosition.Female)
sparrow.map6b$cMPosition.Male <- sparrow.map6b$cMPosition.Male + (sparrow.map6a[which(sparrow.map6b$SNP.Name[1] == sparrow.map6a$SNP.Name),]$cMPosition.Male)
#Adding order within chunks
sparrow.map6b$Order <- sparrow.map6b$Order + (sparrow.map6a[which(sparrow.map6b$SNP.Name[1] == sparrow.map6a$SNP.Name),]$Order)
#Remove overlap
sparrow.map6a <- sparrow.map6a[-which(sparrow.map6a$SNP.Name %in% sparrow.map6b$SNP.Name),]
#Merge
sparrow.map6a <- merge(sparrow.map6a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map6b <- merge(sparrow.map6b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map6 <- rbind(sparrow.map6a,sparrow.map6b)
sparrow.map6[,c(5,6,7,8,10,11,12)] <- NULL
rm(sparrow.map6a)
rm(sparrow.map6b)



#chr 7
sparrow.map7a <- parse_map(mapfile = "crimap/chr7a.map")
sparrow.map7b <- parse_map(mapfile = "crimap/chr7b.map")
#Combining Chr chunks
#adding cM pos within chunks
sparrow.map7b$cMPosition <- sparrow.map7b$cMPosition + (sparrow.map7a[which(sparrow.map7b$SNP.Name[1] == sparrow.map7a$SNP.Name),]$cMPosition)
#adding M&F cM pos within Chunks
sparrow.map7b$cMPosition.Female <- sparrow.map7b$cMPosition.Female + (sparrow.map7a[which(sparrow.map7b$SNP.Name[1] == sparrow.map7a$SNP.Name),]$cMPosition.Female)
sparrow.map7b$cMPosition.Male <- sparrow.map7b$cMPosition.Male + (sparrow.map7a[which(sparrow.map7b$SNP.Name[1] == sparrow.map7a$SNP.Name),]$cMPosition.Male)
#Adding order within chunks
sparrow.map7b$Order <- sparrow.map7b$Order + (sparrow.map7a[which(sparrow.map7b$SNP.Name[1] == sparrow.map7a$SNP.Name),]$Order)
#Remove overlap
sparrow.map7a <- sparrow.map7a[-which(sparrow.map7a$SNP.Name %in% sparrow.map7b$SNP.Name),]
#Merge
sparrow.map7a <- merge(sparrow.map7a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map7b <- merge(sparrow.map7b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map7 <- rbind(sparrow.map7a,sparrow.map7b)
sparrow.map7[,c(5,6,7,8,10,11,12)] <- NULL
rm(sparrow.map7a)
rm(sparrow.map7b)


#Chr 8
sparrow.map8a <- parse_map(mapfile = "crimap/chr8a.map")
sparrow.map8b <- parse_map(mapfile = "crimap/chr8b.map")
#Combining Chr chunks
#adding cM pos within chunks
sparrow.map8b$cMPosition <- sparrow.map8b$cMPosition + (sparrow.map8a[which(sparrow.map8b$SNP.Name[1] == sparrow.map8a$SNP.Name),]$cMPosition)
#adding M&F cM pos within Chunks
sparrow.map8b$cMPosition.Female <- sparrow.map8b$cMPosition.Female + (sparrow.map8a[which(sparrow.map8b$SNP.Name[1] == sparrow.map8a$SNP.Name),]$cMPosition.Female)
sparrow.map8b$cMPosition.Male <- sparrow.map8b$cMPosition.Male + (sparrow.map8a[which(sparrow.map8b$SNP.Name[1] == sparrow.map8a$SNP.Name),]$cMPosition.Male)
#Adding order within chunks
sparrow.map8b$Order <- sparrow.map8b$Order + (sparrow.map8a[which(sparrow.map8b$SNP.Name[1] == sparrow.map8a$SNP.Name),]$Order)
#Remove overlap
sparrow.map8a <- sparrow.map8a[-which(sparrow.map8a$SNP.Name %in% sparrow.map8b$SNP.Name),]
#Merge
sparrow.map8a <- merge(sparrow.map8a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map8b <- merge(sparrow.map8b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map8 <- rbind(sparrow.map8a,sparrow.map8b)
sparrow.map8[,c(5,6,7,8,10,11,12)] <- NULL
rm(sparrow.map8a)
rm(sparrow.map8b)


#Chr 29
sparrow.map29a <- parse_map(mapfile = "crimap/chr29a.map")
sparrow.map29b <- parse_map(mapfile = "crimap/chr29b.map")
sparrow.map29c <- parse_map(mapfile = "crimap/chr29c.map")

#Combining Chr chunks
#adding cM pos within chunks
sparrow.map29b$cMPosition <- sparrow.map29b$cMPosition + (sparrow.map29a[which(sparrow.map29b$SNP.Name[1] == sparrow.map29a$SNP.Name),]$cMPosition)
sparrow.map29c$cMPosition <- sparrow.map29c$cMPosition + (sparrow.map29b[which(sparrow.map29c$SNP.Name[1] == sparrow.map29b$SNP.Name),]$cMPosition)
#adding M&F cM pos within Chunks
sparrow.map29b$cMPosition.Female <- sparrow.map29b$cMPosition.Female + (sparrow.map29a[which(sparrow.map29b$SNP.Name[1] == sparrow.map29a$SNP.Name),]$cMPosition.Female)
sparrow.map29c$cMPosition.Female <- sparrow.map29c$cMPosition.Female + (sparrow.map29b[which(sparrow.map29c$SNP.Name[1] == sparrow.map29b$SNP.Name),]$cMPosition.Female)
sparrow.map29b$cMPosition.Male <- sparrow.map29b$cMPosition.Male + (sparrow.map29a[which(sparrow.map29b$SNP.Name[1] == sparrow.map29a$SNP.Name),]$cMPosition.Male)
sparrow.map29c$cMPosition.Male <- sparrow.map29c$cMPosition.Male + (sparrow.map29b[which(sparrow.map29c$SNP.Name[1] == sparrow.map29b$SNP.Name),]$cMPosition.Male)
#Adding order within chunks
sparrow.map29b$Order <- sparrow.map29b$Order + (sparrow.map29a[which(sparrow.map29b$SNP.Name[1] == sparrow.map29a$SNP.Name),]$Order)
sparrow.map29c$Order <- sparrow.map29c$Order + (sparrow.map29b[which(sparrow.map29c$SNP.Name[1] == sparrow.map29b$SNP.Name),]$Order)
#Remove overlap
sparrow.map29a <- sparrow.map29a[-which(sparrow.map29a$SNP.Name %in% sparrow.map29b$SNP.Name),]
sparrow.map29b <- sparrow.map29b[-which(sparrow.map29b$SNP.Name %in% sparrow.map29c$SNP.Name),]
#Merge
sparrow.map29a <- merge(sparrow.map29a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map29b <- merge(sparrow.map29b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map29c <- merge(sparrow.map29c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map29 <- rbind(sparrow.map29a,sparrow.map29b,sparrow.map29c)
sparrow.map29[,c(5,6,7,8,10,11,12)] <- NULL
rm(sparrow.map29a)
rm(sparrow.map29b)
rm(sparrow.map29c)


#Ch 6
sparrow.map6 <- parse_map(mapfile = "crimap/chr6a.map")
sparrow.map6 <- merge(sparrow.map6, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map6[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 9
sparrow.map9 <- parse_map(mapfile = "crimap/chr9a.map")
sparrow.map9 <- merge(sparrow.map9, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map9[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 10
sparrow.map10 <- parse_map(mapfile = "crimap/chr10a.map")
sparrow.map10 <- merge(sparrow.map10, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map10[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 11
sparrow.map11 <- parse_map(mapfile = "crimap/chr11a.map")
sparrow.map11 <- merge(sparrow.map11, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map11[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 12
sparrow.map12 <- parse_map(mapfile = "crimap/chr12a.map")
sparrow.map12 <- merge(sparrow.map12, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map12[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 13
sparrow.map13 <- parse_map(mapfile = "crimap/chr13a.map")
sparrow.map13 <- merge(sparrow.map13, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map13[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 14
sparrow.map14 <- parse_map(mapfile = "crimap/chr14a.map")
sparrow.map14 <- merge(sparrow.map14, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map14[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 15
sparrow.map15 <- parse_map(mapfile = "crimap/chr15a.map")
sparrow.map15 <- merge(sparrow.map15, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map15[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 17
sparrow.map17 <- parse_map(mapfile = "crimap/chr17a.map")
sparrow.map17 <- merge(sparrow.map17, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map17[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 18
sparrow.map18 <- parse_map(mapfile = "crimap/chr18a.map")
sparrow.map18 <- merge(sparrow.map18, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map18[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 19
sparrow.map19 <- parse_map(mapfile = "crimap/chr19a.map")
sparrow.map19 <- merge(sparrow.map19, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map19[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 20
sparrow.map20 <- parse_map(mapfile = "crimap/chr20a.map")
sparrow.map20 <- merge(sparrow.map20, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map20[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 21
sparrow.map21 <- parse_map(mapfile = "crimap/chr21a.map")
sparrow.map21 <- merge(sparrow.map21, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map21[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 22
sparrow.map22 <- parse_map(mapfile = "crimap/chr22a.map")
sparrow.map22 <- merge(sparrow.map22, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map22[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 23
sparrow.map23 <- parse_map(mapfile = "crimap/chr23a.map")
sparrow.map23 <- merge(sparrow.map23, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map23[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 24
sparrow.map24 <- parse_map(mapfile = "crimap/chr24a.map")
sparrow.map24 <- merge(sparrow.map24, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map24[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 26
sparrow.map26 <- parse_map(mapfile = "crimap/chr26a.map")
sparrow.map26 <- merge(sparrow.map26, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map26[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 27
sparrow.map27 <- parse_map(mapfile = "crimap/chr27a.map")
sparrow.map27 <- merge(sparrow.map27, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map27[,c(5,6,7,8,10,11,12)] <- NULL

#Chr 28
sparrow.map28 <- parse_map(mapfile = "crimap/chr28a.map")
sparrow.map28 <- merge(sparrow.map28, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map28[,c(5,6,7,8,10,11,12)] <- NULL

sparrow.map1$Chr <- "Chr 1"
sparrow.map2$Chr <- "Chr 2"
sparrow.map3$Chr <- "Chr 3"
sparrow.map4$Chr <- "Chr 4"
sparrow.map5$Chr <- "Chr 5"
sparrow.map6$Chr <- "Chr 6"
sparrow.map7$Chr <- "Chr 7"
sparrow.map8$Chr <- "Chr 8"
sparrow.map9$Chr <- "Chr 9"
sparrow.map10$Chr <- "Chr 10"
sparrow.map11$Chr <- "Chr 11"
sparrow.map12$Chr <- "Chr 12"
sparrow.map13$Chr <- "Chr 13"
sparrow.map14$Chr <- "Chr 14"
sparrow.map15$Chr <- "Chr 15"
sparrow.map17$Chr <- "Chr 17"
sparrow.map18$Chr <- "Chr 18"
sparrow.map19$Chr <- "Chr 19"
sparrow.map20$Chr <- "Chr 20"
sparrow.map21$Chr <- "Chr 21"
sparrow.map22$Chr <- "Chr 22"
sparrow.map23$Chr <- "Chr 23" 
sparrow.map24$Chr <- "Chr 24"
sparrow.map26$Chr <- "Chr 26"
sparrow.map27$Chr <- "Chr 27"
sparrow.map28$Chr <- "Chr 28"
sparrow.map29$Chr <- "Chr 29"
sparrow.map.all <- rbind(sparrow.map1,sparrow.map2,sparrow.map3,sparrow.map4,sparrow.map5,sparrow.map6,sparrow.map7,sparrow.map8,
                     sparrow.map9,sparrow.map10,sparrow.map11,sparrow.map12,sparrow.map13,sparrow.map14,sparrow.map15,sparrow.map17,
                     sparrow.map18,sparrow.map19,sparrow.map20,sparrow.map21,sparrow.map22,sparrow.map23,sparrow.map24,sparrow.map26,
                     sparrow.map27,sparrow.map28,sparrow.map29)
sparrow.map.all$Chr <- factor(sparrow.map.all$Chr, levels=c("Chr 1","Chr 2","Chr 3","Chr 4","Chr 5","Chr 6","Chr 7","Chr 8","Chr 9",
                                                               "Chr 10","Chr 11","Chr 12","Chr 13","Chr 14","Chr 15","Chr 17","Chr 18",
                                                               "Chr 19","Chr 20","Chr 21","Chr 22","Chr 23","Chr 24","Chr 26","Chr 27",
                                                               "Chr 28","Chr 29"))
p <-ggplot(data = sparrow.map.all, aes(x = BP/1000000, y = cMPosition)) +
  geom_point(size = 1) +
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")
p + facet_wrap(~Chr, scales = "free")

p <-ggplot(data = sparrow.map.all, aes(x = BP/1000000, y = cMPosition.Female)) +
  geom_point(size = 1) +
  labs(y= "Female Linkage Map Length (cM)", x ="Female Chromosome Length (Mb)")
p + facet_wrap(~Chr, scales = "free")

p <-ggplot(data = sparrow.map.all, aes(x = BP/1000000, y = cMPosition.Male)) +
  geom_point(size = 1) +
  labs(y= "Male Linkage Map Length (cM)", x ="Male Chromosome Length (Mb)")
p + facet_wrap(~Chr, scales = "free")

p <-ggplot(data = sparrow.map.all) +
  geom_point(aes(x = BP/1000000, y = cMPosition.Male), color = "blue", size = 1) +
  geom_point(aes(x = BP/1000000, y = cMPosition.Female), color = "red", size = 1) +
  labs(y= "Linkage Map Length (cM)", x ="Chromosome Length (Mb)")
p + facet_wrap(~Chr, scales = "free")





#Plotting cM x SNP Order & cM x BP
library(ggplot2)

ggplot() +
  geom_point(data = sparrow.map1, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 1: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map2, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 2: Chr length x Map Length')


ggplot() +
  geom_point(data = sparrow.map3, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 3: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map4, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 4: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map5, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 5: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map6, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 6: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map7, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 7: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map8, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 8: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map29, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 29: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map1, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 1: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map2, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 2: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map3, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 3: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map4, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 4: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map5, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 5: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map6, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 6: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map7, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 7: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map8, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 8: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map29, aes(x = BP/1000000, y = cMPosition.Female)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Female') +
  ggtitle('Chromosome 29: Chr length x Map Length')


ggplot() +
  geom_point(data = sparrow.map1, aes(x = BP/1000000, y = cMPosition.Male)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Male') +
  ggtitle('Chromosome 1: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map2, aes(x = BP/1000000, y = cMPosition.Male)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Male') +
  ggtitle('Chromosome 2: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map3, aes(x = BP/1000000, y = cMPosition.Male)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Male') +
  ggtitle('Chromosome 3: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map4, aes(x = BP/1000000, y = cMPosition.Male)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Male') +
  ggtitle('Chromosome 4: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map5, aes(x = BP/1000000, y = cMPosition.Male)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Male') +
  ggtitle('Chromosome 5: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map6, aes(x = BP/1000000, y = cMPosition.Male)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Male') +
  ggtitle('Chromosome 6: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map7, aes(x = BP/1000000, y = cMPosition.Male)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Male') +
  ggtitle('Chromosome 7: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map8, aes(x = BP/1000000, y = cMPosition.Male)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Male') +
  ggtitle('Chromosome 8: Chr length x Map Length')

ggplot() +
  geom_point(data = sparrow.map29, aes(x = BP/1000000, y = cMPosition.Male)) +
  xlab('Chromosome length (Mb)') +
  ylab('cM Position.Male') +
  ggtitle('Chromosome 29: Chr length x Map Length')


#===================================
#===
#===
#===
#Graphing Chr length x Map Length
ChrLxMapL <- vector()
ChrLxMapL$Chr <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,26,27,28,29)
ChrLxMapL <- as.data.frame(ChrLxMapL)
ChrLxMapL$ChrL <- c(max(sparrow.map1$BP),max(sparrow.map2$BP),max(sparrow.map3$BP),max(sparrow.map4$BP),
               max(sparrow.map5$BP),max(sparrow.map6$BP),max(sparrow.map7$BP),max(sparrow.map8$BP),
               max(sparrow.map9$BP),max(sparrow.map10$BP),max(sparrow.map11$BP),max(sparrow.map12$BP),
               max(sparrow.map13$BP),max(sparrow.map14$BP),max(sparrow.map15$BP),max(sparrow.map17$BP),
               max(sparrow.map18$BP),max(sparrow.map19$BP),max(sparrow.map20$BP),max(sparrow.map21$BP),
               max(sparrow.map22$BP),max(sparrow.map23$BP),max(sparrow.map24$BP),max(sparrow.map26$BP),
               max(sparrow.map27$BP),max(sparrow.map28$BP),max(sparrow.map29$BP))
ChrLxMapL$ChrL <- ChrLxMapL$ChrL/1000000
ChrLxMapL$MapL <- c(max(sparrow.map1$cMPosition),max(sparrow.map2$cMPosition),max(sparrow.map3$cMPosition),max(sparrow.map4$cMPosition),
                    max(sparrow.map5$cMPosition),max(sparrow.map6$cMPosition),max(sparrow.map7$cMPosition),max(sparrow.map8$cMPosition),
                    max(sparrow.map9$cMPosition),max(sparrow.map10$cMPosition),max(sparrow.map11$cMPosition),max(sparrow.map12$cMPosition),
                    max(sparrow.map13$cMPosition),max(sparrow.map14$cMPosition),max(sparrow.map15$cMPosition),max(sparrow.map17$cMPosition),
                    max(sparrow.map18$cMPosition),max(sparrow.map19$cMPosition),max(sparrow.map20$cMPosition),max(sparrow.map21$cMPosition),
                    max(sparrow.map22$cMPosition),max(sparrow.map23$cMPosition),max(sparrow.map24$cMPosition),max(sparrow.map26$cMPosition),
                    max(sparrow.map27$cMPosition),max(sparrow.map28$cMPosition),max(sparrow.map29$cMPosition))
ChrLxMapL <- as.data.frame(ChrLxMapL)
ChrLxMapL$Chr <- as.character(ChrLxMapL$Chr)
ggplot(data = ChrLxMapL, aes(x = ChrL, y = MapL, label = Chr)) +
  geom_point() +
  xlab('Chromosome Length (Mb)') +
  ylab('Map Length (cM)') +
  ggtitle('')+
  geom_text(size = 6, hjust = -.5, vjust = -.5) +
  geom_smooth(method = "lm")+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18))
  
#Graphing Female x Male cM
FcmxMcm <- vector()
FcmxMcm$Chr <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,26,27,28,29)
FcmxMcm <- as.data.frame(FcmxMcm)
FcmxMcm$cmF <- c(max(sparrow.map1$cMPosition.Female),max(sparrow.map2$cMPosition.Female),max(sparrow.map3$cMPosition.Female),max(sparrow.map4$cMPosition.Female),
                    max(sparrow.map5$cMPosition.Female),max(sparrow.map6$cMPosition.Female),max(sparrow.map7$cMPosition.Female),max(sparrow.map8$cMPosition.Female),
                    max(sparrow.map9$cMPosition.Female),max(sparrow.map10$cMPosition.Female),max(sparrow.map11$cMPosition.Female),max(sparrow.map12$cMPosition.Female),
                    max(sparrow.map13$cMPosition.Female),max(sparrow.map14$cMPosition.Female),max(sparrow.map15$cMPosition.Female),max(sparrow.map17$cMPosition.Female),
                    max(sparrow.map18$cMPosition.Female),max(sparrow.map19$cMPosition.Female),max(sparrow.map20$cMPosition.Female),max(sparrow.map21$cMPosition.Female),
                    max(sparrow.map22$cMPosition.Female),max(sparrow.map23$cMPosition.Female),max(sparrow.map24$cMPosition.Female),max(sparrow.map26$cMPosition.Female),
                    max(sparrow.map27$cMPosition.Female),max(sparrow.map28$cMPosition.Female),max(sparrow.map29$cMPosition.Female))
FcmxMcm$cmM <- c(max(sparrow.map1$cMPosition.Male),max(sparrow.map2$cMPosition.Male),max(sparrow.map3$cMPosition.Male),max(sparrow.map4$cMPosition.Male),
                 max(sparrow.map5$cMPosition.Male),max(sparrow.map6$cMPosition.Male),max(sparrow.map7$cMPosition.Male),max(sparrow.map8$cMPosition.Male),
                 max(sparrow.map9$cMPosition.Male),max(sparrow.map10$cMPosition.Male),max(sparrow.map11$cMPosition.Male),max(sparrow.map12$cMPosition.Male),
                 max(sparrow.map13$cMPosition.Male),max(sparrow.map14$cMPosition.Male),max(sparrow.map15$cMPosition.Male),max(sparrow.map17$cMPosition.Male),
                 max(sparrow.map18$cMPosition.Male),max(sparrow.map19$cMPosition.Male),max(sparrow.map20$cMPosition.Male),max(sparrow.map21$cMPosition.Male),
                 max(sparrow.map22$cMPosition.Male),max(sparrow.map23$cMPosition.Male),max(sparrow.map24$cMPosition.Male),max(sparrow.map26$cMPosition.Male),
                 max(sparrow.map27$cMPosition.Male),max(sparrow.map28$cMPosition.Male),max(sparrow.map29$cMPosition.Male))
ggplot(data = FcmxMcm, aes(x = cmM, y = cmF, label = Chr)) +
  geom_point() +
  xlab('Male Linakge Map Length (cM)') +
  ylab('Female Linkage Map Length (cM)') +
  ggtitle('')+
  geom_text(size = 8, hjust = -.5, vjust = -.5) +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = .75)+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18))
