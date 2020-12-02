
#~~ Load libraries

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(purrr)    # for reduce() function
library(GenABEL)

#~~ load the data and functions

load("3_Examine_SNP_Clusters.RData", verbose = T)

hist(fullsnp$stdevTheta_AB)

# Prune SNPs based on the SD of Theta. (NB your script was getting rid of rare
# SNPs. Best to use -which and reverse the sign.)

stdev_thresh <- .046

z <- fullsnp[-which(fullsnp$stdevTheta_AB > stdev_thresh),]
z <- z[-which(z$stdevTheta_AA > stdev_thresh),]
z <- z[-which(z$stdevTheta_BB > stdev_thresh),]

# z <- z[1:2200,] JUST FOR TESTING

# Split dataset into chunks. (NB I hit memory issues with just 4, so I had this
# solution that makes a vector z$chunk in the data frame):

chunk_size <- 1000

z$chunk <- rep(1:ceiling(nrow(z)/chunk_size), each = chunk_size, length.out = nrow(z))

table(z$chunk)

# Extract the SNP information and create a genotype table

geno_list <- list()

for(i in 1:max(z$chunk)){
  
  message(paste("Running chunk", i, "of", max(z$chunk)))
  
  # Extract SNPs into an object
  
  xsnps <- snp_extract_func(c(z$cust_id[which(z$chunk == i)]))

  #~~ Run John's script to extract the good SNPs
  
  new_x1 <- list()
  
  counter <- 0
  
  for (j in 1:length(unique(xsnps$cust_id))) {
    
    x1 <- xsnps[which(xsnps$cust_id == unique(xsnps$cust_id)[j]),] 
    aa <- x1[which(x1$AB_geno == "AA"),]
    ab <- x1[which(x1$AB_geno == "AB"),]
    bb <- x1[which(x1$AB_geno == "BB"),]
    
    counter <- counter + 1
    new_x1[[counter]] <- aa[between(aa$Theta ,mean(aa$Theta)-2*sd(aa$Theta),mean(aa$Theta)+2*sd(aa$Theta)),]
    
    counter <- counter + 1
    new_x1[[counter]] <- ab[between(ab$Theta ,mean(ab$Theta)-2*sd(ab$Theta),mean(ab$Theta)+2*sd(ab$Theta)),]
    
    counter <- counter + 1
    new_x1[[counter]] <- bb[between(bb$Theta ,mean(bb$Theta)-2*sd(bb$Theta),mean(bb$Theta)+2*sd(bb$Theta)),]
    
    rm(x1, aa, ab, bb)
  } 
  
  new_x1 <- bind_rows(new_x1)     # much faster than do.call
  
  #~~ Make into a genotype table
  
  new_x1 <- subset(new_x1, select = c(cust_id, ringnr, id, variable, AB_geno))
  str(new_x1)
  
  new_x1 <- dcast(new_x1, ringnr + id + variable ~ cust_id)
  # new_x1[1:10, 1:10]   # JUST TO EXAMINE
  
  geno_list[[i]] <- new_x1
  
  rm(counter, j, new_x1, xsnps)

}

#~~ Merge all the elements of the list

test <- geno_list %>% reduce(full_join, by = c("ringnr", "id", "variable")) 

# NEW LINES HERE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ examine duplicate IDs & select the entry with the highest genotyping success

test <- add_count(test, ringnr, name = "Count")

test$PropNAs <- apply(test[,4:ncol(test)],                                      # select columns 4 to the end
                          1,                                                    # apply the function over rows
                          function(foo) length(which(is.na(foo)))/length(foo))  # calculate % of NAs

test <- test %>%
  group_by(ringnr) %>%
  summarise(MinPropNAs = min(PropNAs)) %>%
  right_join(test)

test <- subset(test, MinPropNAs == PropNAs)

test <- subset(test, select = -c(PropNAs, MinPropNAs, variable, Count))

#~~ Make PLINK files

# Phenotype file

phenotab <- data.frame(Family = 1,
                       id = as.numeric(idnames(sparrowgen)),
                       Dad = 0,
                       Mum = 0,
                       sex = phdata(sparrowgen)$sex,
                       Trait = -9)

phenotab <- subset(phenotab, id %in% test$id)

# Genotype file - need to recode the SNPs

test[1:10, 1:10]

test <- right_join(phenotab, test)
test$ringnr <- NULL

test[test == "AA"] <- "1 1"
test[test == "AB"] <- "1 2"
test[test == "BB"] <- "2 2"
test[is.na(test)] <- "0 0"

#~~ Map file

maptab <- data.frame(Chromosome = chromosome(sparrowgen),
                     SNP.Name = snpnames(sparrowgen),
                     Position = map(sparrowgen), stringsAsFactors = F)

maptab <- subset(maptab, SNP.Name %in% names(test))

map2 <- data.frame(SNP.Name = names(test)[7:ncol(test)])
map2$Order <- 1:nrow(map2)

maptab <- left_join(maptab, map2)
maptab <- arrange(maptab, Order)

maptab$Order <- NULL

#~~ Write to file

phenotab <- subset(phenotab, select = c(id, sex))

if(!dir.exists("clean_genotype_data")) dir.create("clean_genotype_data")

write.table(maptab,   "clean_genotype_data/Clean_Map_SNPs.map", row.names = F, col.names = F, quote = F)
write.table(test,     "clean_genotype_data/Clean_Map_SNPs.ped", row.names = F, col.names = F, quote = F)
write.table(phenotab, "clean_genotype_data/Clean_Map_SNPs.phe", row.names = F, quote = F)

#~~ Convert to GenABEL


convert.snp.ped(pedfile = "clean_genotype_data/Clean_Map_SNPs.ped", 
                mapfile = "clean_genotype_data/Clean_Map_SNPs.map",
                outfile = "clean_genotype_data/Clean_Map_SNPs.gen",
                strand = "u", bcast = 1000, traits = 1, mapHasHeaderLine = F)   #added traits = 1

cleanspar <- load.gwaa.data(phenofile = "clean_genotype_data/Clean_Map_SNPs.phe",
                         genofile = "clean_genotype_data/Clean_Map_SNPs.gen")

nids(cleanspar)
nsnps(cleanspar)


save(cleanspar, file = "clean_genotype_data/Clean_Map_SNPs.RData")
#saveRDS(cleanspar, file = "clean_genotype_data/Clean_Map_SNPs.rds")
