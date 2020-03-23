setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/")
library(GenABEL)
library(crimaptools)
read.table()


#1.       Make bed for the plink file (remember to include --autosome-num 31)

#    plink --file Pdo_200k_n3960_21032017 --maf 0.05 --autosome-num 31 --make-bed --out newplink

#2.       Read .fam into R
  famfile <- read.table("newplink.fam")
  names(famfile) <- c("Fam", "ID","Father", "Mother", "Sex", "Phenotype")
  
  ids <- 1:3960
  new.ids <- data.frame(1, famfile$ID, 1, ids)
  names(new.ids) <- c("OldFam", "OldID", "NewFam", "NewID", rownames = FALSE)
  new.ids$OldID <- as.character(new.ids$OldID)
  row.names(new.ids) <- NULL
  write.table(new.ids, file = "ID-Recode-Key.txt", row.names = FALSE)
  head(new.ids)
  
  
  library(reshape2)
  
  fatherids <-  melt(famfile, id.vars = "ID", measure.vars = "Father")
  fatherids <- fatherids$value
  fatheridslist <- list()
  counter <- 0
  for (i in 1:length(fatherids)) {
    if(fatherids[i] %in% famfile$ID){
    counter <- counter + 1
    fatheridslist[[counter]] <- new.ids$NewID[which(new.ids$OldID == fatherids[i])]
  }else{
    counter <- counter + 1
    fatheridslist[[counter]] <- 0
  }  
  }
  fatheridslist <- do.call(rbind, fatheridslist)
  fatheridslist <- as.vector(fatheridslist)
  
  motherids <-  melt(famfile, id.vars = "ID", measure.vars = "Mother")
  motherids <- motherids$value
  motheridslist <- list()
  counter <- 0
  for (i in 1:length(motherids)) {
    if(motherids[i] %in% famfile$ID){
      counter <- counter + 1
      motheridslist[[counter]] <- new.ids$NewID[which(new.ids$OldID == motherids[i])]
    }else{
      counter <- counter + 1
      motheridslist[[counter]] <- 0
    }  
  }
  motheridslist <- do.call(rbind, motheridslist)
  motheridslist <- as.vector(motheridslist)
  
  fullids <- data.frame()
  fullids$ID <- as.character(fullids$ID)
  idslist <- c(fullids$ID, fullids$value)
  idslist <- unique(idslist)
  idslist <- idslist[-3961]
  which(idslist == "0")

  pedigree.new <- as.data.frame(1:3960)
  names(pedigree.new)[1] <- "FamID"
  pedigree.new$FamID <- 1
  pedigree.new$withinFamID <- 1:3960
  pedigree.new$paternalID <- fatheridslist
  pedigree.new$maternalID <- motheridslist
  
  write.table(pedigree.new, file = "updateParents.txt",sep=" ",row.names=FALSE, col.names = FALSE)
  #Can now update parental ids
#a.       Recode the IDs to numbers

  "plink --file Pdo_200k_n3960_21032017 --autosome-num 31 --maf 0.05 --check-sex --make-bed --out newplinksex"
#b. Make a phenotype file for Genabel (id, sex, knownsex)
  rm(pheno.file)
  sex <- read.table("newplinksex.sexcheck", header = TRUE)
  #knownsex <- read.table("All_sexings_Helgeland_200k_29112017.txt", header = TRUE)
  pheno.file <- data.frame(1, sex$SNPSEX)
  names(pheno.file) <- c("id","sex")
  pheno.file <- as.data.frame(pheno.file)
  pheno.file$id <- 1:3960 
  #Sexs assinged by plink sex check are incorrect because of Birds having alternate heterogametic sex
  #Unknowns coded to males; ind. 136,302,340,348,417,517,537,596,627,646,663,664,706,1000,1704,1902,1940,1993,2143,2292,2380,2434,2486,2493,2538,2654,2723,2816,2978,3158,3272,3314,3388,3538,3911 mothers coded to female to match pedigree, SNP assigned male
  #ind. 537,563,764,2030,3779 fathers coded to male to match pedigree
  pheno.file$sex[pheno.file$sex == 2] <- "m"
  pheno.file$sex[pheno.file$sex == 1] <- "f"
  pheno.file$sex[pheno.file$sex == 0] <- "u"
  pheno.file$sex[pheno.file$sex == "m"] <- 1
  pheno.file$sex[pheno.file$sex == "f"] <- 2
  pheno.file$sex[pheno.file$sex == "u"] <- 1
  write.table(pheno.file,"newfile.pheno",col.names = TRUE, row.names = FALSE)
  pheno.file[2] <- NULL
  pheno.file$sex <- as.numeric(pheno.file$sex)
  write.table(pheno.file, file = "updateSex.txt", sep = " ", row.names = FALSE, col.names = FALSE)
  
  pheno.file$sex <- as.numeric(pheno.file$sex)
  
  pheno.file[2:3] <- pheno.file[1:2]
  pheno.file[1] <- 1
  names(pheno.file) <-  c("FID", "IID", "sex")

  #Unknown sex coded to males!
  #Alternatively can reassign based on F value?
  #pheno.file$sex[which(sex$F < .5 & sex$SNPSEX == 0)] <- 1 #recode to m
  #pheno.file$sex[which(sex$F >= .5 & sex$SNPSEX == 0)] <- 0  #recode to f
  
#c.       Make the file for recoding IDs in PLINK (4 cols – OldFamily, Old ID, NewFamily, NewID) 

  write.table(new.ids, file = "updateIDs.txt",sep=" ",row.names=FALSE, col.names = FALSE)
  
#3.       Recode the IDs in the PLINK file with –update-ids into a new PLINk file

  system("plink --help")
  
#4.       Make the .genabelmap file for making the genabel object
  #Make Binary Files
  system("plink --file Pdo_200k_n3960_21032017 --autosome-num 31 --maf 0.05 --make-bed --out workfiles")
  
  #Update ids
  system("plink --bfile workfiles --update-ids updateIDs.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_ids")
  
  
  system("plink --bfile workfiles_ids --update-parents updateParents.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_parids")
  
  
  system("plink --bfile workfiles_parids --update-sex updateSex.txt --autosome-num 31 --maf 0.05 --make-bed --out finalbinaries")
  
  
  system("plink --bfile finalbinaries --recode 12 --autosome-num 31 --maf 0.05 --out test")

  #a.       Read in PLINK map, remove the cM column, save without headers
  
  mapfile <- read.table("pedfile.map")
  mapfile <- mapfile[,-3]
  write.table(mapfile, file = "mapfile.txt", sep=" ",row.names = FALSE,col.names = FALSE)
  
#5.       Make the GenABEL file
      convert.snp.ped(pedfile = "test.ped", 
                mapfile = "mapfile.txt",
                outfile = "newfile.gen",
                strand = "u", bcast = 1000, traits = 1, mapHasHeaderLine = FALSE)

     ### Error in convert.snp.ped(pedfile = "Pdo_200k_n3960_21032017.ped", mapfile = "mapfile.txt",  : 
     ###                           coding '43' for SNP not recognised !
      
#6.  
      
      
           sparrowabel <- load.gwaa.data(phe = "newfile.pheno", 
                                  gen = "newfile.gen")
#After mendel errors:
           system("plink --bfile finalbinaries --set-me-missing --autosome-num 31 --maf 0.05 --out test_ME-fixed")
           mapfile <- read.table("test_ME-fixed.map")
           mapfile <- mapfile[,-3]
           write.table(mapfile, file = "mapfile_ME-fixed.txt", sep=" ",row.names = FALSE,col.names = FALSE)
           convert.snp.ped(pedfile = "test_ME-fixed.ped", 
                           mapfile = "mapfile_ME-fixed.txt",
                           outfile = "newfile-ME-fixed.gen",
                           strand = "u", bcast = 1000, traits = 1, mapHasHeaderLine = FALSE)
           sparrowabel <- load.gwaa.data(phe = "newfile.pheno", 
                                         gen = "newfile-ME-fixed.gen")

           data1 <- sparrowabel
           #Relaxed callrate to 0.9
           qc1 <- check.marker(data1, callrate = .9, ibs.threshold = .9, maf = .01, p.level = 0)
           
           
           #NEED TO SORT OUT HOW TO PROPERLY CHECK ZZ WZ
           data1 <- data1[qc1$idok, qc1$snpok]
           data1 <- Xfix(data1)
           save(data1, file = "sparrowABEL_QC.RData")
           load("sparrowABEL_QC.RData")
           