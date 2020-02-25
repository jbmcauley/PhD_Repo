setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/")

read.table()


#1.       Make bed for the plink file (remember to include --autosome-num 31)

#    plink --file Pdo_200k_n3960_21032017 --maf 0.05 --autosome-num 31 --make-bed --out newplink

#2.       Read .fam into R
  famfile <- read.table("newplink.fam")
  names(famfile) <- c("Fam", "ID","Father", "Mother", "Sex", "Phenotype")
  library(reshape2)
  
  fatherids <-  melt(famfile, id.vars = "ID", measure.vars = "Father")
  fatherids <- fatherids$value
  fatheridslist <- list()
  counter <- 0
  for (i in 1:length(fatherids)) {
    if(fatherids[i] %in% new.ids$OldID){
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
    if(motherids[i] %in% new.ids$OldID){
      counter <- counter + 1
      motheridslist[[counter]] <- new.ids$NewID[which(new.ids$OldID == motherids[i])]
    }else{
      counter <- counter + 1
      motheridslist[[counter]] <- 0
    }  
  }
  motheridslist <- do.call(rbind, motheridslist)
  motheridslist <- as.vector(motheridslist)
  
  fullids$ID <- as.character(fullids$ID)
  idslist <- c(fullids$ID, fullids$value)
  idslist <- unique(idslist)
  idslist <- idslist[-3961]
  which(idslist == "0")

  pedigree.new <- as.data.frame(1:3960)
  names(pedigree.new)[1] <- "FamID"
  pedigree.new$FamID <- 1
  pedigree.new$withinFamID <- new.ids$NewID[1:3960]
  pedigree.new$paternalID <- fatheridslist
  pedigree.new$maternalID <- motheridslist
  
  write.table(pedigree.new, file = "updateParents.txt",sep=" ",row.names=FALSE, col.names = FALSE)
  #Can now update parental ids
#a.       Recode the IDs to numbers
  ids <- 1:4226
  new.ids <- data.frame(1, idslist, 1, ids)
  new.ids$OldID <- as.character(new.ids$OldID)
  names(new.ids) <- c("OldFam", "OldID", "NewFam", "NewID")
  head(new.ids)
  
 
  
#b. Make a phenotype file for Genabel (id, sex, knownsex)
  rm(pheno.file)
  sex <- read.table("newplinksex.sexcheck", header = TRUE)
  knownsex <- read.table("All_sexings_Helgeland_200k_29112017.txt", header = TRUE)
  pheno.file <- data.frame(famfile$V2, sex$SNPSEX)
  names(pheno.file) <- c("ID","sex")
  pheno.file <- as.data.frame(pheno.file)
  pheno.file$ID <- new.ids$NewID 
  pheno.file$sex[pheno.file$sex == 2] <- "m"
  pheno.file$sex[pheno.file$sex == 1] <- "f"
  pheno.file$sex[pheno.file$sex == 0] <- "u"
  pheno.file$knownsex <- pheno.file$sex
  pheno.file$sex[pheno.file$sex == "m"] <- 1
  pheno.file$sex[pheno.file$sex == "f"] <- 0
  pheno.file$sex[pheno.file$sex == "u"] <- 1
  #Unknown sex coded to males!
  #Alternatively can reassign based on F value?
  #pheno.file$sex[which(sex$F < .5 & sex$SNPSEX == 0)] <- 1 #recode to m
  #pheno.file$sex[which(sex$F >= .5 & sex$SNPSEX == 0)] <- 0  #recode to f
  
#c.       Make the file for recoding IDs in PLINK (4 cols – OldFamily, Old ID, NewFamily, NewID) 

  write.table(new.ids, file = "updateIDs.txt",sep=" ",row.names=FALSE, col.names = FALSE)
  
#3.       Recode the IDs in the PLINK file with –update-ids into a new PLINk file

#4.       Make the .genabelmap file for making the genabel object
#a.       Read in PLINK map, remove the cM column, save without headers

#5.       Make the GenABEL file
      convert.snp.ped(pedfile = "newfile.ped", 
                mapfile = "newfile.genabelmap",
                outfile = "newfile.gen",
                strand = "u", bcast = 1000, traits = 1, mapHasHeaderLine = F)

#6.  
           sparrowabel <- load.gwaa.data(phe = "newfile.pheno", 
                                  gen = "newfile.gen")
