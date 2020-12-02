setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/")
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/Ped_&_map_file_formatting/")
Good_SNPS <- read.table("GoodSNPlist.txt")
#Recoding the Helgeland ID's to the New IDS established.
x <- read.table("ID-Recode-Helgeland.txt", header = FALSE)
y <- read.table("updateParents_Helgeland.txt", header = FALSE)


Key_IDS <- read.table("ID-Recode-Key.txt", header = TRUE)
Key_IDS$OldID <- as.character(Key_IDS$OldID)

Helg_IDS <- as.data.frame(sparrowgen.Helgeland@phdata$id)
names(Helg_IDS) <- "OldID"
Helg_IDS$OldID <- as.character(Helg_IDS$OldID)

Helg_IDS$NewID <- Key_IDS[which(Key_IDS$OldID %in% Helg_IDS$OldID),]$NewID

#Recode Parents
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/Ped_&_map_file_formatting/")
Helg_Ped <- read.table("SNP_pedigree_Helgeland_05122017.txt", header = TRUE)
Helg_Ped <- Helg_Ped[1:3116,]
Helg_IDS$OldFather <- as.character(Helg_Ped$sire)
Helg_IDS$OldMother <- as.character(Helg_Ped$dam)
Helg_IDS[is.na(Helg_IDS)] <- 0
setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/")

#Fathers
faths <- list()
counter <- 0
for (i in 1:length(Helg_IDS$OldFather)) {
counter <- counter +1
  if(Helg_IDS$OldFather[i] != "0" ){
faths[[counter]] <- Key_IDS[which(Key_IDS$OldID %in% Helg_IDS$OldFather[i]),]$NewID}
  else{ 
  faths[[counter]] <- 0
  }
}
idx <- !(sapply(faths, length))
faths[idx] <- 0
faths <- do.call(rbind, faths)
faths <- as.data.frame(faths)
Helg_IDS$NewFather <- faths$V1

#Mothers
moths <- list()
counter <- 0
for (i in 1:length(Helg_IDS$OldMother)) {
  counter <- counter +1
  if(Helg_IDS$OldMother[i] != "0" ){
    moths[[counter]] <- Key_IDS[which(Key_IDS$OldID %in% Helg_IDS$OldMother[i]),]$NewID}
  else{ 
    moths[[counter]] <- 0
  }
}
idx <- !(sapply(moths, length))
moths[idx] <- 0
moths <- do.call(rbind, moths)
moths <- as.data.frame(moths)
Helg_IDS$NewMother <- moths$V1

#Save plink fam and Id fix files:
#update id: Old FAm, Old Individual ID, New Fam, New ID
#Uppdate Parents: Fam ID, ID, New Father, New Mother
ID_Update <- cbind.data.frame(1,Helg_IDS$OldID,1,Helg_IDS$NewID)
write.table(ID_Update, file = "Helgeland_ID_Update.txt", row.names = FALSE, col.names = FALSE)
  
Parent_Update <- cbind.data.frame(1, Helg_IDS$NewID, Helg_IDS$NewFather, Helg_IDS$NewMother)
write.table(Parent_Update, file = "Helgeland_Parent_Update.txt", row.names = FALSE, col.names = FALSE)



#Create sex fix files, if needed. Fam file may be all that is necessary actually...





#In Plink: Update IDs, Update Parents


#Make Binary Files
system("plink --file Helgeland_01_2018 --autosome-num 31 --maf 0.05 --make-bed --out workfiles")

#Update ids
system("plink --bfile workfiles --update-ids Helgeland_ID_Update.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_ids")

#update Parents
system("plink --bfile workfiles_ids --update-parents Helgeland_Parent_Update.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_parids")

#Updating Sex
plink --bfile workfiles_parids --autosome-num 31 --maf 0.05 --check-sex --make-bed --out plinksex_Helgeland

rm(pheno.file)
sex <- read.table("plinksex_Helgeland.sexcheck", header = TRUE)
pheno.file <- data.frame(1, sex$SNPSEX)
names(pheno.file) <- c("id","sex")
pheno.file <- as.data.frame(pheno.file)
pheno.file$id <- 1:length(sex$FID) 
pheno.file$sex[pheno.file$sex == 2] <- "m"
pheno.file$sex[pheno.file$sex == 1] <- "f"
pheno.file$sex[pheno.file$sex == 0] <- "u"
pheno.file$sex[pheno.file$sex == "m"] <- 1
pheno.file$sex[pheno.file$sex == "f"] <- 2
pheno.file$sex[pheno.file$sex == "u"] <- 1

write.table(pheno.file,"newfile_Helgeland.pheno",col.names = TRUE, row.names = FALSE)
pheno.file$sex <- as.numeric(pheno.file$sex)
pheno.file[3] <- pheno.file$sex
pheno.file[2] <- pheno.file$id
pheno.file[1] <- 1
write.table(pheno.file, file = "updateSex_Helgeland.txt", sep = " ", row.names = FALSE, col.names = FALSE)






plink --bfile workfiles_parids --update-sex updateSex_Helgeland.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_sexfix

plink --bfile workfiles_sexfix --mendel --autosome-num 31 --maf 0.05 --make-bed --out workfiles_mendelfix
plink --bfile workfiles_sexfix  --set-me-missing --autosome-num 31 --maf 0.05 --make-bed --out workfiles_mndmiss
plink --bfile workfiles_sexfix --recode 12 --autosome-num 31 --maf 0.05 --out test_ME-fixed

#Recode binaries to .ped
system("plink --bfile workfiles_parids --recode 12 --autosome-num 31 --maf 0.05 --out test_Helgeland")


mapfile_helg <- read.table("test_Helgeland.map")
mapfile_helg <- mapfile_helg[,-3]
write.table(mapfile_helg, file = "mapfile_helg.txt", sep=" ",row.names = FALSE,col.names = FALSE)

#5.       Make the GenABEL file
convert.snp.ped(pedfile = "test_Helgeland.ped", 
                mapfile = "mapfile_helg.txt",
                outfile = "Helgeland.gen",
                strand = "u", bcast = 1000, traits = 1, mapHasHeaderLine = FALSE)



#b. Make a phenotype file for Genabel (id, sex, knownsex)
rm(pheno.file)
sex <- read.table("plinksex_Helgeland.sexcheck", header = TRUE)
#knownsex <- read.table("All_sexings_Helgeland_200k_29112017.txt", header = TRUE)
pheno.file <- data.frame(1, sex$SNPSEX)
names(pheno.file) <- c("id","sex")
pheno.file <- as.data.frame(pheno.file)
pheno.file$id <- sex$IID
#Sexs assinged by plink sex check are incorrect because of Birds having alternate heterogametic sex
#Unknowns coded to males; ind. 136,302,340,348,417,517,537,596,627,646,663,664,706,1000,1704,1902,1940,1993,2143,2292,2380,2434,2486,2493,2538,2654,2723,2816,2978,3158,3272,3314,3388,3538,3911 mothers coded to female to match pedigree, SNP assigned male
#ind. 537,563,764,2030,3779 fathers coded to male to match pedigree
pheno.file$sex[pheno.file$sex == 2] <- "m"
pheno.file$sex[pheno.file$sex == 1] <- "f"
pheno.file$sex[pheno.file$sex == 0] <- "u"
pheno.file$sex[pheno.file$sex == "m"] <- 1
pheno.file$sex[pheno.file$sex == "f"] <- 0
pheno.file$sex[pheno.file$sex == "u"] <- 1
pheno.file$sex <- as.numeric(pheno.file$sex)
write.table(pheno.file,"Helgeland.pheno",col.names = TRUE, row.names = FALSE)

write.table(pheno.file, file = "updateSex_Helgeland.txt", sep = " ", row.names = FALSE, col.names = FALSE)

pheno.file$sex <- as.numeric(pheno.file$sex)

pheno.file[2:3] <- pheno.file[1:2]
pheno.file[1] <- 1
names(pheno.file) <-  c("FID", "IID", "sex")



sparrowgen.Helgeland <- load.gwaa.data(phe = "Helgeland.pheno",
                                       gen = "Helgeland.gen")
sparrow.abel <- sparrowgen.Helgeland
