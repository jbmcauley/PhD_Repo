
#Make Binary Files
system("plink --file Pdo_200k_n3960_21032017 --autosome-num 31 --maf 0.05 --make-bed --out workfiles")

#Update ids
system("plink --bfile workfiles --update-ids updateIDs.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_ids")

#update Parents
system("plink --bfile workfiles_ids --update-parents updateParents.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_parids")

#Update sex
system("plink --bfile workfiles_parids --update-sex updateSex.txt --autosome-num 31 --maf 0.05 --make-bed --out finalbinaries")
#Recode binaries to .ped
system("plink --bfile finalbinaries --recode 12 --autosome-num 31 --maf 0.05 --out test")


#Fixing mendelian errors
system("plink --bfile finalbinaries --mendel --autosome-num 31 --maf 0.05 --make-bed --out test_mendelfix")
system("plink --bfile finalbinaries --set-me-missing --autosome-num 31 --maf 0.05 --make-bed --out test_mendelerrorssettomissing")
system("plink --bfile test_mendelerrorssettomissing --recode 12 --autosome-num 31 --maf 0.05 --out test_ME-fixed")

#Sex check file set up

"plink --file Pdo_200k_n3960_21032017 --autosome-num 31 --maf 0.05 --check-sex --make-bed --out newplinksex"










#checking paretnage; check mendelian errors per pair:

setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/")

imendel <- read.table("test_mendelfix.imendel", header = TRUE)
fmendel <- read.table("test_mendelfix.fmendel", header = TRUE)

length(which(as.numeric(imendel$N) > 10000))
length(which(as.numeric(fmendel$N) > 10000))






#Plink for Helgeland files; correcting pedigree; Updating IDs; Updating Parents; Updating Sex (May not be necessary as the .famped file in crymap will set up male/female, but can be done just in case).

#Make Binary Files
system("plink --file Helgeland_01_2018 --autosome-num 31 --maf 0.05 --make-bed --out workfiles")

#Update ids
system("plink --bfile workfiles --update-ids ID-Recode-Helgeland.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_ids")

#update Parents
system("plink --bfile workfiles_ids --update-parents updateParents_Helgeland.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_parids")


plink --bfile workfiles_parids --update-sex updateSex_Helgeland.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_sexfix

plink --bfile workfiles_sexfix --mendel --autosome-num 31 --maf 0.05 --make-bed --out workfiles_mendelfix
plink --bfile workfiles_sexfix  --set-me-missing --autosome-num 31 --maf 0.05 --make-bed --out workfiles_mndmiss
plink --bfile workfiles_sexfix --recode 12 --autosome-num 31 --maf 0.05 --out test_ME-fixed

#Recode binaries to .ped
system("plink --bfile workfiles_mendelfix --recode 12 --autosome-num 31 --maf 0.05 --out test_Helgeland")



### If recoding sex necessary see below
###
###

#Update Sex

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

system("plink --bfile workfiles_parids --update-sex updateSex_Helgeland.txt --autosome-num 31 --maf 0.05 --make-bed --out finalbinaries")







#Setting up LEPMAP Plink file. By individual chromsome

#example chromsomes

setwd("C:/Users/s1945757/Downloads/Sparrow/input")
sparrowgen8 <- sparrowgen[,sparrowgen@gtdata@chromosome == 8] 

export.plink(sparrowgen2, filebasename = "chr2", phenotypes = NULL, transpose = FALSE, export012na = TRUE)
export.plink(sparrowgen3, filebasename = "chr3", phenotypes = NULL, transpose = FALSE, export012na = TRUE)
export.plink(sparrowgen4, filebasename = "chr4", phenotypes = NULL, transpose = FALSE, export012na = TRUE)
export.plink(sparrowgen5, filebasename = "chr5", phenotypes = NULL, transpose = FALSE, export012na = TRUE)
export.plink(sparrowgen6, filebasename = "chr6", phenotypes = NULL, transpose = FALSE, export012na = TRUE)
export.plink(sparrowgen7, filebasename = "chr7", phenotypes = NULL, transpose = FALSE, export012na = TRUE)
export.plink(sparrowgen8, filebasename = "chr8", phenotypes = NULL, transpose = FALSE, export012na = TRUE)



