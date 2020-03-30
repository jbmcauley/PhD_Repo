
#Recoding the Helgeland ID's to the New IDS established.


Key_IDS <- read.table("ID-Recode-Key.txt", header = TRUE)
Key_IDS$OldID <- as.character(Key_IDS$OldID)

Helg_IDS <- as.data.frame(sparrowgen.Helgeland@phdata$id)
names(Helg_IDS) <- "OldID"
Helg_IDS$OldID <- as.character(Helg_IDS$OldID)

Helg_IDS$NewID <- Key_IDS[which(Key_IDS$OldID %in% Helg_IDS$OldID),]$NewID

#Recode Parents

Helg_IDS$OldFather <- sparrowgen.Helgeland@phdata$Father
Helg_IDS$OldMother <- sparrowgen.Helgeland@phdata$Mother


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
write.table(ID_Update, file = "Helgeland_ID_Update.txt", row.names = FALSE)
  
Parent_Update <- cbind.data.frame(1, Helg_IDS$NewID, Helg_IDS$NewFather, Helg_IDS$NewMother)
write.table(Parent_Update, file = "Helgeland_Parent_Update.txt", row.names = FALSE)



#Create sex fix files, if needed. Fam file may be all that is necessary actually...





#In Plink: Update IDs, Update Parents


#Make Binary Files
system("plink --file Helgeland_01_2018 --autosome-num 31 --maf 0.05 --make-bed --out workfiles")

#Update ids
system("plink --bfile workfiles --update-ids Helgeland_ID_Update.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_ids")

#update Parents
system("plink --bfile workfiles_ids --update-parents Helgeland_Parent_Update.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_parids")

#Recode from binary to .ped .map
system("plink --bfile workfiles_parids --recode 12 --autosome-num 31 --maf 0.05 --out test")




mapfile_helg <- read.table("test.map")
mapfile_helg <- mapfile_helg[,-3]
write.table(mapfile_helg, file = "mapfile_helg.txt", sep=" ",row.names = FALSE,col.names = FALSE)

#5.       Make the GenABEL file
convert.snp.ped(pedfile = "test.ped", 
                mapfile = "mapfile_helg.txt",
                outfile = "Helgeland.gen",
                strand = "u", bcast = 1000, traits = 1, mapHasHeaderLine = FALSE)




