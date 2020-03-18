
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



