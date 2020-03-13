
#Make Binary Files
system("plink --file Pdo_200k_n3960_21032017 --autosome-num 31 --maf 0.05 --make-bed --out workfiles")

#Update ids
system("plink --bfile workfiles --update-ids updateIDs.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_ids")


system("plink --bfile workfiles_ids --update-parents updateParents.txt --autosome-num 31 --maf 0.05 --make-bed --out workfiles_parids")


system("plink --bfile workfiles_parids --update-sex updateSex.txt --autosome-num 31 --maf 0.05 --make-bed --out finalbinaries")

system("plink --bfile finalbinaries --recode 12 --autosome-num 31 --maf 0.05 --out test")

system("plink --bfile finalbinaries --mendel --autosome-num 31 --maf 0.05 --out test_mendelfix")


#Sex checks


"plink --file Pdo_200k_n3960_21032017 --autosome-num 31 --maf 0.05 --check-sex --make-bed --out newplinksex"
