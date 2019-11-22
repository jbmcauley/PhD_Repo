library(GenABEL)

dir()

system("plink --file Pdo_200k_n3960_21032017 --horse --make-bed --out Pdo_200k_n3960_21032017")


famfile <- read.table("Pdo_200k_n3960_21032017.fam", stringsAsFactors = F)
famfile <- famfile[,2:5]
names(famfile) <- c("id", "Father", "Mother", "sex")


write.table(famfile, "Pdo_200k_n3960_21032017.phe", row.names = F, sep = "\t", quote = F)

#~~ Create map file

mapfile <- read.table("Pdo_200k_n3960_21032017.map")
head(mapfile)
write.table(mapfile[,c(1, 2, 4)], "Pdo_200k_n3960_21032017.genabelmap", row.names = F, col.names = F, quote = F, sep = "\t")

#~~ Make GenAbel files

convert.snp.ped(pedfile = "Pdo_200k_n3960_21032017.ped", 
                mapfile = "Pdo_200k_n3960_21032017.genabelmap",
                outfile = "Pdo_200k_n3960_21032017.genabel",
                strand = "u", bcast = 10000, traits = 1, mapHasHeaderLine = F)

###### PROBLEM HERE - FIX WITH LETTERS????

sparrowgen <- load.gwaa.data(phenofile = "Pdo_200k_n3960_21032017.phe",
                             genofile  = "Pdo_200k_n3960_21032017.genabel")



save(sparrowgen, file = "sparrowgen.RData")


