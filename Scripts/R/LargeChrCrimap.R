library(GenABEL)
library(crimaptools)


#chromosomes 1, 2, 3, 4, 5, 6, 7, 8, 29

length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)]) # 23366
a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][7800:11800]
d <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][11700:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)])]
#e <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][15600:17477]
#f <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][19500:23366]





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr1a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1a.pre", genfile = "crimap/chr1a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr1a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1a.pre", genfile = "crimap/chr1a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr1a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr1a.gen", crimap.path =  "C:/PathApps/crimap.exe")








create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr1b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1b.pre", genfile = "crimap/chr1b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr1b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1b.pre", genfile = "crimap/chr1b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr1b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr1b.gen", crimap.path =  "C:/PathApps/crimap.exe")




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr1c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1c.pre", genfile = "crimap/chr1c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr1c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1c.pre", genfile = "crimap/chr1c.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr1c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr1c.gen", crimap.path =  "C:/PathApps/crimap.exe")






create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1d", 
                    snplist = d,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr1d.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1d.pre", genfile = "crimap/chr1d.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1d", 
                    snplist = d,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr1d.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1d.pre", genfile = "crimap/chr1d.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr1d.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr1d.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1e", 
                    snplist = e,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr1e.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1e.pre", genfile = "crimap/chr1e.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "1e", 
                    snplist = e,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr1e.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr1e.pre", genfile = "crimap/chr1e.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr1e.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr1e.gen", crimap.path =  "C:/PathApps/crimap.exe")






#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "1f", 
#                    snplist = f,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE)
#run_crimap_prepare(genfile = "crimap/chr1f.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr1f.pre", genfile = "crimap/chr1f.gen", familyPedigree = sparrow.famped)
#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "1f", 
#                    snplist = f,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE,
#                    use.mnd = TRUE)
#run_crimap_prepare(genfile = "crimap/chr1f.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr1f.pre", genfile = "crimap/chr1f.gen", familyPedigree = sparrow.famped)


#run_crimap_map(genfile = "crimap/chr1f.gen", crimap.path = "C:/PathApps/crimap.exe")
#run_crimap_chrompic(genfile = "crimap/chr1f.gen", crimap.path =  "C:/PathApps/crimap.exe")





#dir("crimap")
#sparrow.map1a <- parse_map(mapfile = "crimap/chr1a.map")

#plot(sparrow.map1a$Order,sparrow.map1a$cMPosition.Female, xlab = "", ylab = "", main = "Chr 1a",cex.main = .75)
#plot(sparrow.map1a$Order,sparrow.map1a$cMPosition.Male, xlab = "", ylab = "", main = "Chr 1a",cex.main = .75)




#chromosomes 1, 2 (32016), 3 (23580), 4 (15319), 5 (12216), 6 (6537), 7 (7122), 8 (8191), 29 (14950)


length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)]) # 23366

a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][7800:11800]
d <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][11700:15700]
e <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][15600:19600]
f <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][19500:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)])]
#g <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][23400:24242]
#h <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][27400:31400]
#i <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][31300:32016]

create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr2a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2a.pre", genfile = "crimap/chr2a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr2a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2a.pre", genfile = "crimap/chr2a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr2a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr2b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2b.pre", genfile = "crimap/chr2b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr2b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2b.pre", genfile = "crimap/chr2b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr2b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2b.gen", crimap.path =  "C:/PathApps/crimap.exe")




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr2c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2c.pre", genfile = "crimap/chr2c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr2c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2c.pre", genfile = "crimap/chr2c.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr2c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2c.gen", crimap.path =  "C:/PathApps/crimap.exe")






create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2d", 
                    snplist = d,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr2d.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2d.pre", genfile = "crimap/chr2d.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2d", 
                    snplist = d,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr2d.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2d.pre", genfile = "crimap/chr2d.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr2d.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2d.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2e", 
                    snplist = e,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr2e.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2e.pre", genfile = "crimap/chr2e.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2e", 
                    snplist = e,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr2e.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2e.pre", genfile = "crimap/chr2e.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr2e.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2e.gen", crimap.path =  "C:/PathApps/crimap.exe")






create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2f", 
                    snplist = f,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr2f.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2f.pre", genfile = "crimap/chr2f.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2f", 
                    snplist = f,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr2f.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2f.pre", genfile = "crimap/chr2f.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr2f.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2f.gen", crimap.path =  "C:/PathApps/crimap.exe")



create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2g", 
                    snplist = g,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr2g.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2g.pre", genfile = "crimap/chr2g.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "2g", 
                    snplist = g,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr2g.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr2g.pre", genfile = "crimap/chr2g.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr2g.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr2g.gen", crimap.path =  "C:/PathApps/crimap.exe")



#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "2h", 
#                    snplist = h,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE)
#run_crimap_prepare(genfile = "crimap/chr2h.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr2h.pre", genfile = "crimap/chr2h.gen", familyPedigree = sparrow.famped)
#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "2h", 
#                    snplist = h,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE,
#                    use.mnd = TRUE)
#run_crimap_prepare(genfile = "crimap/chr2h.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr2h.pre", genfile = "crimap/chr2h.gen", familyPedigree = sparrow.famped)


#run_crimap_map(genfile = "crimap/chr2h.gen", crimap.path = "C:/PathApps/crimap.exe")
#run_crimap_chrompic(genfile = "crimap/chr2h.gen", crimap.path =  "C:/PathApps/crimap.exe")



#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "2i", 
#                    snplist = i,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE)
#run_crimap_prepare(genfile = "crimap/chr2i.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr2i.pre", genfile = "crimap/chr2i.gen", familyPedigree = sparrow.famped)
#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "2i", 
#                    snplist = i,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE,
#                    use.mnd = TRUE)
#run_crimap_prepare(genfile = "crimap/chr2i.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr2i.pre", genfile = "crimap/chr2i.gen", familyPedigree = sparrow.famped)


#run_crimap_map(genfile = "crimap/chr2i.gen", crimap.path = "C:/PathApps/crimap.exe")
#run_crimap_chrompic(genfile = "crimap/chr2i.gen", crimap.path =  "C:/PathApps/crimap.exe")










#chromosomes 1, 2 (32016), 3 (23580), 4 (15319), 5 (12216), 6 (6537), 7 (7122), 8 (8191), 29 (14950)

length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)]) 


a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][7800:11800]
d <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][11700:15700]
e <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][15600:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)])]
#f <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][19500:23500]
#g <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][23400:23580]




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr3a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3a.pre", genfile = "crimap/chr3a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr3a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3a.pre", genfile = "crimap/chr3a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr3a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr3a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr3b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3b.pre", genfile = "crimap/chr3b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr3b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3b.pre", genfile = "crimap/chr3b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr3b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr3b.gen", crimap.path =  "C:/PathApps/crimap.exe")




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr3c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3c.pre", genfile = "crimap/chr3c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr3c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3c.pre", genfile = "crimap/chr3c.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr3c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr3c.gen", crimap.path =  "C:/PathApps/crimap.exe")






create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3d", 
                    snplist = d,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr3d.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3d.pre", genfile = "crimap/chr3d.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3d", 
                    snplist = d,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr3d.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3d.pre", genfile = "crimap/chr3d.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr3d.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr3d.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3e", 
                    snplist = e,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr3e.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3e.pre", genfile = "crimap/chr3e.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "3e", 
                    snplist = e,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr3e.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr3e.pre", genfile = "crimap/chr3e.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr3e.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr3e.gen", crimap.path =  "C:/PathApps/crimap.exe")






#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "3f", 
#                    snplist = f,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE)
#run_crimap_prepare(genfile = "crimap/chr3f.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr3f.pre", genfile = "crimap/chr3f.gen", familyPedigree = sparrow.famped)
#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "3f", 
#                    snplist = f,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE,
#                    use.mnd = TRUE)
#run_crimap_prepare(genfile = "crimap/chr3f.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr3f.pre", genfile = "crimap/chr3f.gen", familyPedigree = sparrow.famped)


#run_crimap_map(genfile = "crimap/chr3f.gen", crimap.path = "C:/PathApps/crimap.exe")
#run_crimap_chrompic(genfile = "crimap/chr3f.gen", crimap.path =  "C:/PathApps/crimap.exe")



#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "3g", 
#                    snplist = g,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE)
#run_crimap_prepare(genfile = "crimap/chr3g.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr3g.pre", genfile = "crimap/chr3g.gen", familyPedigree = sparrow.famped)
#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "3g", 
#                    snplist = g,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE,
#                    use.mnd = TRUE)
#run_crimap_prepare(genfile = "crimap/chr3g.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr3g.pre", genfile = "crimap/chr3g.gen", familyPedigree = sparrow.famped)


#run_crimap_map(genfile = "crimap/chr3g.gen", crimap.path = "C:/PathApps/crimap.exe")
#run_crimap_chrompic(genfile = "crimap/chr3g.gen", crimap.path =  "C:/PathApps/crimap.exe")










#chromosomes 1, 2 (32,016), 3 (23,580), 4 (15,319), 5 (12,216), 6 (6,537), 7 (7,122), 8 (8,191), 29 (14,950)

length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)]) # 23366


a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)][7800:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)])]
#d <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 4)][11700:15319]





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr4a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4a.pre", genfile = "crimap/chr4a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr4a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4a.pre", genfile = "crimap/chr4a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr4a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr4a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr4b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4b.pre", genfile = "crimap/chr4b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr4b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4b.pre", genfile = "crimap/chr4b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr4b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr4b.gen", crimap.path =  "C:/PathApps/crimap.exe")




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr4c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4c.pre", genfile = "crimap/chr4c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "4c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr4c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr4c.pre", genfile = "crimap/chr4c.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr4c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr4c.gen", crimap.path =  "C:/PathApps/crimap.exe")






#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "4d", 
#                    snplist = d,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE)
#run_crimap_prepare(genfile = "crimap/chr4d.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr4d.pre", genfile = "crimap/chr4d.gen", familyPedigree = sparrow.famped)
#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "4d", 
#                    snplist = d,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE,
#                    use.mnd = TRUE)
#run_crimap_prepare(genfile = "crimap/chr4d.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr4d.pre", genfile = "crimap/chr4d.gen", familyPedigree = sparrow.famped)


#run_crimap_map(genfile = "crimap/chr4d.gen", crimap.path = "C:/PathApps/crimap.exe")
#run_crimap_chrompic(genfile = "crimap/chr4d.gen", crimap.path =  "C:/PathApps/crimap.exe")








#chromosomes 1, 2 (32016), 3 (23580), 4 (15319), 5 (12216), 6 (6537), 7 (7122), 8 (8191), 29 (14950)

length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 5)]) # 23366


a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 5)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 5)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 5)][7800:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 5)])]
#d <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 5)][11700:12216]


create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "5a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr5a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr5a.pre", genfile = "crimap/chr5a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "5a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr5a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr5a.pre", genfile = "crimap/chr5a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr5a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr5a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "5b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr5b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr5b.pre", genfile = "crimap/chr5b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "5b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr5b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr5b.pre", genfile = "crimap/chr5b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr5b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr5b.gen", crimap.path =  "C:/PathApps/crimap.exe")




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "5c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr5c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr5c.pre", genfile = "crimap/chr5c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "5c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr5c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr5c.pre", genfile = "crimap/chr5c.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr5c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr5c.gen", crimap.path =  "C:/PathApps/crimap.exe")






#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "5d", 
#                    snplist = d,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE)
#run_crimap_prepare(genfile = "crimap/chr5d.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr5d.pre", genfile = "crimap/chr5d.gen", familyPedigree = sparrow.famped)
#create_crimap_input(gwaa.data = sparrow.abel, 
#                    familyPedigree = sparrow.famped, 
#                    analysisID = "5d", 
#                    snplist = d,
#                    outdir = "crimap", 
#                    clear.existing.analysisID = TRUE,
#                    use.mnd = TRUE)
#run_crimap_prepare(genfile = "crimap/chr5d.gen", crimap.path = "C:/PathApps/crimap.exe")
#parse_mend_err(prefile = "crimap/chr5d.pre", genfile = "crimap/chr5d.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr5d.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr5d.gen", crimap.path =  "C:/PathApps/crimap.exe")












#chromosomes 1, 2 (32016), 3 (23580), 4 (15319), 5 (12216), 6 (6537), 7 (7122), 8 (8191), 29 (14950)
length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 6)]) # 23366


a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 6)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 6)][3900:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 6)])]


create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "6a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr6a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr6a.pre", genfile = "crimap/chr6a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "6a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr6a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr6a.pre", genfile = "crimap/chr6a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr6a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr6a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "6b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr6b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr6b.pre", genfile = "crimap/chr6b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "6b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr6b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr6b.pre", genfile = "crimap/chr6b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr6b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr6b.gen", crimap.path =  "C:/PathApps/crimap.exe")











#chromosomes 1, 2 (32016), 3 (23580), 4 (15319), 5 (12216), 6 (6537), 7 (7122), 8 (8191), 29 (14950)
length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 7)]) # 23366



a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 7)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 7)][3900:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 7)])]


create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "7a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr7a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr7a.pre", genfile = "crimap/chr7a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "7a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr7a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr7a.pre", genfile = "crimap/chr7a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr7a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr7a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "7b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr7b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr7b.pre", genfile = "crimap/chr7b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "7b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr7b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr7b.pre", genfile = "crimap/chr7b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr7b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr7b.gen", crimap.path =  "C:/PathApps/crimap.exe")









#chromosomes 1, 2 (32016), 3 (23580), 4 (15319), 5 (12216), 6 (6537), 7 (7122), 8 (8191), 29 (14950)
length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 8)]) # 23366



a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 8)][1:3000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 8)][2900:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 8)])]
#c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 8)][5800:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 8)]) # 23366
#]


create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "8a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr8a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr8a.pre", genfile = "crimap/chr8a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "8a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr8a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr8a.pre", genfile = "crimap/chr8a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr8a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr8a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "8b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr8b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr8b.pre", genfile = "crimap/chr8b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "8b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr8b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr8b.pre", genfile = "crimap/chr8b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr8b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr8b.gen", crimap.path =  "C:/PathApps/crimap.exe")




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "8c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr8c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr8c.pre", genfile = "crimap/chr8c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "8c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr8c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr8c.pre", genfile = "crimap/chr8c.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr8c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr8c.gen", crimap.path =  "C:/PathApps/crimap.exe")









#chromosomes 1, 2 (32016), 3 (23580), 4 (15319), 5 (12216), 6 (6537), 7 (7122), 8 (8191), 29 (14950)
length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)]) # 23366



a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)][7800:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)])]
#d <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 29)][11700:14950]





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr29a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29a.pre", genfile = "crimap/chr29a.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29a", 
                    snplist = a,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =  TRUE)
run_crimap_prepare(genfile = "crimap/chr29a.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29a.pre", genfile = "crimap/chr29a.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr29a.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr29a.gen", crimap.path =  "C:/PathApps/crimap.exe")





create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)

run_crimap_prepare(genfile = "crimap/chr29b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29b.pre", genfile = "crimap/chr29b.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29b", 
                    snplist = b,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd =TRUE)
run_crimap_prepare(genfile = "crimap/chr29b.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29b.pre", genfile = "crimap/chr29b.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr29b.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr29b.gen", crimap.path =  "C:/PathApps/crimap.exe")




create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr29c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29c.pre", genfile = "crimap/chr29c.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29c", 
                    snplist = c,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr29c.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29c.pre", genfile = "crimap/chr29c.gen", familyPedigree = sparrow.famped)

run_crimap_map(genfile = "crimap/chr29c.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr29c.gen", crimap.path =  "C:/PathApps/crimap.exe")






create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29d", 
                    snplist = d,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE)
run_crimap_prepare(genfile = "crimap/chr29d.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29d.pre", genfile = "crimap/chr29d.gen", familyPedigree = sparrow.famped)
create_crimap_input(gwaa.data = sparrow.abel, 
                    familyPedigree = sparrow.famped, 
                    analysisID = "29d", 
                    snplist = d,
                    outdir = "crimap", 
                    clear.existing.analysisID = TRUE,
                    use.mnd = TRUE)
run_crimap_prepare(genfile = "crimap/chr29d.gen", crimap.path = "C:/PathApps/crimap.exe")
parse_mend_err(prefile = "crimap/chr29d.pre", genfile = "crimap/chr29d.gen", familyPedigree = sparrow.famped)


run_crimap_map(genfile = "crimap/chr29d.gen", crimap.path = "C:/PathApps/crimap.exe")
run_crimap_chrompic(genfile = "crimap/chr29d.gen", crimap.path =  "C:/PathApps/crimap.exe")


setwd("C:/Users/s1945757/PhD_Repo/PLINK-files 200k SNP-data/Dblxoversremoved/GenABEL_QCs/QC2")

sparrow.map1a <- parse_map(mapfile = "crimap/chr1a.map")
sparrow.map1b <- parse_map(mapfile = "crimap/chr1b.map")
sparrow.map1c <- parse_map(mapfile = "crimap/chr1c.map")
sparrow.map1d <- parse_map(mapfile = "crimap/chr1d.map")
sparrow.map1e <- parse_map(mapfile = "crimap/chr1e.map")
sparrow.map2a <- parse_map(mapfile = "crimap/chr2a.map")
sparrow.map2b <- parse_map(mapfile = "crimap/chr2b.map")
sparrow.map2c <- parse_map(mapfile = "crimap/chr2c.map")
sparrow.map2d <- parse_map(mapfile = "crimap/chr2d.map")
sparrow.map2e <- parse_map(mapfile = "crimap/chr2e.map")
sparrow.map2f <- parse_map(mapfile = "crimap/chr2f.map")
sparrow.map2g <- parse_map(mapfile = "crimap/chr2g.map")
sparrow.map3a <- parse_map(mapfile = "crimap/chr3a.map")
sparrow.map3b <- parse_map(mapfile = "crimap/chr3b.map")
sparrow.map3c <- parse_map(mapfile = "crimap/chr3c.map")
sparrow.map3d <- parse_map(mapfile = "crimap/chr3d.map")
sparrow.map3e <- parse_map(mapfile = "crimap/chr3e.map")
sparrow.map4a <- parse_map(mapfile = "crimap/chr4a.map")
sparrow.map4b <- parse_map(mapfile = "crimap/chr4b.map")
sparrow.map4c <- parse_map(mapfile = "crimap/chr4c.map")
sparrow.map5a <- parse_map(mapfile = "crimap/chr5a.map")
sparrow.map5b <- parse_map(mapfile = "crimap/chr5b.map")
sparrow.map5c <- parse_map(mapfile = "crimap/chr5c.map")
sparrow.map6a <- parse_map(mapfile = "crimap/chr6a.map")
sparrow.map6b <- parse_map(mapfile = "crimap/chr6b.map")
sparrow.map7a <- parse_map(mapfile = "crimap/chr7a.map")
sparrow.map7b <- parse_map(mapfile = "crimap/chr7b.map")
sparrow.map8a <- parse_map(mapfile = "crimap/chr8a.map")
sparrow.map8b <- parse_map(mapfile = "crimap/chr8b.map")
sparrow.map8c <- parse_map(mapfile = "crimap/chr8c.map")
sparrow.map29a <- parse_map(mapfile = "crimap/chr29a.map")
sparrow.map29b <- parse_map(mapfile = "crimap/chr29b.map")
sparrow.map29c <- parse_map(mapfile = "crimap/chr29c.map")


sparrow.cmp1a <- parse_map_chrompic(chrompicfile = "crimap/chr1a.cmp")
sparrow.cmp1b <- parse_map_chrompic(chrompicfile = "crimap/chr1b.cmp")
sparrow.cmp1c <- parse_map_chrompic(chrompicfile = "crimap/chr1c.cmp")
sparrow.cmp1d <- parse_map_chrompic(chrompicfile = "crimap/chr1d.cmp")
sparrow.cmp1e <- parse_map_chrompic(chrompicfile = "crimap/chr1e.cmp")
sparrow.cmp2a <- parse_map_chrompic(chrompicfile = "crimap/chr2a.cmp")
sparrow.cmp2b <- parse_map_chrompic(chrompicfile = "crimap/chr2b.cmp")
sparrow.cmp2c <- parse_map_chrompic(chrompicfile = "crimap/chr2c.cmp")
sparrow.cmp2d <- parse_map_chrompic(chrompicfile = "crimap/chr2d.cmp")
sparrow.cmp2e <- parse_map_chrompic(chrompicfile = "crimap/chr2e.cmp")
sparrow.cmp2f <- parse_map_chrompic(chrompicfile = "crimap/chr2f.cmp")
sparrow.cmp2g <- parse_map_chrompic(chrompicfile = "crimap/chr2g.cmp")
sparrow.cmp3a <- parse_map_chrompic(chrompicfile = "crimap/chr3a.cmp")
sparrow.cmp3b <- parse_map_chrompic(chrompicfile = "crimap/chr3b.cmp")
sparrow.cmp3c <- parse_map_chrompic(chrompicfile = "crimap/chr3c.cmp")
sparrow.cmp3d <- parse_map_chrompic(chrompicfile = "crimap/chr3d.cmp")
sparrow.cmp3e <- parse_map_chrompic(chrompicfile = "crimap/chr3e.cmp")
sparrow.cmp4a <- parse_map_chrompic(chrompicfile = "crimap/chr4a.cmp")
sparrow.cmp4b <- parse_map_chrompic(chrompicfile = "crimap/chr4b.cmp")
sparrow.cmp4c <- parse_map_chrompic(chrompicfile = "crimap/chr4c.cmp")
sparrow.cmp5a <- parse_map_chrompic(chrompicfile = "crimap/chr5a.cmp")
sparrow.cmp5b <- parse_map_chrompic(chrompicfile = "crimap/chr5b.cmp")
sparrow.cmp5c <- parse_map_chrompic(chrompicfile = "crimap/chr5c.cmp")
sparrow.cmp6a <- parse_map_chrompic(chrompicfile = "crimap/chr6a.cmp")
sparrow.cmp6b <- parse_map_chrompic(chrompicfile = "crimap/chr6b.cmp")
sparrow.cmp7a <- parse_map_chrompic(chrompicfile = "crimap/chr7a.cmp")
sparrow.cmp7b <- parse_map_chrompic(chrompicfile = "crimap/chr7b.cmp")
sparrow.cmp8a <- parse_map_chrompic(chrompicfile = "crimap/chr8a.cmp")
sparrow.cmp8b <- parse_map_chrompic(chrompicfile = "crimap/chr8b.cmp")
sparrow.cmp8c <- parse_map_chrompic(chrompicfile = "crimap/chr8c.cmp")
sparrow.cmp29a <- parse_map_chrompic(chrompicfile = "crimap/chr29a.cmp")
sparrow.cmp29b <- parse_map_chrompic(chrompicfile = "crimap/chr29b.cmp")
sparrow.cmp29c <- parse_map_chrompic(chrompicfile = "crimap/chr29c.cmp")



#Combining Chr chunks

#adding cM pos within chunks
sparrow.map1b$cMPosition.add <- sparrow.map1b$cMPosition + (sparrow.map1a[which(sparrow.map1b$SNP.Name[1] == sparrow.map1a$SNP.Name),]$cMPosition)
sparrow.map1c$cMPosition.add <- sparrow.map1c$cMPosition + (sparrow.map1b[which(sparrow.map1c$SNP.Name[1] == sparrow.map1b$SNP.Name),]$cMPosition.add)
sparrow.map1d$cMPosition.add <- sparrow.map1d$cMPosition + (sparrow.map1c[which(sparrow.map1d$SNP.Name[1] == sparrow.map1c$SNP.Name),]$cMPosition.add)
#sparrow.map1e$cMPosition.add <- sparrow.map1e$cMPosition + (sparrow.map1d[which(sparrow.map1e$SNP.Name[1] == sparrow.map1d$SNP.Name),]$cMPosition.add)
#sparrow.map1f$cMPosition.add <- sparrow.map1f$cMPosition + (sparrow.map1e[which(sparrow.map1f$SNP.Name[1] == sparrow.map1e$SNP.Name),]$cMPosition.add)

#Adding order within chunks (3899)
sparrow.map1b$Order.add <- sparrow.map1b$Order + (sparrow.map1a[which(sparrow.map1b$SNP.Name[1] == sparrow.map1a$SNP.Name),]$Order)
sparrow.map1c$Order.add <- sparrow.map1c$Order + (sparrow.map1b[which(sparrow.map1c$SNP.Name[1] == sparrow.map1b$SNP.Name),]$Order.add)
sparrow.map1d$Order.add <- sparrow.map1d$Order + (sparrow.map1c[which(sparrow.map1d$SNP.Name[1] == sparrow.map1c$SNP.Name),]$Order.add)
#sparrow.map1e$Order.add <- sparrow.map1e$Order + (sparrow.map1d[which(sparrow.map1e$SNP.Name[1] == sparrow.map1d$SNP.Name),]$Order.add)
#sparrow.map1f$Order.add <- sparrow.map1f$Order + (sparrow.map1e[which(sparrow.map1f$SNP.Name[1] == sparrow.map1e$SNP.Name),]$Order.add)


#sparrow.cmp1b$cMPosition.add <- sparrow.cmp1b$cMPosition + (sparrow.cmp1a[which(sparrow.cmp1b$SNP.Name[1] == sparrow.cmp1a$SNP.Name),]$cMPosition)
#sparrow.cmp1c$cMPosition.add <- sparrow.cmp1c$cMPosition + (sparrow.cmp1b[which(sparrow.cmp1c$SNP.Name[1] == sparrow.cmp1b$SNP.Name),]$cMPosition.add)
#sparrow.cmp1d$cMPosition.add <- sparrow.cmp1d$cMPosition + (sparrow.cmp1c[which(sparrow.cmp1d$SNP.Name[1] == sparrow.cmp1c$SNP.Name),]$cMPosition.add)
#sparrow.cmp1e$cMPosition.add <- sparrow.cmp1e$cMPosition + (sparrow.cmp1d[which(sparrow.cmp1e$SNP.Name[1] == sparrow.cmp1d$SNP.Name),]$cMPosition.add)
#sparrow.cmp1f$cMPosition.add <- sparrow.cmp1f$cMPosition + (sparrow.cmp1e[which(sparrow.cmp1f$SNP.Name[1] == sparrow.cmp1e$SNP.Name),]$cMPosition.add)

#Adding order within chunks (3899)
#sparrow.cmp1b$Order.add <- sparrow.cmp1b$Order + (sparrow.cmp1a[which(sparrow.cmp1b$SNP.Name[1] == sparrow.cmp1a$SNP.Name),]$Order)
#sparrow.cmp1c$Order.add <- sparrow.cmp1c$Order + (sparrow.cmp1b[which(sparrow.cmp1c$SNP.Name[1] == sparrow.cmp1b$SNP.Name),]$Order.add)
#sparrow.cmp1d$Order.add <- sparrow.cmp1d$Order + (sparrow.cmp1c[which(sparrow.cmp1d$SNP.Name[1] == sparrow.cmp1c$SNP.Name),]$Order.add)
#sparrow.cmp1e$Order.add <- sparrow.cmp1e$Order + (sparrow.cmp1d[which(sparrow.cmp1e$SNP.Name[1] == sparrow.cmp1d$SNP.Name),]$Order.add)
#sparrow.cmp1f$Order.add <- sparrow.cmp1f$Order + (sparrow.cmp1e[which(sparrow.cmp1f$SNP.Name[1] == sparrow.cmp1e$SNP.Name),]$Order.add)



#Chr 2b (3900): Female 463.275; Male 525.665; Combined 494.491
#Chr 2c (3901): Female 451.251    492.581   472.457
#Chr 2d (3901): 
#Chr 2e (3901): 
#Chr 2f (3901): 
#Chr 2g (3901): 
#Chr 2h (3901):
#Chr 2i (3901): 
sparrow.map2b$cMPosition.add <- sparrow.map2b$cMPosition + (sparrow.map2a[which(sparrow.map2b$SNP.Name[1] == sparrow.map2a$SNP.Name),]$cMPosition)
sparrow.map2c$cMPosition.add <- sparrow.map2c$cMPosition + (sparrow.map2b[which(sparrow.map2c$SNP.Name[1] == sparrow.map2b$SNP.Name),]$cMPosition.add)
sparrow.map2d$cMPosition.add <- sparrow.map2d$cMPosition + (sparrow.map2c[which(sparrow.map2d$SNP.Name[1] == sparrow.map2c$SNP.Name),]$cMPosition.add)
sparrow.map2e$cMPosition.add <- sparrow.map2e$cMPosition + (sparrow.map2d[which(sparrow.map2e$SNP.Name[1] == sparrow.map2d$SNP.Name),]$cMPosition.add)
sparrow.map2f$cMPosition.add <- sparrow.map2f$cMPosition + (sparrow.map2e[which(sparrow.map2f$SNP.Name[1] == sparrow.map2e$SNP.Name),]$cMPosition.add)
#sparrow.map2g$cMPosition.add <- sparrow.map2g$cMPosition + (sparrow.map2f[which(sparrow.map2g$SNP.Name[1] == sparrow.map2f$SNP.Name),]$cMPosition.add)
#sparrow.map2h$cMPosition.add <- sparrow.map2h$cMPosition + (sparrow.map2g[which(sparrow.map2h$SNP.Name[1] == sparrow.map2g$SNP.Name),]$cMPosition.add)
#sparrow.map2i$cMPosition.add <- sparrow.map2i$cMPosition + (sparrow.map2h[which(sparrow.map2i$SNP.Name[1] == sparrow.map2h$SNP.Name),]$cMPosition.add)

#Adding order within chunks (3899)
sparrow.map2b$Order.add <- sparrow.map2b$Order + (sparrow.map2a[which(sparrow.map2b$SNP.Name[1] == sparrow.map2a$SNP.Name),]$Order)
sparrow.map2c$Order.add <- sparrow.map2c$Order + (sparrow.map2b[which(sparrow.map2c$SNP.Name[1] == sparrow.map2b$SNP.Name),]$Order.add)
sparrow.map2d$Order.add <- sparrow.map2d$Order + (sparrow.map2c[which(sparrow.map2d$SNP.Name[1] == sparrow.map2c$SNP.Name),]$Order.add)
sparrow.map2e$Order.add <- sparrow.map2e$Order + (sparrow.map2d[which(sparrow.map2e$SNP.Name[1] == sparrow.map2d$SNP.Name),]$Order.add)
sparrow.map2f$Order.add <- sparrow.map2f$Order + (sparrow.map2e[which(sparrow.map2f$SNP.Name[1] == sparrow.map2e$SNP.Name),]$Order.add)
#sparrow.map2g$Order.add <- sparrow.map2g$Order + (sparrow.map2f[which(sparrow.map2g$SNP.Name[1] == sparrow.map2f$SNP.Name),]$Order.add)
#sparrow.map2h$Order.add <- sparrow.map2h$Order + (sparrow.map2g[which(sparrow.map2h$SNP.Name[1] == sparrow.map2g$SNP.Name),]$Order.add)
#sparrow.map2i$Order.add <- sparrow.map2i$Order + (sparrow.map2h[which(sparrow.map2i$SNP.Name[1] == sparrow.map2h$SNP.Name),]$Order.add)


#Chr 3a:
#Chr 3b: 
#Chr 3c: 
#Chr 3d: 
#Chr 3e: 
#Chr 3f: 
#Chr 3g:
sparrow.map3b$cMPosition.add <- sparrow.map3b$cMPosition + (sparrow.map3a[which(sparrow.map3b$SNP.Name[1] == sparrow.map3a$SNP.Name),]$cMPosition)
sparrow.map3c$cMPosition.add <- sparrow.map3c$cMPosition + (sparrow.map3b[which(sparrow.map3c$SNP.Name[1] == sparrow.map3b$SNP.Name),]$cMPosition.add)
sparrow.map3d$cMPosition.add <- sparrow.map3d$cMPosition + (sparrow.map3c[which(sparrow.map3d$SNP.Name[1] == sparrow.map3c$SNP.Name),]$cMPosition.add)
sparrow.map3e$cMPosition.add <- sparrow.map3e$cMPosition + (sparrow.map3d[which(sparrow.map3e$SNP.Name[1] == sparrow.map3d$SNP.Name),]$cMPosition.add)
#sparrow.map3f$cMPosition.add <- sparrow.map3f$cMPosition + (sparrow.map3e[which(sparrow.map3f$SNP.Name[1] == sparrow.map3e$SNP.Name),]$cMPosition.add)
#sparrow.map3g$cMPosition.add <- sparrow.map3g$cMPosition + (sparrow.map3f[which(sparrow.map3g$SNP.Name[1] == sparrow.map3f$SNP.Name),]$cMPosition.add)

#Adding order within chunks (3899)
sparrow.map3b$Order.add <- sparrow.map3b$Order + (sparrow.map3a[which(sparrow.map3b$SNP.Name[1] == sparrow.map3a$SNP.Name),]$Order)
sparrow.map3c$Order.add <- sparrow.map3c$Order + (sparrow.map3b[which(sparrow.map3c$SNP.Name[1] == sparrow.map3b$SNP.Name),]$Order.add)
sparrow.map3d$Order.add <- sparrow.map3d$Order + (sparrow.map3c[which(sparrow.map3d$SNP.Name[1] == sparrow.map3c$SNP.Name),]$Order.add)
sparrow.map3e$Order.add <- sparrow.map3e$Order + (sparrow.map3d[which(sparrow.map3e$SNP.Name[1] == sparrow.map3d$SNP.Name),]$Order.add)
#sparrow.map3f$Order.add <- sparrow.map3f$Order + (sparrow.map3e[which(sparrow.map3f$SNP.Name[1] == sparrow.map3e$SNP.Name),]$Order.add)
#sparrow.map3g$Order.add <- sparrow.map3g$Order + (sparrow.map3f[which(sparrow.map3g$SNP.Name[1] == sparrow.map3f$SNP.Name),]$Order.add)


#Chr 4a: 
#Chr 4b:
#Chr 4c:
#Chr 4d:
sparrow.map4b$cMPosition.add <- sparrow.map4b$cMPosition + (sparrow.map4a[which(sparrow.map4b$SNP.Name[1] == sparrow.map4a$SNP.Name),]$cMPosition)
sparrow.map4c$cMPosition.add <- sparrow.map4c$cMPosition + (sparrow.map4b[which(sparrow.map4c$SNP.Name[1] == sparrow.map4b$SNP.Name),]$cMPosition.add)
#sparrow.map4d$cMPosition.add <- sparrow.map4d$cMPosition + (sparrow.map4c[which(sparrow.map4d$SNP.Name[1] == sparrow.map4c$SNP.Name),]$cMPosition.add)

#Adding order within chunks (3899)
sparrow.map4b$Order.add <- sparrow.map4b$Order + (sparrow.map4a[which(sparrow.map4b$SNP.Name[1] == sparrow.map4a$SNP.Name),]$Order)
sparrow.map4c$Order.add <- sparrow.map4c$Order + (sparrow.map4b[which(sparrow.map4c$SNP.Name[1] == sparrow.map4b$SNP.Name),]$Order.add)
#sparrow.map4d$Order.add <- sparrow.map4d$Order + (sparrow.map4c[which(sparrow.map4d$SNP.Name[1] == sparrow.map4c$SNP.Name),]$Order.add)



#Chr 5

sparrow.map5b$cMPosition.add <- sparrow.map5b$cMPosition + (sparrow.map5a[which(sparrow.map5b$SNP.Name[1] == sparrow.map5a$SNP.Name),]$cMPosition)
sparrow.map5c$cMPosition.add <- sparrow.map5c$cMPosition + (sparrow.map5b[which(sparrow.map5c$SNP.Name[1] == sparrow.map5b$SNP.Name),]$cMPosition.add)
#sparrow.map5d$cMPosition.add <- sparrow.map5d$cMPosition + (sparrow.map5c[which(sparrow.map5d$SNP.Name[1] == sparrow.map5c$SNP.Name),]$cMPosition.add)
#Adding order within chunks (3899)
sparrow.map5b$Order.add <- sparrow.map5b$Order + (sparrow.map5a[which(sparrow.map5b$SNP.Name[1] == sparrow.map5a$SNP.Name),]$Order)
sparrow.map5c$Order.add <- sparrow.map5c$Order + (sparrow.map5b[which(sparrow.map5c$SNP.Name[1] == sparrow.map5b$SNP.Name),]$Order.add)
#sparrow.map5d$Order.add <- sparrow.map5d$Order + (sparrow.map5c[which(sparrow.map5d$SNP.Name[1] == sparrow.map5c$SNP.Name),]$Order.add)


#chr6
sparrow.map6b$cMPosition.add <- sparrow.map6b$cMPosition + (sparrow.map6a[which(sparrow.map6b$SNP.Name[1] == sparrow.map6a$SNP.Name),]$cMPosition)
#Adding order within chunks (3899)
sparrow.map6b$Order.add <- sparrow.map6b$Order + (sparrow.map6a[which(sparrow.map6b$SNP.Name[1] == sparrow.map6a$SNP.Name),]$Order)



#chr 7
sparrow.map7b$cMPosition.add <- sparrow.map7b$cMPosition + (sparrow.map7a[which(sparrow.map7b$SNP.Name[1] == sparrow.map7a$SNP.Name),]$cMPosition)
#Adding order within chunks (3899)
sparrow.map7b$Order.add <- sparrow.map7b$Order + (sparrow.map7a[which(sparrow.map7b$SNP.Name[1] == sparrow.map7a$SNP.Name),]$Order)


#Chr 8

sparrow.map8b$cMPosition.add <- sparrow.map8b$cMPosition + (sparrow.map8a[which(sparrow.map8b$SNP.Name[1] == sparrow.map8a$SNP.Name),]$cMPosition)
#sparrow.map8c$cMPosition.add <- sparrow.map8c$cMPosition + (sparrow.map8b[which(sparrow.map8c$SNP.Name[1] == sparrow.map8b$SNP.Name),]$cMPosition.add)
#Adding order within chunks (3899)
sparrow.map8b$Order.add <- sparrow.map8b$Order + (sparrow.map8a[which(sparrow.map8b$SNP.Name[1] == sparrow.map8a$SNP.Name),]$Order)
#sparrow.map8c$Order.add <- sparrow.map8c$Order + (sparrow.map8b[which(sparrow.map8c$SNP.Name[1] == sparrow.map8b$SNP.Name),]$Order.add)


#Chr 29

sparrow.map29b$cMPosition.add <- sparrow.map29b$cMPosition + (sparrow.map29a[which(sparrow.map29b$SNP.Name[1] == sparrow.map29a$SNP.Name),]$cMPosition)
sparrow.map29c$cMPosition.add <- sparrow.map29c$cMPosition + (sparrow.map29b[which(sparrow.map29c$SNP.Name[1] == sparrow.map29b$SNP.Name),]$cMPosition.add)
#sparrow.map29d$cMPosition.add <- sparrow.map29d$cMPosition + (sparrow.map29c[which(sparrow.map29d$SNP.Name[1] == sparrow.map29c$SNP.Name),]$cMPosition.add)
#Adding order within chunks (3899)
sparrow.map29b$Order.add <- sparrow.map29b$Order + (sparrow.map29a[which(sparrow.map29b$SNP.Name[1] == sparrow.map29a$SNP.Name),]$Order)
sparrow.map29c$Order.add <- sparrow.map29c$Order + (sparrow.map29b[which(sparrow.map29c$SNP.Name[1] == sparrow.map29b$SNP.Name),]$Order.add)
#sparrow.map29d$Order.add <- sparrow.map29d$Order + (sparrow.map29c[which(sparrow.map29d$SNP.Name[1] == sparrow.map29c$SNP.Name),]$Order.add)







#Plotting cM x SNP Order
library(ggplot2)

dev.new()
par(mfrow=c(3,7))

ggplot() +
  geom_point(data = sparrow.map1a, aes(x = Order, y = cMPosition)) +
  geom_point(data = sparrow.map1b, aes(x = Order.add, y = cMPosition.add)) +
  geom_point(data = sparrow.map1c, aes(x = Order.add, y = cMPosition.add)) +
  geom_point(data = sparrow.map1d, aes(x = Order.add, y = cMPosition.add)) +
  #geom_point(data = sparrow.map1e, aes(x = Order.add, y = cMPosition.add)) +
 # geom_point(data = sparrow.map1f, aes(x = Order.add, y = cMPosition.add)) +
  xlab('Order') +
  ylab('cM Position') +
  ggtitle('Chromosome 1: SNP Order x cM Position')

ggplot() +
  geom_point(data = sparrow.map2a, aes(x = Order, y = cMPosition)) +
  geom_point(data = sparrow.map2b, aes(x = Order.add, y = cMPosition.add)) +
  geom_point(data = sparrow.map2c, aes(x = Order.add, y = cMPosition.add)) +
  geom_point(data = sparrow.map2d, aes(x = Order.add, y = cMPosition.add)) +
  geom_point(data = sparrow.map2e, aes(x = Order.add, y = cMPosition.add)) +
  geom_point(data = sparrow.map2f, aes(x = Order.add, y = cMPosition.add)) +
  #geom_point(data = sparrow.map2g, aes(x = Order.add, y = cMPosition.add)) +
  #geom_point(data = sparrow.map2h, aes(x = Order.add, y = cMPosition.add)) +
  #geom_point(data = sparrow.map2i, aes(x = Order.add, y = cMPosition.add)) +
  xlab('Order') +
  ylab('cM Position') +
  ggtitle('Chromosome 2: SNP Order x cM Position')


ggplot() +
  geom_point(data = sparrow.map3a, aes(x = Order, y = cMPosition)) +
  geom_point(data = sparrow.map3b, aes(x = Order.add, y = cMPosition.add)) +
  geom_point(data = sparrow.map3c, aes(x = Order.add, y = cMPosition.add)) +
  geom_point(data = sparrow.map3d, aes(x = Order.add, y = cMPosition.add)) +
  geom_point(data = sparrow.map3e, aes(x = Order.add, y = cMPosition.add)) +
  #geom_point(data = sparrow.map3f, aes(x = Order.add, y = cMPosition.add)) +
  #geom_point(data = sparrow.map3g, aes(x = Order.add, y = cMPosition.add)) +
  xlab('Order') +
  ylab('cM Position') +
  ggtitle('Chromosome 3: SNP Order x cM Position')

ggplot() +
  geom_point(data = sparrow.map4a, aes(x = Order, y = cMPosition)) +
  geom_point(data = sparrow.map4b, aes(x = Order.add, y = cMPosition.add)) +
  geom_point(data = sparrow.map4c, aes(x = Order.add, y = cMPosition.add)) +
  #geom_point(data = sparrow.map4d, aes(x = Order.add, y = cMPosition.add)) +
  xlab('Order') +
  ylab('cM Position') +
  ggtitle('Chromosome 4: SNP Order x cM Position')

ggplot() +
  geom_point(data = sparrow.map5a, aes(x = Order, y = cMPosition)) +
  geom_point(data = sparrow.map5b, aes(x = Order.add, y = cMPosition.add)) +
  geom_point(data = sparrow.map5c, aes(x = Order.add, y = cMPosition.add)) +
  #geom_point(data = sparrow.map5d, aes(x = Order.add, y = cMPosition.add)) +
  xlab('Order') +
  ylab('cM Position') +
  ggtitle('Chromosome 5: SNP Order x cM Position')

ggplot() +
  geom_point(data = sparrow.map6a, aes(x = Order, y = cMPosition)) +
  geom_point(data = sparrow.map6b, aes(x = Order.add, y = cMPosition.add)) +
  xlab('Order') +
  ylab('cM Position') +
  ggtitle('Chromosome 6: SNP Order x cM Position')

ggplot() +
  geom_point(data = sparrow.map7a, aes(x = Order, y = cMPosition)) +
  geom_point(data = sparrow.map7b, aes(x = Order.add, y = cMPosition.add)) +
  xlab('Order') +
  ylab('cM Position') +
  ggtitle('Chromosome 7: SNP Order x cM Position')

ggplot() +
  geom_point(data = sparrow.map8a, aes(x = Order, y = cMPosition)) +
  geom_point(data = sparrow.map8b, aes(x = Order.add, y = cMPosition.add)) +
  #geom_point(data = sparrow.map8c, aes(x = Order.add, y = cMPosition.add)) +
  xlab('Order') +
  ylab('cM Position') +
  ggtitle('Chromosome 8: SNP Order x cM Position')

ggplot() +
  geom_point(data = sparrow.map29a, aes(x = Order, y = cMPosition)) +
  geom_point(data = sparrow.map29b, aes(x = Order.add, y = cMPosition.add)) +
  geom_point(data = sparrow.map29c, aes(x = Order.add, y = cMPosition.add)) +
 # geom_point(data = sparrow.map29d, aes(x = Order.add, y = cMPosition.add)) +
  xlab('Order') +
  ylab('cM Position') +
  ggtitle('Chromosome 29: SNP Order x cM Position')

ggarrange(a,b,c,d,e,f,g,h,i,
          labels = c("Chr1", "Chr2","Chr3","Chr4","Chr5","Chr6","Chr7",
                     "Chr8","Chr29"),
          ncol = 3, nrow = 3)






#Plotting BP vs cM

#Set up file with BP
map <- as.data.frame(sparrow.abel@gtdata@map)
names(map) <- "BP"
map$SNP.Name <- row.names(map)
row.names(map) <- NULL


#Add BP to map or cmp map files (or both)


sparrow <- merge(sparrow.map1a, map, by.x = "SNP.Name",by.y = "SNP.Name")



sparrow.map1a <- merge(sparrow.map1a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map1b <- merge(sparrow.map1b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map1c <- merge(sparrow.map1c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map1d <- merge(sparrow.map1d, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map1e <- merge(sparrow.map1e, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map1f <- merge(sparrow.map1f, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map2a <- merge(sparrow.map2a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2b <- merge(sparrow.map2b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2c <- merge(sparrow.map2c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2d <- merge(sparrow.map2d, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2e <- merge(sparrow.map2e, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2f <- merge(sparrow.map2f, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2g <- merge(sparrow.map2g, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2h <- merge(sparrow.map2h, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map2i <- merge(sparrow.map2i, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map3a <- merge(sparrow.map3a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map3b <- merge(sparrow.map3b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map3c <- merge(sparrow.map3c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map3d <- merge(sparrow.map3d, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map3e <- merge(sparrow.map3e, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map3f <- merge(sparrow.map3f, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map3g <- merge(sparrow.map3g, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map4a <- merge(sparrow.map4a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map4b <- merge(sparrow.map4b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map4c <- merge(sparrow.map4c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map4d <- merge(sparrow.map4d, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map5a <- merge(sparrow.map5a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map5b <- merge(sparrow.map5b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map5c <- merge(sparrow.map5c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map5d <- merge(sparrow.map5d, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map6a <- merge(sparrow.map6a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map6b <- merge(sparrow.map6b, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map7a <- merge(sparrow.map7a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map7b <- merge(sparrow.map7b, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map8a <- merge(sparrow.map8a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map8b <- merge(sparrow.map8b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map8c <- merge(sparrow.map8c, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map9 <- merge(sparrow.map9, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map10 <- merge(sparrow.map10, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map11 <- merge(sparrow.map11, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map12 <- merge(sparrow.map12, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map13 <- merge(sparrow.map13, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map14 <- merge(sparrow.map14, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map15 <- merge(sparrow.map15, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map17 <- merge(sparrow.map17, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map18 <- merge(sparrow.map18, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map19 <- merge(sparrow.map19, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map20 <- merge(sparrow.map20, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map21 <- merge(sparrow.map21, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map22 <- merge(sparrow.map22, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map23 <- merge(sparrow.map23, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map24 <- merge(sparrow.map24, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map25 <- merge(sparrow.map25, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map26 <- merge(sparrow.map26, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map27 <- merge(sparrow.map27, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map28 <- merge(sparrow.map28, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map29a <- merge(sparrow.map29a, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map29b <- merge(sparrow.map29b, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map29c <- merge(sparrow.map29c, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map29d <- merge(sparrow.map29d, map, by.x = "SNP.Name",by.y = "SNP.Name")

sparrow.map30 <- merge(sparrow.map30, map, by.x = "SNP.Name",by.y = "SNP.Name")
sparrow.map32 <- merge(sparrow.map32, map, by.x = "SNP.Name",by.y = "SNP.Name")



