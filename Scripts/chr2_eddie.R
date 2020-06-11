#Script for Eddie Large Chrs

.libPaths(c("/exports/csce/eddie/biology/groups/johnston/john/userLibrary",.libPaths()))

library(crimaptools)
library(GenABEL)
library(GenABEL.data)

sparrow.famped <- read.table("FamPed_20200414.txt", header = TRUE)
load("sparrow_abel_EDDIE.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Chr 2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][7800:11800]
d <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][11700:15700]
e <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][15600:19600]
f <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][19500:23500]
g <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)][23400:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 2)])]


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

