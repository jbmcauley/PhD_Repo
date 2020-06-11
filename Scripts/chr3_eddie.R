#Script for Eddie Large Chrs

.libPaths(c("/exports/csce/eddie/biology/groups/johnston/john/userLibrary",.libPaths()))

library(crimaptools)
library(GenABEL)
library(GenABEL.data)

sparrow.famped <- read.table("FamPed_20200414.txt", header = TRUE)
load("sparrow_abel_EDDIE.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Chr 3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)]) 


a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][7800:11800]
d <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][11700:15700]
e <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)][15600:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 3)])]




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


