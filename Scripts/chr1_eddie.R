#Script for Eddie Large Chrs

.libPaths(c("/exports/csce/eddie/biology/groups/johnston/john/userLibrary",.libPaths()))

library(crimaptools)
library(GenABEL)
library(GenABEL.data)

sparrow.famped <- read.table("FamPed_20200414.txt", header = TRUE)
load("sparrow_abel_EDDIE.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Chr 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

a <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][1:4000]
b <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][3900:7900]
c <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][7800:11800]
d <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][11700:15700]
e <- names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)][15600:length(names(sparrow.abel@gtdata@chromosome)[which(sparrow.abel@gtdata@chromosome == 1)])]






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




