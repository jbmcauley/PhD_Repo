library(GenABEL)
data("ge03d2ex")

an0 <- qtscore(dm2, ge03d2ex, trait = "binomial")


#evidence for inflation; The 'se' produced by estlambda can not be used to test if inflation is sig. and make conclustions about the presence of sig or insig stratification
lambda(an0)


#check.marker is used to conduct 1st round of QC. If there is possibility of heterogenity of the study data it is a good idea to skip HW checks in this 
#first round, this can be done with the p.level argument set to 0
qc1 <- check.marker(ge03d2ex, p.level = 0)
#This simple QC function will normally NEVER be run as you need to provide a number of thresholds specific to the platform you used for genotyping.


#Call Rate
#The call rate for a given SNP is defined as the proportion of individualsin the study for which the corresponding SNP information is not missing. 
#If we filter using a call rate of 95%, meaning we retain SNPs for which there is less than 5% missingdata. More stringent cut points (e.g., less than 5%) may be employed in smaller sample settings.


#To generate new data set with only individuals who have passed your QC
data1 <- ge03d2ex[qc1$idok, qc1$snpok]
#Any residual sporadic X-errors (male heterozygosity), these can and should be fixed (set to NA) by:
data1 <- Xfix(data1)

detach(phdata(ge03d2ex))
attach(phdata(data1))
descriptives.marker(data1)[2]
descriptives.marker(data1[dm2==1,])[2]
descriptives.marker(data1[dm2==0,])[2]
#Fit to HWE has improved after QC however cases still have excess number of markers. This may be due to genetic sub-structure.

#Finding Genetic Sub-structure:
#2nd round of QC - detection of genetic outliers which may contaminate results.
#1st compute a matrix of genomic kinship btw all pairs of ind, using only autosomal markers by:
data1.gkin <- ibs(data1[,autosomal(data1)], weight="freq")
#Numbers below the diagonal show the genomic estimate of kinship (aka 'genomic kinship' or 'genome-wide IBD'), numbers above the diagonal correspond
#to 0.5 plus the genomic homozygosity, and the numbers above the diagonal 