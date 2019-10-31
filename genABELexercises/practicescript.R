#GWAS is a quest to identify any association between a genetic polymorphism and the value of a trait of interest. Best scenario - observed association
#results from causation, that is the polymorphism studied is functionally involved in the control of the trait.
#Association has no direction, and making causal inference in epidemiology in general and in genetic epidemiology in particular is usually not possible 
#based on statistical analysis only.


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
#to 0.5 plus the genomic homozygosity, and the numbers above the diagonal tell how many SNPs were typed successfully for both subjects
#(thus UBD estimate is derived using this number of SNPs)

#2nd transform the matric to distance matrix using standard R
data1.dist <- as.dist(0.5-data1.gkin)
#finally perform Classical Multidimensinal Scaling by 
data1.mds <- cmdscale(data1.dist)
plot(data1.mds)

#identify points belonging to cluseters by:
km <- kmeans(data1.mds, centers=2, nstart = 100)
cl1 <- names(which(km$cluster==1))
cl2 <- names(which(km$cluster==2))
if (length(cl1) > length(cl2)) {x <- cl2; cl2<-cl1; cl1 <- x}
cl1
cl2
#One can form new dataset without the outliers by:
data2 <- data1[cl2,]

#With outliers dropped QC must be repeated
qc2 <- check.marker(data2, hweids=(phdata(data2)$dm2==0), fdr=0.2)
#drop the markers which failed QC by
data2 <- data2[qc2$idok, qc2$snpok]
#In this example we now have our final analysis dataset so we can attach the phenotypica data to the search path
detach(phdata(data1))
attach(phdata(data2))
#ge03d2ex.clean <- data2
#saveRDS(ge03d2ex.clean, file = "ge03d2ex.clean.RData")

#Check if complete QC improved the fit of genetic data to HWE:
descriptives.marker(data2)[2]
descriptives.marker(data2[phdata(data2)$dm2==1,])[2]
descriptives.marker(data2[phdata(data2)$dm2==0,])[2]
#No longer is there an excessive number of ind or SNPs out of HWE, can move on.

data2.qt <- qtscore(dm2, data2, trait="binomial")
lambda(data2.qt)
plot(data2.qt, df="Pc1df")
#scan summary 
descriptives.scan(data2.qt, sort="Pc1df")
#asses GW significance by
data2.qte <- qtscore(dm2, data2, times = 200, quiet=TRUE, trait="binomial")
descriptives.scan(data2.qte, sort="Pc1df")
#to adjust for sex and age we can...
data2.qtae <- qtscore(dm2~sex+age, data2, times = 200, quiet=TRUE, trait="binomial")
descriptives.scan(data2.qtae, sort="Pc1df")
#In this case there is little difference btw adjusted and unadjusted analysis, but this is not always the case; adjustment may make a study much 
#more powerful when covariates explain a large proportion of environmental trait variation.

#Finally a stratified (by BMI) analysis 
data2.qtse <- qtscore(dm2~sex+age, data2, ids=((bmi>30 &dm2==1) | dm2==0), times=200, quiet = TRUE, trait="binomial")





