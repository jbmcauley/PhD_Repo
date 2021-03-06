

rGLSadj <-
  function(formula.FixedEffects = y ~ 1, genabel.data, phenotype.data, id.name = "id", GRM = NULL, V = NULL, memory=1e8) {
    #Check input data
    if (class(genabel.data) != "gwaa.data") stop("The input of genabel.data is not a GenABEL object")
    if (is.null(genabel.data@phdata$id)) stop("IDs not given as id in the phdata list")
    if (!is.null(GRM)) { if(!isSymmetric(GRM)) warning("The given GRM must be a symmetric matrix!") }
    if (!is.null(V)) { if(!isSymmetric(V)) warning("The given V must be a symmetric matrix!") }
    V.input  <- V
    #require(hglm)
    #require(GenABEL)
    #Get trait name
    trait <- all.vars(formula.FixedEffects)[1]
    #Remove NAs from phenotypic data
    y.all <- phenotype.data[,names(phenotype.data)%in%trait]
    phenotype.data <- phenotype.data[!is.na(y.all),]
    #Connect IDs in GenABEL data set with IDs in the phenotype file
    id1 <- phenotype.data[,names(phenotype.data) %in% id.name] #ID for phenotype data
    id2 <- genabel.data@phdata$id #ID for genotype data
    test1 <- id1 %in% id2
    test2 <- id2 %in% id1
    genabel.data <- genabel.data[test2,] #Exclude individuals having no phenotype information
    phenotype.data <- phenotype.data[test1,] #Exclude individuals having no genotype information
    id1 <- phenotype.data[,names(phenotype.data) %in% id.name] #ID for phenotype data for cleaned data
    id2 <- genabel.data@phdata$id #ID for genotype data for cleaned data
    #####################
    #Construct incidence matrix for repeated observations
    N=length(id2)
    n=length(id1)
    indx <- numeric(n)
    for (i in 1:N) {
      indx <- indx + i * (id1 %in% id2[i])
    }
    Z.indx <- diag(N)[indx,]
    #Construct response and design matrix
    y <- phenotype.data[ , names(phenotype.data) %in% trait] #Create the response variable
    X <- model.matrix(formula.FixedEffects, data = phenotype.data) #Fixed effect design matrix
    #####################
    if (is.null(V)) {
      #Construct GRM
      if (is.null(GRM)) {
        autosomalMarkers <- which(chromosome(genabel.data)!= "X")
        GRM <- compute.GRM(genabel.data[ , snpnames(genabel.data)[autosomalMarkers]])
      }
      eig <- eigen(GRM)
      if (max(diag(GRM)) > 1.6) print("There seems to be highly inbred individuals in your data")
      if (min(eig$values < -0.5)) print("The genetic relationship matrix is far from positive definite")
      non_zero.eigenvalues <- eig$values>(1e-6) #Put numerically small eigenvalues to zero
      eig$values[ !non_zero.eigenvalues ] <- 0
      print("GRM ready")
      #####################
      #Fit hglm
      Z.GRM <- ( eig$vectors %*% diag(sqrt(eig$values)) )[indx, ]
      Z <- (cbind(Z.GRM, Z.indx))
      mod1 <- hglm(y=y, X=X, Z=Z, RandC = c(ncol(Z.GRM), ncol(Z.indx)), maxit = 200)
      if (mod1$Converge != "converged") stop("The variance component estimation did not converge in 200 iterations. Try to estimate them separately and provide the estimated (co)variance matrix V as input. \n\n")
      print("Variance component estimation ready")
      #####################
      #Construct rotation matrix
      ratio <- mod1$varRanef/mod1$varFix
      V <- constructV(Z=Z, RandC = c(ncol(Z.GRM), ncol(Z.indx)), ratio)
    }
    eig.V <- eigen(V)
    transf.matrix <- diag(1/sqrt(eig.V$values)) %*% t(eig.V$vectors)
    y.new <- transf.matrix %*% y
    X.new <- transf.matrix %*% X
    print("Rotation matrix ready")
    #####################
    #Fit a linear model for each SNP
    SNP.matrix <- as.double(genabel.data)
    if (sum(is.na(SNP.matrix)) > 0) {
      SNP.matrix <- SmoothSNPmatrix(SNP.matrix)
    }
    m <- ncol(SNP.matrix)
    p.val <- chi.sq.val <- SNP.est <- rep(1, m)
    colnames(X.new) <- as.character(1:ncol(X.new)) #To avoid columns having strange names
    print("Rotate LMM started")
    #Fit using QR factorization
    #Null model
    qr0 <- qr(X.new)
    est0 <- qr.coef(qr0, y.new)
    res <- y.new - X.new %*% est0
    n <- length(y.new)
    RSS.0 <- sum(res^2)/n
    #Split computations into reasonably sized blocks
    if (memory < n) memory <- n
    step.size <- floor(memory/n)
    steps <- ceiling(m/step.size)
    jj=1
    kk=0
    for (step.i in 1:steps) {
      if (step.i == steps) kk <- kk + m %% step.size else kk <- kk + step.size
      markers.to.fit <- jj:kk
      snp.new <- transf.matrix %*% SNP.matrix[indx, markers.to.fit]
      mm <- 0
      for (j in markers.to.fit) {
        mm <- mm+1
        X1 <- cbind(snp.new[, mm], X.new)
        qr1 <- qr(X1)
        est1 <- qr.coef(qr1, y.new)
        res <- y.new - X1 %*% est1
        RSS.1 <- sum(res^2)/n
        SNP.est[j] <- est1[1]
        LRT <- -n * (log(RSS.1) - log(RSS.0))
        p.val[j] <- 1 - pchisq(LRT, df=1)
        chi.sq.val[j] <- LRT
      }
      jj <- jj + step.size
    }
    print("Rotate LMM ready")
    #####################
    qt.results <- Create_gwaa_scan(genabel.data, p.val, SNP.est)
    qt.results@results$chi2.1df <- chi.sq.val
    if (is.null(V.input)) qt.results@call$hglm <- mod1
    return(qt.results)
  }


process_rGLSadj_results <- function(qt.results, gwaa.object){
  require(dplyr)
  if(!is.data.frame(qt.results)){
    qt.results <- results(qt.results)
    
  }
  qt.results$SNP.Name <- row.names(qt.results)
  zzz <- which(qt.results$P1df == 0)
  if(length(zzz) > 0) qt.results$P1df[zzz] <- pchisq(qt.results$chi2.1df[zzz],  1, lower.tail = F)
  qt.results <- arrange(qt.results, -chi2.1df)
  qt.results$ExpP <- seq(1/nrow(qt.results), 1, 1/nrow(qt.results))
  
  lambda <- Create_gwaa_scan(gwaa.object, qt.results$P1df, qt.results$effB)@lambda$estimate
  
  message(paste("lambda is", lambda))
  
  if(lambda > 1){
    qt.results$Pc1df <- pchisq(qt.results$chi2.1df/lambda, df = 1, lower.tail = F) 
  } else {
    qt.results$Pc1df <- qt.results$P1df
  }
  
  qt.results
}


plot_rGLSadj_results <- function(rGLSadj_results, PP = FALSE){
  
  require(dplyr)
  require(ggplot2)
  
  #rGLSadj_results <- qt.results
  rGLSadj_results$Chromosome <- as.character(rGLSadj_results$Chromosome)
  
  if("X" %in% rGLSadj_results$Chromosome){
    suppressWarnings(rGLSadj_results$Chromosome[which(rGLSadj_results$Chromosome == "X")] <- max(as.numeric(rGLSadj_results$Chromosome), na.rm = T) + 1)
    print(paste0("X recoded as ", max(as.numeric(rGLSadj_results$Chromosome))))
  }
  
  rGLSadj_results$Chromosome <- as.numeric(rGLSadj_results$Chromosome)
  rGLSadj_results <- arrange(rGLSadj_results, Chromosome, Position)
  
  rGLSadj_results$Diff <- c(diff(rGLSadj_results$Position), 0)
  rGLSadj_results$Diff[which(rGLSadj_results$Diff < 0)] <- 1000
  rGLSadj_results$Cumu <- cumsum(rGLSadj_results$Diff)
  
  chrinfo <- NULL
  
  for(j in unique(rGLSadj_results$Chromosome)){
    temp1 <- subset(rGLSadj_results, Chromosome == j)
    temp2 <- data.frame(Chromosome = j,
                        Start = temp1[1,"Cumu"],
                        Stop = temp1[nrow(temp1),"Cumu"],
                        ChrLength = temp1[nrow(temp1), "Position"])
    chrinfo <- rbind(chrinfo, temp2)
  }
  chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)
  
  if(PP){
    
    ggplot(rGLSadj_results, aes(-log10(ExpP), -log10(Pc1df))) +
      geom_point(alpha = 0.5) +
      theme_bw() +
      theme(legend.position = "none") +
      labs(x = "Expected -log10 P", y = "Observed -log10 P") +
      geom_hline(yintercept = -log10(0.05/nrow(rGLSadj_results)), linetype = "dashed") +
      geom_abline(slope = 1, intercept = 0)
    
  } else {
    
    ggplot(rGLSadj_results, aes(Cumu, -log10(Pc1df), colour = factor(Chromosome %% 2))) +
      geom_point(alpha = 0.5) +
      scale_colour_brewer(palette = "Set1") +
      theme_bw() +
      theme(legend.position = "none") +
      scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$Chromosome) +
      labs(x = "Chromosome", y = "-log10 P") +
      geom_hline(yintercept = -log10(0.05/nrow(rGLSadj_results)), linetype = "dashed")
  }
  
  
}
