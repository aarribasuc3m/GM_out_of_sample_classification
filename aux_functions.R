####### Code for "Classification of nutritional status of children from 
####### Geometric Morphometrics: an approach to out-of-sample data" (2024)
####### Medialdea L., Arribas-Gil A., Pérez-Romero A., Gómez A. 

####### Auxiliary functions
####### Requires libraries shapes, roahd, MASS, class, fdasrvf

### Out-of-sample LOOCV classification with GPA or wGPA alignment
LOO_CV_OOS_classification <- function(n,k,M,true_class, method="GPA", all_ef=T, knnk=5){
  
  predicted_class_LDA <- data.frame(mean_registration=rep("",n), median_registration=rep("",n), true_class=true_class) # LDA classification
  predicted_class_LR <- data.frame(mean_registration=rep("",n), median_registration=rep("",n), true_class=true_class) # LR classification
  predicted_class_knn <- data.frame(mean_registration=rep("",n), median_registration=rep("",n), true_class=true_class) # knn classification

  g1 <- unique(true_class)[1]  # We assume 2 groups, for LR
  g2 <- unique(true_class)[2]
  
  size=rep(0,n)
  for(i in 1:n){
    mm = matrix(rep(apply(M[,,i],2,mean),k),ncol=2,byrow=T)
    size[i] = sqrt((sum((M[,,i]-mm)^2)))  #Calculate centroid size for individual h
  }
  
  for (i in 1:n) { #Loop over all individuals in the sample
    
    #(w)Procrustes configuration removing individual i
    M_ =  M[,,-i]
    if(method=="wGPA"){
      Proc_ = procWGPA_silent(M_)
      # If Sigmak is singular, we find a non-singular small perturbation of it (required later for pairwise wOPA alignment, rows 69 and 75)
      Sigmak=Proc_$Sigmak
      eig=eigen(Sigmak, symmetric = TRUE)
      while (min(eig$values)==0){
        Sigmak = Sigmak+10^-5*diag(rep(rnorm(1),k))
        eig=eigen(Sigmak, symmetric = TRUE)
      }
      
    }else{
      Proc_ = shapes::procGPA(M_)
    }
      
    # mean shape
    ms = Proc_$mshape
    # median shape
    Proc_mfData = roahd::mfData(1:k,list(t(Proc_$rotated[,1,]),t(Proc_$rotated[,2,])))
    med_mfData = roahd::median_mfData(Proc_mfData) #median shape
    median_shape = mfData2matrix(med_mfData,k,2) #pasar a objeto matricial
    
    if (all_ef==T) { #We work with the residuals of the allometric effect regression
      reg_allom <- geomorph::procD.lm(rotated~log(size), data=Proc_)
      X_proc_ <- reg_allom$residuals
    }else{ #We work with the Procrustes coordinates
      X_proc_ = matrix(unlist(lapply(1:(n-1),function(x){as.numeric(t(Proc_$rotated[,,x]))})),nrow=n-1,byrow=T) # Proc Coord in matrix form
    }
                                                                                                              # x1,y1,x2,y2,....,xk,yk

    # Training LDA
    X_class <- data.frame(X_proc_, cl=true_class[-i])   
    mod <- MASS::lda(cl ~ ., X_class) 
    
    # Training LR (clases must be 0,1)
    X_class$cl <- (X_class$cl==g2)*1   #g1=0, g2=1
    lr <- glm(cl ~. , family=binomial(link="logit"), data=X_class)

    # (w)fOPA of new individual (ind i) to the sample mean 
    if(method=="wGPA"){
      testProc_new_mean = matrix(t(WOPA(ms,M[,,i],Sigma=Sigmak)), nrow=1)  # vector format length 2k: x1,y1,x2,y2,....,xk,yk
    }else{
      testProc_new_mean = matrix(t(shapes::procOPA(ms,M[,,i])$Bhat), nrow=1)  # vector format length 2k: x1,y1,x2,y2,....,xk,yk
    }
    # (w)fOPA of new individual (ind i) to the functional median 
    if(method=="wGPA"){
      testProc_new_median = matrix(t(WOPA(median_shape,M[,,i],Sigma=Sigmak)), nrow=1)  # vector format length 2k: x1,y1,x2,y2,....,xk,yk
    }else{
      testProc_new_median = matrix(t(shapes::procOPA(median_shape,M[,,i])$Bhat), nrow=1)  # vector format length 2k: x1,y1,x2,y2,....,xk,yk
    }
    
    #If Allometric Effect
    if (all_ef==T) { #We work with the residuals of the allometric effect regression
      testProc_new_mean <- testProc_new_mean - reg_allom$coefficients[1,] - reg_allom$coefficients[2,]*log(size[i])
      testProc_new_median <- testProc_new_median - reg_allom$coefficients[1,] - reg_allom$coefficients[2,]*log(size[i])
    }
    
    
    # check classification of new individual: (w)fOPA to sample mean
    new_lr <- as.data.frame(testProc_new_mean)
    colnames(new_lr) <- names(lr$coefficients)[-1]
    
    # LDA
    P.LDAmean <- predict(mod, new_lr)
    
    #LR
    lr_mean <- predict(lr, newdata =  new_lr, type="response") 
    prob_lr_mean <- ifelse(lr_mean < 0.5, g1, g2)
    
    #k-NN 
    train = X_proc_
    test = testProc_new_mean
    cl = factor(true_class[-i])
    pred_class_meanknn <- class::knn(train, test, cl, k=knnk)
    
    # check classification of new individual: (w)fOPA to functional median
    new_lr <- as.data.frame(testProc_new_median)
    colnames(new_lr) <- names(lr$coefficients)[-1]
    
    # LDA
    P.LDAmedian <- predict(mod, new_lr)
    
    #LR
    lr_median <- predict(lr, newdata =  new_lr, type="response") 
    prob_lr_median <- ifelse(lr_median < 0.5, g1, g2)
    
    #k-NN 
    test = testProc_new_median
    pred_class_medianknn <- class::knn(train, test, cl, k=knnk)
    
    predicted_class_LDA[i,1:2] <- c(P.LDAmean$class, P.LDAmedian$class)
    predicted_class_LR[i,1:2] <- c(prob_lr_mean, prob_lr_median)
    predicted_class_knn[i,1:2] <- c(pred_class_meanknn, pred_class_medianknn)
    # 
  }
  
  return(list(predicted_class_LDA = predicted_class_LDA, predicted_class_LR=predicted_class_LR, predicted_class_knn = predicted_class_knn))
  
  }
  
### Silent wGPA function (remove prints and plots)

procWGPA_silent <- function (x, fixcovmatrix = FALSE, initial = "Identity", maxiterations = 10, 
                          scale = TRUE, reflect = FALSE, prior = "Exponential", diagonal = TRUE, 
                          sampleweights = "Equal")  #just removing print and plots from shapes::procWGPA
{
  X <- x
  priorargument <- prior
  alpha <- "not estimated"
  gamma <- "not estimated"
  k <- dim(X)[1]
  n <- dim(X)[3]
  m <- dim(X)[2]
  if (initial[1] == "Identity") {
    Sigmak <- diag(k)
  }
  else {
    if (initial[1] == "Rawdata") {
      tol <- 1e-10
      if (m == 2) {
        Sigmak <- diag(diag(var(t(X[, 1, ]))) + diag(var(t(X[, 
                                                             2, ]))))/2 + tol
      }
      if (m == 3) {
        Sigmak <- diag(diag(var(t(X[, 1, ]))) + diag(var(t(X[, 
                                                             2, ]))) + diag(var(t(X[, 3, ]))))/3 + tol
      }
    }
    else {
      Sigmak <- initial
    }
  }
  mu <- shapes::procGPA(X, scale = scale)$mshape
  if (fixcovmatrix[1] != FALSE) {
    Sigmak <- fixcovmatrix
  }
  ans <- shapes:::procWGPA1(X, mu, metric = Sigmak, scale = scale, reflect = reflect, 
                   sampleweights = sampleweights)
  if ((maxiterations > 1) && (fixcovmatrix[1] == FALSE)) {
    ans0 <- ans
    dif <- 999999
    it <- 1
    while ((dif > 1e-05) && (it < maxiterations)) {
      it <- it + 1
      # if (it == 2) {
      #   cat("Differences in norm of Sigma estimates... \n ")
      # }
      if (prior[1] == "Identity") {
        prior <- diag(k)
      }
      if (prior[1] == "Inversegamma") {
        lam <- eigen(ans$Sigmak)$values
        nlam <- min(c(n * m - m - 3, k - 3))
        mu <- mean(1/lam[1:(nlam)])
        alpha <- 1/mu
        out <- nlm(iglogl, p = c(1), lam = lam, nlam = nlam)
        gamma <- abs(out$estimate[1])
        alpha <- gamma/mean(1/lam[1:nlam])
        newmetric <- n * m/(n * m + 2 * (1 + gamma)) * 
          (ans$Sigmak + (2 * alpha/(n * m)) * diag(k))
      }
      if (prior[1] == "Exponential") {
        lam <- eigen(ans$Sigmak)$values
        nlam <- min(c(n * m - m - 2, k - 2))
        mu <- mean(1/lam[1:(nlam)])
        alpha <- 1/mu
        gamma <- 1
        newmetric <- n * m/(n * m + 2 * (2)) * (ans$Sigmak + 
                                                  (2 * alpha/(n * m)) * diag(k))
      }
      if (is.double(prior[1])) {
        newmetric <- (ans$Sigmak + prior)
      }
      if (diagonal == TRUE) {
        newmetric <- diag(diag(newmetric))
      }
      if (fixcovmatrix[1] != FALSE) {
        newmetric <- fixcovmatrix
      }
      ans2 <- shapes:::procWGPA1(X, ans$mshape, metric = newmetric, 
                        scale = scale, sampleweights = sampleweights)
      # plotshapes(ans2$rotated)
      dif <- shapes:::Enorm((ans$Sigmak - ans2$Sigmak))
      ans <- ans2
      # cat(c(it, " ", dif, " \n"))
    }
  }
  if ((priorargument[1] == "Exponential") || (priorargument[1] == 
                                              "Inversegamma")) {
    ans$alpha <- alpha
    ans$gamma <- gamma
  }
  # cat(" \n")
  ans
}

### Pairwise wPA alignment 
WOPA<-function(ms,x,Sigma){ #ms is the target to align x to, Sigma is a given covariance matrix
    
    k=nrow(x) #number of landmarks
    M1_ms=array(dim=c(k,ncol(ms),2))
    M1_ms[,,1]=ms
    M1_ms[,,2]=x
    
    wOPA <- procWGPA_silent(M1_ms,fixcovmatrix = Sigma)
    x_reg <- wOPA$rotated[,,2]       #aligned x coordinates
    ms_reg <- wOPA$rotated[,,1]      #aligned ms coordinates
    
    OPAms <- shapes::procOPA(ms, ms_reg)  #ms_reg should remain equal to ms: estimate rotation and scaling from ms
    ans <- OPAms$s*(x_reg - rep(1,k) %*% t(apply(ms_reg,2,mean))) %*% OPAms$R + rep(1,k) %*%t(apply(ms,2,mean)) #rotate, scale and translate x_reg to align it to ms
    
    return(ans)
}


### Out-of-sample LOOCV classification with SRVF alignment
SRVF_LOO_CV_OOS_classification <- function(n,k,M,true_class, all_ef=T,  knnk=5){
  
  predicted_class_LDA <- data.frame(mean_registration=rep("",n), median_registration=rep("",n), true_class=true_class) # LDA classification
  predicted_class_LR <- data.frame(mean_registration=rep("",n), median_registration=rep("",n), true_class=true_class) # LR classification
  predicted_class_knn <- data.frame(mean_registration=rep("",n), median_registration=rep("",n), true_class=true_class) # knn classification
  
  g1 <- unique(true_class)[1]  # We assume 2 groups, for LR
  g2 <- unique(true_class)[2]
  
  size=rep(0,n)
  for(i in 1:n){
    mm = matrix(rep(apply(M[,,i],2,mean),k),ncol=2,byrow=T)
    size[i] = sqrt((sum((M[,,i]-mm)^2)))  #Calculate centroid size for individual h
  }
  
  for (i in 1:n) { #Loop over all individuals in the sample
    
    #Raw coordinates configuration removing individual i
    M_ =  M[,,-i]
    
    #Get global SRVF alignment
    beta=array(0,dim=c(2,k,n-1))
    for ( h in 1:(n-1)){   #transpose to get configurations in 2*k*(n-1) 3d array
      beta[,,h]=t(M_[,,h])  
    }
    SRVF <- fdasrvf::curve_srvf_align(beta)  #SRVF alignment
    SRVFs=M_
    for ( h in 1:(n-1)){   #get back to 3d array with dim k*2*(n-1)
      SRVFs[,,h]=t(SRVF$betan[,,h])
    }
    
    #SRVF aligned coordinates into (n-1)*2k matrix
    X_proc_ = matrix(unlist(lapply(1:(n-1),function(x){as.numeric(t(SRVFs[,,x]))})),nrow=n-1,byrow=T) # x1,y1,x2,y2,....,xk,yk

    #Mean shape and median shape of aligned coordinates
    ms = apply(SRVFs[,,],1:2, mean)
    # ms = t(fdasrvf::curve_karcher_mean(beta)$mu) #an alternative could be the Karcher mean of SRVFs
    Proc_mfData = roahd::mfData(1:k,list(t(SRVFs[,1,]),t(SRVFs[,2,])))
    med_mfData = roahd::median_mfData(Proc_mfData) #median shape
    median_shape = mfData2matrix(med_mfData,k,2) #pasar a objeto matricial
    # median_shape = t(curve_karcher_mean(beta,ms="median")$mu) #an alternative could be the Karcher median of SRVFs
  
    
    # If Allometric Regression: We obtain the new set of 2*k variables
    if (all_ef==T) { #We work with the residuals of the allometric effect regression
        reg_allom = X_proc_
        coefs = matrix(0,nrow=2,ncol=2*k)
        for(j in 1:(2*k)){
            lm_regall <- lm(X_proc_[,j]~ log(size[-i]))
            reg_allom[,j] = lm_regall$res
            coefs[,j] = lm_regall$coefficients
        }
        X_proc_ <- reg_allom
      }                                                                                                           # x1,y1,x2,y2,....,xk,yk
    
    
    # Training LDA
    X_class <- data.frame(X_proc_, cl=true_class[-i])   
    mod <- MASS::lda(cl ~ ., X_class) 
    
    # Training LR (clases must be 0,1)
    X_class$cl <- (X_class$cl==g2)*1   #g1=0, g2=1
    lr <- glm(cl ~. , family=binomial(link="logit"), data=X_class)
    
    # Pairwise SRVF alignment of new individual (ind i) to the sample mean of SRVFs
    testProc_new_mean =  matrix( curve_pair_align2(t(ms), t(M[,,i]))$beta2n, nrow=1) # vector format length 2k: x1,y1,x2,y2,....,xk,yk
    # SRVF of new individual (ind i) to the functional median of SRVFs
    testProc_new_median = matrix( curve_pair_align2(t(median_shape), t(M[,,i]))$beta2n, nrow=1) # vector format length 2k: x1,y1,x2,y2,....,xk,yk
    
    
    # If Allometric Regression
    if (all_ef==T) { #We work with the residuals of the allometric effect regression
      testProc_new_mean <- testProc_new_mean - coefs[1,] - coefs[2,]*log(size[i])
      testProc_new_median <- testProc_new_median - coefs[1,] - coefs[2,]*log(size[i])
    }
    
    
    # check classification of new individual: alignment to sample mean
    new_lr <- as.data.frame(testProc_new_mean)
    colnames(new_lr) <- names(lr$coefficients)[-1]
    
    # LDA
    P.LDAmean <- predict(mod, new_lr)
    
    #LR
    lr_mean <- predict(lr, newdata =  new_lr, type="response") 
    prob_lr_mean <- ifelse(lr_mean < 0.5, g1, g2)
    
    #k-NN 
    train = X_proc_
    test = testProc_new_mean
    cl = factor(true_class[-i])
    pred_class_meanknn <- class::knn(train, test, cl, k=knnk)
    
    # check classification of new individual: alignment to functional median
    new_lr <- as.data.frame(testProc_new_median)
    colnames(new_lr) <- names(lr$coefficients)[-1]
    
    # LDA
    P.LDAmedian <- predict(mod, new_lr)
    
    #LR
    lr_median <- predict(lr, newdata =  new_lr, type="response") 
    prob_lr_median <- ifelse(lr_median < 0.5, g1, g2)
    
    #k-NN 
    test = testProc_new_median
    pred_class_medianknn <- class::knn(train, test, cl, k=knnk)
    
    predicted_class_LDA[i,1:2] <- c(P.LDAmean$class, P.LDAmedian$class)
    predicted_class_LR[i,1:2] <- c(prob_lr_mean, prob_lr_median)
    predicted_class_knn[i,1:2] <- c(pred_class_meanknn, pred_class_medianknn)
    
    
    # 
  }
  
  return(list(predicted_class_LDA = predicted_class_LDA, predicted_class_LR=predicted_class_LR, predicted_class_knn = predicted_class_knn))
  
}

### Pairwise SRVF alignment 
curve_pair_align2 <- function (beta1, beta2) #just a revised version of fdasrvf::curve_pair_align 
{
  T1 = ncol(beta1)
  # centroid1 = fdasrvf:::calculatecentroid(beta1)  #we do not want to translate curve 1, which is the mean/median pattern
  # dim(centroid1) = c(length(centroid1), 1)
  # beta1 = beta1 - fdasrvf:::repmat(centroid1, 1, T1)
  centroid2 = fdasrvf:::calculatecentroid(beta2)
  dim(centroid2) = c(length(centroid2), 1)
  beta2 = beta2 - fdasrvf:::repmat(centroid2, 1, T1)
  q1 = fdasrvf::curve_to_q(beta1)$q
  out = fdasrvf:::reparam_curve(beta1, beta2)
  beta2n = out$R %*% fdasrvf:::shift_f(beta2, out$tau)
  gamI = fdasrvf:::invertGamma(out$gam)
  beta2n = fdasrvf:::group_action_by_gamma_coord(beta2n, gamI)
  q2n = fdasrvf:::curve_to_q(beta2n)$q
  # return(list(beta2n = out$beta2new, q2n = q2n, gam = gamI, q1 = q1))   #this is the original return line in the package, there's a typo, it should be beta2n = beta2n 
  
  return(list(beta2n = beta2n, q2n = q2n, gam = gamI, q1 = q1))
}


### In-sample LOOCV classification with (w)GPA alignment
LOO_CV_IS_classification <- function(n,k,M,true_class, method="GPA", all_ef=T, knnk=5){ 
  
  predicted_class_LDA <- data.frame(pred_class=rep("",n), true_class=true_class) # LDA classification
  predicted_class_LR <- data.frame(pred_class=rep("",n), true_class=true_class) # LR classification
  predicted_class_knn <- data.frame(pred_class=rep("",n), true_class=true_class) # knn classification
  
  g1 <- unique(true_class)[1]  # We assume 2 groups, for LR
  g2 <- unique(true_class)[2]
  
  #Get global (w)GPA
  if(method=="wGPA"){
    Proc = procWGPA_silent(M)
  }else{
    Proc = shapes::procGPA(M)
  }
  
  if (all_ef==T) { #We work with the residuals of the allometric effect regression
    reg_allom <- geomorph::procD.lm(rotated~log(size), data=Proc)
    X_proc <- reg_allom$residuals
  }else{ #We work with the Procrustes coordinates
    X_proc = matrix(unlist(lapply(1:n,function(x){as.numeric(t(Proc$rotated[,,x]))})),nrow=n,byrow=T) # Proc Coord in matrix form
  }
                                                                                                             # x1,y1,x2,y2,....,xk,yk
  #LOO classification
  X_class <- data.frame(X_proc, cl=true_class) 
  # Training LR (clases must be 0,1)
  X_classLR <- X_class
  X_classLR$cl <- (X_classLR$cl==g2)*1   #g1=0, g2=1
  
  for (i in 1:n){  
    # Training LDA
    mod <- MASS::lda(cl ~ ., X_class[-i,]) 
    
    lr <- glm(cl ~. , family=binomial(link="logit"), data=X_classLR[-i,])
    
    # LDA
    P.LDA <- predict(mod, X_class[i,1:(2*k)])
    
    #LR
    lr_p <- predict(lr, newdata =  X_classLR[i,1:(2*k)], type="response") 
    prob_lr <- ifelse(lr_p < 0.5, g1, g2)
    
    #k-NN 
    train = X_proc[-i,]
    test = X_proc[i,]
    cl = factor(true_class[-i])
    pred_class_knn <- class::knn(train, test, cl, k=knnk)
    
    predicted_class_LDA[i,1] <- P.LDA$class
    predicted_class_LR[i,1] <- prob_lr
    predicted_class_knn[i,1] <- pred_class_knn
    # 
  }
  
  return(list(predicted_class_LDA = predicted_class_LDA, predicted_class_LR=predicted_class_LR, predicted_class_knn = predicted_class_knn))
  
}


### In-sample LOOCV classification with SRVF alignment
SRVF_LOO_CV_IS_classification <- function(n,k,M,true_class, all_ef=T, knnk=5){
  
  predicted_class_LDA <- data.frame(pred_class=rep("",n), true_class=true_class) # LDA classification
  predicted_class_LR <- data.frame(pred_class=rep("",n), true_class=true_class) # LR classification
  predicted_class_knn <- data.frame(pred_class=rep("",n), true_class=true_class) # knn classification
  
  g1 <- unique(true_class)[1]  # We assume 2 groups, for LR
  g2 <- unique(true_class)[2]
  
  size=rep(0,n)
  for(i in 1:n){
    mm = matrix(rep(apply(M[,,i],2,mean),k),ncol=2,byrow=T)
    size[i] = sqrt((sum((M[,,i]-mm)^2)))  #Calculate centroid size for individual h
  }

  #Get global SRVF alignment
  beta=array(0,dim=c(2,k,n))
  for ( h in 1:n){   #transpose to get configurations in 2*k*n 3d array
    beta[,,h]=t(M[,,h])  
  }
  SRVF <- fdasrvf::curve_srvf_align(beta)  #SRVF alignment
  SRVFs=M
  for ( h in 1:n){   #get back to 3d array with dim k*2*n
    SRVFs[,,h]=t(SRVF$betan[,,h])
  }
    
  #SRVF aligned coordinates into n*2k matrix
  X_proc = matrix(unlist(lapply(1:(n-1),function(x){as.numeric(t(SRVFs[,,x]))})),nrow=n,byrow=T) # x1,y1,x2,y2,....,xk,yk
    
  # If Allometric Regression: We obtain the new set of 2*k variables
  if (all_ef==T) { #We work with the residuals of the allometric effect regression
      reg_allom = X_proc
      coefs = matrix(0,nrow=2,ncol=2*k)
      for(j in 1:(2*k)){
        lm_regall <- lm(X_proc[,j]~ log(size))
        reg_allom[,j] = lm_regall$res
        coefs[,j] = lm_regall$coefficients
      }
      X_proc <- reg_allom
  }                                                                                                           # x1,y1,x2,y2,....,xk,yk
    
  #LOO classification
  X_class <- data.frame(X_proc, cl=true_class) 
  # Training LR (clases must be 0,1)
  X_classLR <- X_class
  X_classLR$cl <- (X_classLR$cl==g2)*1   #g1=0, g2=1
  
  for (i in 1:n){  
    # Training LDA
    mod <- MASS::lda(cl ~ ., X_class[-i,]) 
    
    lr <- glm(cl ~. , family=binomial(link="logit"), data=X_classLR[-i,])
    
    # LDA
    P.LDA <- predict(mod, X_class[i,1:(2*k)])
    
    #LR
    lr_p <- predict(lr, newdata =  X_classLR[i,1:(2*k)], type="response") 
    prob_lr <- ifelse(lr_p < 0.5, g1, g2)
    
    #k-NN 
    train = X_proc[-i,]
    test = X_proc[i,]
    cl = factor(true_class[-i])
    pred_class_knn <- class::knn(train, test, cl, k=knnk)
    
    predicted_class_LDA[i,1] <- P.LDA$class
    predicted_class_LR[i,1] <- prob_lr
    predicted_class_knn[i,1] <- pred_class_knn
    # 
  }
  
  return(list(predicted_class_LDA = predicted_class_LDA, predicted_class_LR=predicted_class_LR, predicted_class_knn = predicted_class_knn))
  
}


### Plot function
plot_shapes <- function(M, k, n, joinline, maintext,lines=T,xxlim=NULL,yylim=NULL,labs=NULL,color=NULL,mean=T,med_pw=T, median=T, ld=NULL){
  if(n==1){
    MM <- array(dim=c(k,2,1))
    MM[,,1] <- M
    M <- MM
  }
  if(is.null(xxlim)){xxlim=range(M[,1,])}
  if(is.null(yylim)){yylim=range(M[,2,])}
  if(is.null(labs)){labs=1:k}
  if(is.null(color)){color=rep("gray90",n)}
  
  plot(M[,1,], M[,2,], axes=F, asp=1, cex=0.25, xlab="", ylab="", pch=20, 
       ylim=yylim, xlim=xxlim, main= maintext ) 
  
  if(lines==T){
     for (i in 1:n) {
        lines(M[joinline,,i],col=color[i])
      }
  }
  
  points(M[,1,],M[,2,],cex=0.25, pch=20,col="gray60")
  if(mean==T){
    ms <- apply(M, 1:2, mean) #meanshape
    lines(ms[joinline,],col="black",lwd=ifelse(lines,2,1))
  }
  
  if(median==T){
    # median shape
    M_mfData = roahd::mfData(1:k,list(t(M[,1,]),t(M[,2,]))) #define bi-variate functional object
    med_mfData = roahd::median_mfData(M_mfData) #median shape
    median_shape = mfData2matrix(med_mfData,k,2) #write in 3d array form
    
    lines(median_shape[joinline,],col="red3",lty=1,lwd=ifelse(lines,2,1))
  }
  
  if(med_pw==T){
    # median pointwise
    md <- apply(M, 1:2, median) #mediann pointwise
    c1<-c()
    c2<-c()
    for ( j in 1:(k-1)){
      c1<- c(c1,seq(md[joinline[j],1],md[joinline[j+1],1],,5))
      c2<- c(c2,seq(md[joinline[j],2],md[joinline[j+1],2],,5))
    }
    md_seq <- cbind(c1,c2)
    lines(md[joinline,],lty=2,col="green")
  }
  
  if(mean==T){
  #Landmarks representation (ld represents landmarks, the remaining points in joinline are semilandmarks)
  if (!is.null(ld)){
    id=which(joinline%in% ld)
  }else{
    id=1:length(joinline)
  }
  pc=rep(1,length(joinline))
  pc[id] = 16
  sz=rep(1.5,length(joinline))
  sz[id]=2
  
  if(lines==T){
    points(ms[joinline,],col="black",pch=pc,cex=sz)
  }
  
  }

}


### Convert mfData object to matrix (k rows, 2 columns)
mfData2matrix <- function(mfData,k,d){
  x <- matrix(0,k,2)  
  for (l in 1:d) {
    x[,l]=mfData$fDList[[l]]$values
  }
  return(x)
}


### Classification performance metrics from 2*2 confusion matrices
acc2<- function(tb){
  #tb is a confussion matrix, 2 rows, 2 columns
  
  acc <- round(sum(diag(tb))/sum(tb),4)
  
  sens <- round(tb[1,1]/sum(tb[,1]),4)
  spec <- round(tb[2,2]/sum(tb[,2]),4)
  
  return(list(acc=acc, sens=sens,spec=spec))
}

### Summary of in-sample classification results, 2*2 confusion matrices
summary_results_IS<- function(class_results){
  
  # 1 LDA 
  tr <- table(class_results$predicted_class_LDA$pred_class,class_results$predicted_class_LDA$true_class)
  cat("LDA")
  print(tr)
  print(acc2(tr))
  # 2 LR 
  tr <- table(class_results$predicted_class_LR$pred_class,class_results$predicted_class_LR$true_class)
  cat("LR")
  print(tr)
  print(acc2(tr))
  # 3 k-NN
  tr <- table(class_results$predicted_class_knn$pred_class,class_results$predicted_class_knn$true_class)
  cat("kNN")
  print(tr)
  print(acc2(tr))
}


### Summary of out-of-sample classification results, 2*2 confusion matrices
summary_results_OOS<- function(class_results){

# 1 LDA with alignment to sample mean
tr <- table(class_results$predicted_class_LDA$mean_registration,class_results$predicted_class_LDA$true_class)
cat("LDA mean alignment")
print(tr)
print(acc2(tr))
# 2 LDA with alignment to sample median
tr <- table(class_results$predicted_class_LDA$median_registration,class_results$predicted_class_LDA$true_class)
cat("LDA median alignment")
print(tr)
print(acc2(tr))
# 3 LR with alignment to sample mean
tr <- table(class_results$predicted_class_LR$mean_registration,class_results$predicted_class_LR$true_class)
cat("LR mean alignment")
print(tr)
print(acc2(tr))
# 4 LR with alignment to sample median
tr <- table(class_results$predicted_class_LR$median_registration,class_results$predicted_class_LR$true_class)
cat("LR median alignment")
print(tr)
print(acc2(tr))
# 5 k-NN with alignment to sample mean
tr <- table(class_results$predicted_class_knn$mean_registration,class_results$predicted_class_knn$true_class)
cat("k-NN mean alignment")
print(tr)
print(acc2(tr))
# 6 k-NN with alignment to sample median
tr <- table(class_results$predicted_class_knn$median_registration,class_results$predicted_class_knn$true_class)
cat("k-NN median alignment")
print(tr)
print(acc2(tr))
}

table2 <- function(x,y){ #to get a 2*2 table even when all individuals are classified in the same group
  
  ta <- table(x,y)
  if (nrow(ta)<2){
    R <- matrix(0,ncol=2,nrow=2)
    colnames(R) <- c("G1","G2")
    rownames(R) <- c("G1","G2")
    r<- match(rownames(ta),rownames(R))
    R[r,]<-ta
    ta <- R
  }
  
  return(ta)
} 

### Append out-of-sample classification results, 2*2 confusion matrices
append_tables_OOS<- function(class_results,tables){
  
  # 4.1 LDA with fOPA to sample mean
  tables[[1]] <- tables[[1]] + table2(class_results$predicted_class_LDA$mean_registration,class_results$predicted_class_LDA$true_class)
  
  # 4.2 LDA with fOPA to sample median
  tables[[2]] <- tables[[2]] + table2(class_results$predicted_class_LDA$median_registration,class_results$predicted_class_LDA$true_class)
  
  # 4.3 LR with fOPA to sample mean
  tables[[3]] <- tables[[3]] + table2(class_results$predicted_class_LR$mean_registration,class_results$predicted_class_LR$true_class)
  
  # 4.4 LR with fOPA to sample median
  tables[[4]] <- tables[[4]] + table2(class_results$predicted_class_LR$median_registration,class_results$predicted_class_LR$true_class)
  
  # 4.5 k-NN with fOPA to sample mean
  tables[[5]] <- tables[[5]] + table2(class_results$predicted_class_knn$mean_registration,class_results$predicted_class_knn$true_class)
  
  # 4.6 k-NN with fOPA to sample median
  tables[[6]] <- tables[[6]] + table2(class_results$predicted_class_knn$median_registration,class_results$predicted_class_knn$true_class)
  
  return(tables)
}

### Classification performance metrics from list of 2*2 confusion matrices
acc2_list <- function(tables, n_runs){
  
  res=vector("list",length=length(tables))
  
  for (i in 1:length(tables)){
    res[[i]] <- acc2(tables[[i]]/n_runs)
  }
  
  return(res)
  
}


### Rotate gorilla skulls so that line between landmarks 1 and 3 is parallel to horizontal axis
rotate_skulls <- function(Proc_rot,ms){
n = dim(Proc_rot)[3]
alpha=ms[3,]-ms[1,]
theta = atan(alpha[2]/alpha[1])
theta = ifelse(alpha[1]<0, pi+theta,theta)

Rot = Proc_rot
for (i in 1:n){
  Rot[,,i] = Rot[,,i]  %*% matrix( c(cos(theta), -sin(theta), sin(theta), cos(theta)), ncol = 2, byrow = TRUE)
}

return(Rot)
}
