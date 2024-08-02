####### Code for "Classification of nutritional status of children from 
####### Geometric Morphometrics: an approach to out-of-sample data" (2024)
####### Medialdea L., Arribas-Gil A., Pérez-Romero A., Gómez A. 

####### Auxiliary functions
####### Requires libraries shapes, roahd, MASS, class


LOO_CV_classification <- function(n,k,M,true_class, all_ef=T, size=NULL, knnk=5){
  
  predicted_class_LDA <- data.frame(mean_registration=rep("",n), median_registration=rep("",n), true_class=true_class) # LDA classification
  predicted_class_LR <- data.frame(mean_registration=rep("",n), median_registration=rep("",n), true_class=true_class) # LR classification
  predicted_class_knn <- data.frame(mean_registration=rep("",n), median_registration=rep("",n), true_class=true_class) # knn classification

  g1 <- unique(true_class)[1]  # We assume 2 groups, for LR
  g2 <- unique(true_class)[2]
  
  for (i in 1:n) { #Loop over all individuals in the sample
    
    #Procrustes configuration removing individual i
    M_ =  M[,,-i]
    Proc_ = procGPA(M_)
    # mean shape
    ms = Proc_$mshape
    # median shape
    Proc_mfData = roahd::mfData(1:k,list(t(Proc_$rotated[,1,]),t(Proc_$rotated[,2,])))
    med_mfData = roahd::median_mfData(Proc_mfData) #median shape
    median_shape = mfData2matrix(med_mfData,k,2) #pasar a objeto matricial
    
    if (all_ef==T) { #We work with the residuals of the allometric effect regression
      reg_allom <- procD.lm(rotated~log(size), data=Proc_)
      X_proc_ <- reg_allom$residuals
    }else{ #We work with the Procrustes coordinates
      X_proc_ = matrix(unlist(lapply(1:(n-1),function(x){as.numeric(t(Proc_$rotated[,,x]))})),nrow=n-1,byrow=T) # Proc Coord in matrix form

    }                                                                                                           # x1,y1,x2,y2,....,xk,yk



    # Training LDA
    X_class <- data.frame(X_proc_, cl=true_class[-i])   
    mod <- MASS::lda(cl ~ ., X_class) 
    
    # Training LR (clases must be 0,1)
    X_class$cl <- (X_class$cl==g2)*1   #g1=0, g2=1
    lr <- glm(cl ~. , family=binomial(link="logit"), data=X_class)

    # fOPA of new individual (ind i) to the sample mean 
    testProc_new_mean = matrix(t(shapes::procOPA(ms,M[,,i],scale=T)$Bhat), nrow=1)  # vector format length 2k: x1,y1,x2,y2,....,xk,yk
    # fOPA of new individual (ind i) to the functional median 
    testProc_new_median = matrix(t(procOPA(median_shape,M[,,i],scale=T)$Bhat), nrow=1)  # vector format length 2k: x1,y1,x2,y2,....,xk,yk
    
    if (all_ef==T) { #We work with the residuals of the allometric effect regression
      mm = matrix(rep(apply(M[,,i],2,mean),k),ncol=2,byrow=T)
      size = sqrt((sum((M[,,i]-mm)^2)))  #Calculate centroid size for individual i
      testProc_new_mean <- testProc_new_mean - reg_allom$coefficients[1,] - reg_allom$coefficients[2,]*log(size)
      testProc_new_median <- testProc_new_median - reg_allom$coefficients[1,] - reg_allom$coefficients[2,]*log(size)
    }
    
    
    # check classification of new individual: fOPA to sample mean
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
    
    # check classification of new individual: fOPA to functional median
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

plot_shapes <- function(M, k, n, joinline, maintext,xxlim=NULL,yylim=NULL,labs=NULL,color=NULL,mean=T,med_pw=T, median=T, ld=NULL){
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
  for (i in 1:n) {
    lines(M[joinline,,i],col=color[i])
  }
  points(M[,1,],M[,2,],cex=0.25, pch=20,col="gray60")
  if(mean==T){
    ms <- apply(M, 1:2, mean) #meanshape
    lines(ms[joinline,],col="black",lwd=2)
  }
  
  if(median==T){
    # median shape
    M_mfData = roahd::mfData(1:k,list(t(M[,1,]),t(M[,2,]))) #define bi-variate functional object
    med_mfData = roahd::median_mfData(M_mfData) #median shape
    median_shape = mfData2matrix(med_mfData,k,2) #write in 3d array form
    
    lines(median_shape[joinline,],col="red3",lty=1,lwd=2)
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
  
  points(ms[joinline,],col="black",pch=pc,cex=sz)

}

### Convert mfData object to matrix (k rows, 2 columns)
mfData2matrix <- function(mfData,k,d){
  x <- matrix(0,k,2)  
  for (l in 1:d) {
    x[,l]=mfData$fDList[[l]]$values
  }
  return(x)
}
