####### Code for "Classification of nutritional status of children from 
####### Geometric Morphometrics: an approach to out-of-sample data" (2024)
####### Medialdea L., Arribas-Gil A., Pérez-Romero A., Gómez A. 

####### Simulation study - visualization of the different scenarios
####### Requires libraries shapes, roahd, MASS, class, fdasrvf
source("aux_functions.R") 



## 1. Loading female and male gorilla skull data + simulation settings

sk_F=shapes::gorf.dat
sk_M=shapes::gorm.dat

k=dim(sk_F)[1]           #number of landmarks = 8

Proc_F <- shapes::procGPA(sk_F)
Proc_M <- shapes::procGPA(sk_M)

n1=n2=50

## 2. Data generation

# 2.a Parameter settings for 5 different scenarios

###################################################################
# Scenario 1: different mean, same dispersion

ms1 <- Proc_F$mshape
ms2 <- Proc_M$mshape

s1 <- s2 <- 5*mean(abs(ms1+ms2)/2)/n1
###################################################################


###################################################################
# # Scenario 2: more distant group means, same dispersion 
# 
# ms1 <- Proc_F$mshape
# 
# lndk7 <- matrix(0,nrow=k,ncol=2)
# lndk7[7,1] <- 20
# ms2 <- Proc_M$mshape + lndk7 
# 
# s1 <- s2 <- 5*mean(abs(ms1+ms2)/2)/n1
###################################################################


###################################################################
# # Scenario 3: same mean - different dispersion (dispersion ratio between groups = 5) 
# 
# ms1 <- ms2 <- Proc_F$mshape
# 
# s1=5*mean(abs(ms1))/n1
# s2=5*s1
###################################################################


###################################################################
# # Scenario 4: same mean - different dispersion (dispersion ratio between groups = 3) 
# 
# ms1 <- ms2 <- Proc_F$mshape
# 
# s1=5*mean(abs(ms1))/n1
# s2=3*s1
###################################################################


###################################################################
# # Scenario 5: same mean - different dispersion (dispersion ratio between groups = 1.5) 
# 
# ms1 <- ms2 <- Proc_F$mshape
# 
# s1=5*mean(abs(ms1))/n1
# s2=1.5*s1
###################################################################


# 2.b Generate aligned data

G1 <- array(0,dim=c(k,2,n1))
G2 <- array(0,dim=c(k,2,n2))
for (i in 1:n1){
  
  G1[,,i] = ms1 + matrix(s1*rnorm(2*k),nrow=k,ncol=2)
  G2[,,i]= ms2 + matrix(s2*rnorm(2*k),nrow=k,ncol=2)
  
}


# 2.c Add rotation, translation and scaling to generate missaligned shapes
G = abind::abind(G1,G2,along=3)
n = dim(G)[3]  #total number of observations

for (i in 1:n){
  
  #rotation matrix
  theta <- 2*pi*runif(1)
  Rot <- matrix( c(cos(theta), -sin(theta), sin(theta), cos(theta)), ncol = 2, byrow = TRUE)
  
  #scaling factor
  alpha <- runif(1,0,10)
  
  #translation matrix
  beta <- matrix( rep(runif(2,-5,5),k),ncol=2,byrow=T)
  
  G[,,i] <- beta + alpha* G[,,i] %*% Rot   
  
}


# 2.d Full Alignment (GPA, wGPA, SRVF) of the whole sample, for visualization of both groups

# 2.d.1 GPA
Proc <- shapes::procGPA(G)
Proc <- rotate_skulls(Proc$rotated, Proc$mshape)  #Rotate so that landmarks 1 and 3 are parallel to horizontal axis

par(mfrow=c(1,2))
plot_shapes(Proc[,,1:n1], k, n1, joinline=c(1,6:8,2:5,1), maintext="Group 1 after GPA",lines=F,med_pw=F, median=T, mean=T)
plot_shapes(Proc[,,(n1+1):n], k, n2, joinline=c(1,6:8,2:5,1), maintext="Group 2 after GPA",lines=F, med_pw=F, median=T, mean=T)

par(mfrow=c(1,1))
plot_shapes(Proc, k, n, joinline=c(1,6:8,2:5,1), maintext="Both groups (G1 + G2) after GPA",med_pw=T, median=T)

# 2.d.2 wGPA
WProc <- procWGPA_silent(G)
WProc <- rotate_skulls(WProc$rotated, WProc$mshape)     #Rotate so that landmarks 1 and 3 are parallel to horizontal axis

par(mfrow=c(1,2))
plot_shapes(WProc[,,1:n1], k, n1, joinline=c(1,6:8,2:5,1), maintext="Group 1 after wGPA",lines=F,med_pw=F, median=T, mean=T)
plot_shapes(WProc[,,(n1+1):n], k, n2, joinline=c(1,6:8,2:5,1), maintext="Group 2 after wGPA",lines=F, med_pw=F, median=T, mean=T)

par(mfrow=c(1,1))
plot_shapes(WProc, k, n, joinline=c(1,6:8,2:5,1), maintext="Both groups (G1 + G2) after wGPA",med_pw=T, median=T)

# 2.d.3 SRVF
#Get global SRVF alignment
beta=array(0,dim=c(2,k,n))
for ( h in 1:n){   #transpose to get configurations in 2*k*n 3d array
  beta[,,h]=t(G[,,h])  
}
SRVF <- fdasrvf::curve_srvf_align(beta, scale=T)  #SRVF alignment
SRVFs=G
for ( h in 1:n){   #get back to 3d array with dim k*2*n
  SRVFs[,,h]=t(SRVF$betan[,,h])
}
SRVFs <- rotate_skulls(SRVFs, apply(SRVFs,1:2,mean))    #Rotate so that landmarks 1 and 3 are parallel to horizontal axis

par(mfrow=c(1,2))
plot_shapes(SRVFs[,,1:n1], k, n1, joinline=c(1,6:8,2:5,1), maintext="Group 1 after SRVF",lines=F,med_pw=F, median=T, mean=T)
plot_shapes(SRVFs[,,(n1+1):n], k, n2, joinline=c(1,6:8,2:5,1), maintext="Group 2 after SRVF",lines=F, med_pw=F, median=T, mean=T)

par(mfrow=c(1,1))
plot_shapes(SRVFs, k, n, joinline=c(1,6:8,2:5,1), maintext="Both groups (G1 + G2) after SRVF",med_pw=T, median=T)


## 3. Classification with GPA alignment, in-sample and out-of-sample LOOCV classification

true_class=c(rep("G1",n1),rep("G2",n2))

#In-sample classification
class_results_IS <- LOO_CV_IS_classification(n,k,G,true_class, all_ef=T, knnk=5)
summary_results_IS(class_results_IS)  #Summary of classification results for LDA, LR and kNN 

#Out-of-sample classification
class_results_OOS <- LOO_CV_OOS_classification(n,k,G,true_class, all_ef=T, knnk=5)
summary_results_OOS(class_results_OOS) #Summary of classification results for LDA, LR and kNN, with sample mean and funct. median as target for registration

## 4. Classification with wGPA alignment, in-sample and out-of-sample LOOCV classification

#In-sample classification
class_results_wGPA_IS <- LOO_CV_IS_classification(n,k,G,true_class, method="wGPA", all_ef=T, knnk=5)
summary_results_IS(class_results_wGPA_IS) #Summary of classification results for LDA, LR and kNN 

#Out-of-sample classification
class_results_wGPA_OOS <- LOO_CV_OOS_classification(n,k,G,true_class, method="wGPA", all_ef=T, knnk=5)
summary_results_OOS(class_results_wGPA_OOS) #Summary of classification results for LDA, LR and kNN, with sample mean and funct. median as target for registration

## 5. Classification with SRVF alignment, out-of-sample LOOCV classification
joinline=c(1,6:8,2:5,1)
jl=joinline[-length(joinline)]  #because last ldk in joinline is the same as first

#In sample
class_results_SRVF_IS <- SRVF_LOO_CV_IS_classification(n,k,G[jl,,],true_class, all_ef=T,  knnk=5) 
summary_results_IS(class_results_SRVF_IS)

#Out-of-sample
class_results_SRVF_OOS <- SRVF_LOO_CV_OOS_classification(n,k,G[jl,,],true_class, all_ef=T, knnk=5) 
summary_results_OOS(class_results_SRVF_OOS)


