####### Code for "Classification of nutritional status of children from 
####### Geometric Morphometrics: an approach to out-of-sample data" (2024)
####### Medialdea L., Arribas-Gil A., Pérez-Romero A., Gómez A. 

####### Full simulation study
####### Requires libraries shapes, roahd, MASS, class, fdasrvf
source("aux_functions.R") 

n_runs = 100  #number of simulations runs for each setting
joinline=c(1,6:8,2:5,1)
jl=joinline[-length(joinline)] 



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


tables_GPA <- vector("list",length=6)  #to store confussion matrices for LDA,LR,knn with alignment to sample mean and to functional median
for(i in 1:6){
  tables_GPA[[i]] <- matrix(0,ncol=2,nrow=2)
}
tables_wGPA <- tables_SRVF <- tables_GPA

for ( j in 1:n_runs){  #loop over simulation runs
  
  # 2.b Generate aligned data
  G1 <- array(0,dim=c(k,2,n1))
  G2 <- array(0,dim=c(k,2,n2))
  for (i in 1:n1){
  
    G1[,,i] = ms1 + matrix(s*rnorm(2*k),nrow=k,ncol=2)
    G2[,,i]= ms2 + matrix(s*rnorm(2*k),nrow=k,ncol=2)
  
  }


  # 2.c Add rotation, translation and scaling to generate missaligned shapes
  G = abind::abind(G1,G2,along=3)
  n = dim(G)[3]

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



  ## 3. Out-of-sample classification with leave-one-out cross validation: alignment methods GPA, wGPA and SRVF

  true_class=c(rep("G1",n1),rep("G2",n2))

  class_results_GPA <- LOO_CV_OOS_classification(n,k,G,true_class, all_ef=T,  knnk=5)
  class_results_wGPA <- LOO_CV_OOS_classification(n,k,G,true_class,method="wGPA", all_ef=T,  knnk=5)
  class_results_SRVF <- SRVF_LOO_CV_OOS_classification(n,k,G[jl,,],true_class, all_ef=T,  knnk=5)

  ## 4. Classification results

  tables_GPA <- append_tables_OOS(class_results_GPA,tables_GPA)
  tables_wGPA <- append_tables_OOS(class_results_wGPA,tables_wGPA)
  tables_SRVF <- append_tables_OOS(class_results_SRVF,tables_SRVF)
}

res_GPA <- acc2_list(tables_GPA, n_runs)
res_wGPA <- acc2_list(tables_wGPA, n_runs)
res_SRVF <- acc2_list(tables_SRVF, n_runs)

