####### Code for "Classification of nutritional status of children from 
####### Geometric Morphometrics: an approach to out-of-sample data" (2024)
####### Medialdea L., Arribas-Gil A., Pérez-Romero A., Gómez A. 

####### Main code and example using data for the analysis of sexual dimorphism in gorillas skulls (available from the shapes package)
####### Requires libraries shapes, roahd, MASS, class

source("aux_functions.R") 


## 1. Loading gorilla skull data - or any other raw coordinates dataset in the form of 3d arrays of dimension (k,2,n) k:nb.land, n:nb.indv.

sk_F=shapes::gorf.dat    #female skulls
sk_M=shapes::gorm.dat    #male skulls

k=dim(sk_F)[1]           #number of landmarks = 8
n1 = dim(sk_F)[3]        #number of females = 30
n2 = dim(sk_M)[3]        #number of males = 29

## 2. Full GPA
# 2.a Full GPA on each group separately, for visualization of each group

# Female skulls Procrustes coordinates
Proc_F <- shapes::procGPA(sk_F)
Proc_F <- rotate_skulls(Proc_F$rotated, Proc_F$mshape)  #Rotate so that landmarks 1 and 3 are parallel to horizontal axis

# Graphical visualization of the coordinates, the mean (black solid), functional median (red solid) and pointwise median (green dashed)
plot_shapes(Proc_F, k, n1, joinline=c(1,6:8,2:5,1), maintext="Female Gorillas skuls",med_pw=T, median=T)

# Male skulls Procrustes coordinates
Proc_M <- shapes::procGPA(sk_M)
Proc_M <- rotate_skulls(Proc_M$rotated, Proc_M$mshape)  #Rotate so that landmarks 1 and 3 are parallel to horizontal axis

# Graphical visualization of the coordinates, the mean (black solid), functional median (red solid) and pointwise median (green dashed)
plot_shapes(Proc_M, k, n2, joinline=c(1,6:8,2:5,1), maintext="Male Gorillas skuls",med_pw=T, median=T)

# 2.b Full GPA of the whole sample, for joint visualization of both groups
sk = abind::abind(sk_F,sk_M,along=3)
Proc <- shapes::procGPA(sk)
n = dim(sk)[3]
Proc <- rotate_skulls(Proc$rotated, Proc$mshape)  #Rotate so that landmarks 1 and 3 are parallel to horizontal axis

plot_shapes(Proc, k, n, joinline=c(1,6:8,2:5,1), maintext="Gorilla skulls (F+M)",med_pw=T, median=T)


## 3. Classification with leave-one-out cross validation

true_class=c(rep("F",n1),rep("M",n2))

class_results <- LOO_CV_OOS_classification(n,k,sk,true_class, all_ef=T, knnk=5)

## 4. Classification results

# 4.1 LDA with fOPA to sample mean
table(class_results$predicted_class_LDA$mean_registration,class_results$predicted_class_LDA$true_class)

# 4.2 LDA with fOPA to sample median
table(class_results$predicted_class_LDA$median_registration,class_results$predicted_class_LDA$true_class)

# 4.3 LR with fOPA to sample mean
table(class_results$predicted_class_LR$mean_registration,class_results$predicted_class_LR$true_class)

# 4.4 LR with fOPA to sample median
table(class_results$predicted_class_LR$median_registration,class_results$predicted_class_LR$true_class)

# 4.5 k-NN with fOPA to sample mean
table(class_results$predicted_class_knn$mean_registration,class_results$predicted_class_knn$true_class)

# 4.6 k-NN with fOPA to sample median
table(class_results$predicted_class_knn$median_registration,class_results$predicted_class_knn$true_class)




