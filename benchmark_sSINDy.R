setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
import::here(sSINDy,sSINDy_lm_vectorized,order_coef,
             error_compute,.from = 'sSINDy.R')
library(microbenchmark)


size_data=1000
interval=100
#Creat a data from unifom [-100,100]
x1=runif(size_data,min=-interval,max=interval)
x2=rnorm(size_data,sd=0.1)
x3=runif(size_data,min=-interval,max=interval)
x4=runif(size_data,min=-interval,max=interval)
x5=runif(size_data,min=-interval,max=interval)
x6=runif(size_data,min=-interval,max=interval)
x7=rnorm(size_data)
X=cbind(x1,x2,x3,x4,x5,x6,x7)
# The model is a linear combination of first five variables
# x6 and x7 are data which is not used
y=10*x1+10*x2+0.1*x3+0.1*x4+.1*x5
y_noise=y+rnorm(size_data)

data_X=X
data_y=y
k_fold=5
cv_group=sample(rep(1:k_fold,size_data/k_fold),size_data)


microbenchmark(
  sSINDy(X,y,k_fold),
  sSINDy_lm_vectorized(X,y,k_fold)
)


microbenchmark(
  sSINDy(X,y,2),
  sSINDy_lm(X,y,2),
  sSINDy_lm_vectorized(X,y,2)
)



microbenchmark(
  error_compute(data_X,data_y,cv_group,k_fold),
  error_compute_MASS(data_X,data_y,cv_group,k_fold) 
)



