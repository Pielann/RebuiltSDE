import::here(sSINDy,order_coef,error_compute,
             .from = 'sSINDy.R')

library(testthat)

#Set up fake data
size_data=100
interval=100
#Creat a data from unifom [-100,100]
x1=runif(size_data,min=-interval,max=interval)
x2=rnorm(size_data,sd=0.1)
x3=runif(size_data,min=-interval,max=interval)
x4=runif(size_data,min=-interval,max=interval)
x5=runif(size_data,min=-interval,max=interval)
x6=runif(size_data,min=-interval,max=interval)
x7=rnorm(size_data)

test_that("The choice of coefficents should have the same order", {
  X=cbind(x1,x2,x3)
  y=10*x1+1*x2
  
  expect_equal(order_coef(X,y), c(3,2,1))
})


test_that("The errors of the last four should be equal or closed to zero ", {
  X=cbind(x1,x2,x3,x4,x5,x6,x7)
  y=10*x1+0.1*x2
  k_fold=5
  cv_group=sample(rep(1:k_fold,size_data/k_fold))

  expect_true(error_compute(X,y,cv_group,k_fold)<1e-10)
})


test_that("sSINDy should at least closed to the correct model", {
  X=cbind(x1,x2,x3,x4,x5,x6,x7)
  y=10*x1+0.1*x2
  k_fold=5
  mod=sSINDy(X,y,k_fold)
  
  for(i in 3:7)
    expect_lt(mod$cv_error[i],1e-10)
  
})


