setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
import::here(sSINDy,order_coef,error_compute,gene_poly_data,
             .from = 'sSINDy.R')



simulate_GBM_drift=function(fun_sim,r=1,sigma=1,end_T,
                            n_times,k_fold,order_poly=3){
  t_diff=end_T/n_times
  X=fun_sim(x=1, r,sigma, T=end_time, N=n_times)
  y_data=(X[2:(n_times+1)]-X[1:n_times])/t_diff
  X_data=X[1:n_times]
  return(lm(y_data~X_data-1)%>%coef)
}

coef_sim=sapply(c(100,1000,10000),function(nstep) 
  replicate(10000,simulate_GBM_drift(sde::GBM,1,1,1,nstep,10)))
apply(coef_sim,2,quantile,prob=c(0.025,0.975))


# Use Weighted Least Square
simulate_GBM_drift=function(fun_sim,r=1,sigma=1,end_T,
                            n_times,k_fold,order_poly=3){
  t_diff=end_T/n_times
  X=fun_sim(x=1, r,sigma, T=end_time, N=n_times)
  y_data=(X[2:(n_times+1)]-X[1:n_times])/t_diff
  X_data=X[1:n_times]
  return(lm(y_data~X_data-1,weights=(1/X_data^2))%>%coef)
}

coef_sim=sapply(c(100,1000),function(nstep) 
  replicate(1000,simulate_GBM_drift(sde::GBM,1,1,1,nstep,10)))
apply(coef_sim,2,quantile,prob=c(0.025,0.975))

