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

coef_sim=replicate(10000,simulate_GBM_drift(sde::GBM,1,1,1,1000,10))
quantile(coef_sim,prob=c(0.025,0.975))




# Try to fit Geometric-Brownian motion which start at x0=1
# dX = rX dt + sigma X dB
end_time=1 #The end time of the Brownian motion
n_times=10000 # How many time step b/t [0,end_time]

# Set up the initial conditions of SDE
drft=expression (3*(6-x))
dfs=expression(3)
x0=1

simulate_OU_drift=function(fun_sim,r=1,sigma=1,end_T,x0,
                           n_times,k_fold,order_poly=3,drft,dfs){
  t_diff=end_T/n_times
  X=fun_sim(t0=0,T=end_T,X0=x0,N=n_times,drift=drft,sigma=dfs)
  y_data=(X[2:(n_times+1)]-X[1:n_times])/t_diff
  X_data=X[1:n_times]
  return(lm(y_data~X_data)%>%coef)
}

coef_sim=replicate(1000,simulate_OU_drift(sde::sde.sim,1,1,1,x0,
                                         10000,10,3,drft,dfs))
apply(coef_sim,1,quantile,prob=c(0.025,0.975))
plot(coef_sim[1,],coef_sim[2,])
points(18,-3,col="red")
