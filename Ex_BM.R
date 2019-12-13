setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
import::here(sSINDy,order_coef,error_compute,gene_poly_data,
             .from = 'sSINDy.R')


# Try to fit Brownian motion which start at x0=10
end_time=1 #The end time of the Brownian motion
n_times=100000 # How many time step b/t [0,end_time]

X=sde::BM(x=10, t0=0, T=end_time, N=n_times) # Generate BM
t_diff=end_time/n_times # Length of time interval
y_data=(X[2:(n_times+1)]-X[1:n_times]) # Change of BM b/t time
X_data=X[1:n_times] 

#The relation b/t X and dX
plot(X_data,y_data)

# Some facts for y_data
plot(density(y_data)) #This plot should be similar to normal N(0,var=dt)
mean(y_data)
var(y_data)

#############################################################
## Model fitting
n_fold=10 # folds of CV

# Fit drift part
# The real model is b(X)=0
sSINDy(gene_poly_data(X_data,3),y_data,10)

#Compute the diffsion term
sd(y_data/sqrt(t_diff))



#############################################################
# Drift Term
# Simuation to check if we can reconstruct the drift term
simulate_BM_drift=function(fun_sim,end_T,n_times,k_fold,order_poly=3){
  t_diff=end_T/n_times
  X=fun_sim(x=10, t0=0, end_T, n_times)
  y_data=(X[2:(n_times+1)]-X[1:n_times])/t_diff
  X_data=X[1:n_times]
  return(sSINDy(gene_poly_data(X_data,order_poly),
                    y_data,k_fold)$coef)
}

# Replicate some times and compute the CI for Drift Term
coef_sim=replicate(1000,simulate_BM_drift(sde::BM,1,1000,10))
apply(coef_sim,1,quantile,prob=c(0.025,0.975))
sum(coef_sim!=0) # really few components is zero

# Diffusion Term
# Simuation to check if we can reconstruct the diffusion term
simulate_BM_diffusion=function(fun_sim,end_T,n_times){
  t_diff=end_T/n_times
  X=fun_sim(x=10, t0=0, end_T, n_times)
  y_data=(X[2:(n_times+1)]-X[1:n_times])
  return("Diffusion"=sd(y_data/sqrt(t_diff)))
}

# Replicate some times and compute the CI for Diffusion Term
coef_sim=replicate(100,simulate_BM_diffusion(sde::BM,1,10000))
quantile(coef_sim,prob=c(0.025,0.975))



