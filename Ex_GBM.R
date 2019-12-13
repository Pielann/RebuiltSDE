setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
import::here(sSINDy,order_coef,error_compute,gene_poly_data,
             .from = 'sSINDy.R')



# Try to fit Geometric-Brownian motion which start at x0=1
# dX = rX dt + sigma X dB
end_time=1 #The end time of the Brownian motion
n_times=10000 # How many time step b/t [0,end_time]

X=sde::GBM(x=1, r=1,sigma=.1, T=end_time, N=n_times) # Generate BM
t_diff=end_time/n_times # Length of time interval
y_data=(X[2:(n_times+1)]-X[1:n_times])/t_diff # Change of GBM b/t time
X_data=X[1:n_times] 

# Some facts for y_data
plot(X_data,y_data) #y has bigger variance as X becomes bigger
plot(X_data,type="l")


#############################################################
## Model fitting
n_fold=10 # folds of CV

# Fit drift part
# The real model is b(X)=X
sSINDy(gene_poly_data(X_data,3),y_data,10)

lm(y_data~X_data-1)%>%coef

#Fit diffsion term
variation=(X[2:(n_times+1)]-X[1:n_times])*(X[2:(n_times+1)]-X[1:n_times])/t_diff
mean(variation)
plot(variation)
sSINDy_new(gene_poly_data(X_data,3),variation,10)

#By lasso
cv_output <- glmnet::cv.glmnet(gene_poly_data(X_data,3), 
                            y_data,intercept=FALSE)
coef(cv_output,cv_output$lambda.min)

#By relaxed lasso
cv_relaxo=relaxo::cvrelaxo(gene_poly_data(X_data,3), 
               y_data,K=5)
cv_relaxo$beta
relaxo::plot.relaxo(cv_relaxo)

#############################################################
# Simuation to check if we can reconstruct the drift term
simulate_GBM_drift=function(fun_sim,r=1,sigma=1,end_T,
                            n_times,k_fold,order_poly=3){
  t_diff=end_T/n_times
  X=fun_sim(x=1, r,sigma, T=end_time, N=n_times)
  y_data=(X[2:(n_times+1)]-X[1:n_times])/t_diff
  X_data=X[1:n_times]
  return(sSINDy(gene_poly_data(X_data,order_poly),
                    y_data,k_fold)$coef)
}

coef_sim=replicate(100,simulate_GBM_drift(sde::GBM,1,1,1,1000,10))
apply(coef_sim,1,summary)
sum(coef_sim!=0) # really few components is zero


# Simuation to check if we can reconstruct the diffusion term
simulate_GBM_diffusion=function(fun_sim,r=1,sigma=1,end_T,
                            n_times,k_fold,order_poly=3){  
  t_diff=end_T/n_times
  X=fun_sim(x=1, r,sigma, T=end_time, N=n_times)
  variation=(X[2:(n_times+1)]-X[1:n_times])*(X[2:(n_times+1)]-X[1:n_times])/t_diff
  X_data=X[1:n_times]
  return(sSINDy(gene_poly_data(X_data,order_poly),
                    variation,k_fold)$coef)
}

coef_sim=replicate(1000,simulate_BM_diffusion(sde::BM,1,1000,10))
apply(coef_sim,1,summary)

#The numbers of bad result
sum(abs(coef_sim[1,])>1e-1)
sum(abs(coef_sim[2,]-1)>1e-1)

plot(abs(coef_sim[1,]))
plot(abs(coef_sim[2,]-1))
