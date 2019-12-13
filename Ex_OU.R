setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
import::here(sSINDy,order_coef,error_compute,gene_poly_data,
             .from = 'sSINDy.R')



# Try to fit Geometric-Brownian motion which start at x0=1
# dX = rX dt + sigma X dB
end_time=1 #The end time of the Brownian motion
n_times=1000 # How many time step b/t [0,end_time]

# Set up the initial conditions of SDE
drft=expression (3*(6-x))
dfs=expression(2)
x0=1

X=sde::sde.sim(t0=0,T=end_time,X0=x0,N=n_times,drift=drft,sigma=dfs)
t_diff=end_time/n_times # Length of time interval
y_data=(X[2:(n_times+1)]-X[1:n_times])/t_diff # Change of GBM b/t time
X_data=X[1:n_times] 


# Some facts for y_data
plot(X_data,y_data) #y has bigger variance as X becomes bigger
plot(X_data,type="l")
abline(18,-3)


#############################################################
## Model fitting
n_fold=10 # folds of CV

# Fit drift part
# The real model is b(X)=X
sSINDy(gene_poly_data(X_data,3),y_data,10)

lm(y_data~X_data)%>%coef

#Fit diffsion term
variation=(X[2:(n_times+1)]-X[1:n_times])*(X[2:(n_times+1)]-X[1:n_times])/t_diff
mean(variation)
plot(variation)
sSINDy(gene_poly_data(X_data,3),variation,10)

#By lasso
cv_output <- glmnet::cv.glmnet(gene_poly_data(X_data,3), 
                            y_data,intercept=FALSE)
coef(cv_output,cv_output$lambda.min)

#By relaxed lasso
cv_relaxo=relaxo::cvrelaxo(gene_poly_data(X_data,3), 
               y_data,K=5)
cv_relaxo$beta

#############################################################
# Simuation to check if we can reconstruct the drift term
simulate_OU_drift=function(fun_sim,r=1,sigma=1,end_T,x0,
                            n_times,k_fold,order_poly=3,drft,dfs){
  t_diff=end_T/n_times
  X=fun_sim(t0=0,T=end_T,X0=x0,N=n_times,drift=drft,sigma=dfs)
  y_data=(X[2:(n_times+1)]-X[1:n_times])/t_diff
  X_data=X[1:n_times]
  return(sSINDy(gene_poly_data(X_data,order_poly),
                    y_data,k_fold)$coef)
}

coef_sim=replicate(100,simulate_OU_drift(sde::sde.sim,1,1,1,x0,
                                         10000,10,3,drft,dfs))
apply(coef_sim,1,quantile,prob=c(0.025,0.975))
apply(coef_sim,1,hist)

sum(coef_sim!=0) # really few components is zero


# Simuation to check if we can reconstruct the diffusion term
simulate_OU_diffusion=function(fun_sim,r=1,sigma=1,end_T,x0,
                            n_times,k_fold,order_poly=3,drft,dfs){  
  t_diff=end_T/n_times
  X=fun_sim(t0=0,T=end_T,X0=x0,N=n_times,drift=drft,sigma=dfs)
  variation=(X[2:(n_times+1)]-X[1:n_times])*(X[2:(n_times+1)]-X[1:n_times])/t_diff
  X_data=X[1:n_times]
  return(mean(variation))
}

coef_sim=replicate(100,simulate_OU_diffusion(sde::sde.sim,1,1,1,10000,10,3,drft,dfs))
apply(coef_sim,1,summary)

#The numbers of bad result
sum(abs(coef_sim[1,])>1e-1)
sum(abs(coef_sim[2,]-1)>1e-1)

plot(abs(coef_sim[1,]))
plot(abs(coef_sim[2,]-1))
