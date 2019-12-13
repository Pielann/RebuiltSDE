# --- functions --- #

# To know the order for variable removing #
#' @param data_X (numeric) data of X
#' @param data_y (numeric) data of y
#' @return (numeric) How many time a variable is chosen

order_coef=function(data_X,data_y){
  n=ncol(data_X)
  order_col=rep(n,n)
  coef_dot=rep(TRUE,n)
  Xi=rep(0,n)
  
  for(i in 1:(n-1)){
    Xi[coef_dot]= MASS::ginv(data_X[,coef_dot],0) %*% data_y
    order_col[coef_dot][which.min(abs(Xi[coef_dot]))]=i
    coef_dot[coef_dot][which.min(abs(Xi[coef_dot]))]=FALSE
  }
  return(order_col)
}


# Compute error for different testing and training data #

#' @param data_X (numeric) data of X
#' @param data_y (numeric) data of y
#' @param cv_group (numeric) cv_group
#' @param k_fold (numeric) k_fold cv
#' @return (numeric) errors

error_compute=function(data_X,data_y,cv_group,k_fold){
  sum_error=0
  
  # For one dimension data_X, the generalized inverse can't be computed
  if(ncol(data_X%>%as.matrix)>1){
    for (i in 1:k_fold) {
      data_X_train=data_X[cv_group==i,]
      data_X_test=data_X[!cv_group==i,]
      data_y_train=data_y[cv_group==i]
      data_y_test=data_y[!cv_group==i]
      Xi= MASS::ginv(data_X_train,0) %*% data_y_train
      sum_error=sum_error+sum((data_y_test-data_X_test%*%Xi)^2)
    }
  }
  else{
    for (i in 1:k_fold) {
      data_X_train=data_X[cv_group==i]
      data_X_test=data_X[!cv_group==i]
      data_y_train=data_y[cv_group==i]
      data_y_test=data_y[!cv_group==i]
      if(cov(data_X_train,data_X_train)==0){
        sum_error=sum_error+sum((data_y_test-mean(data_y_train))^2)
      }
      else{
        Xi= cov(data_X_train,data_y_train)/cov(data_X_train,data_X_train)
        sum_error=sum_error+sum((data_y_test-data_X_test*Xi)^2)
      }
    }
  }
  return(sum_error)
}


# Stochastic sSINDy #

#' @param data_X (numeric) data of X
#' @param data_y (numeric) data of y
#' @param k_fold (numeric) k_fold cv
#' @return (list) list of some results

sSINDy=function(data_X,data_y,k_fold){
  library(magrittr)
  
  n_row=nrow(data_X)
  n_col=ncol(data_X)
  cv_group=sample(rep(1:k_fold,ceiling(n_row/k_fold)),n_row)
  errors=rep(0,n_col+1)
  coef_order=order_coef(data_X,data_y)
  coef_final=rep(0,n_col)
  
  #Compute errros via CV
  errors[2:(n_col+1)]=sapply((1:n_col),function(i)
  {error_compute(data_X[,coef_order>=(n_col-i+1)],data_y,cv_group,k_fold)})
  errors[1]=sum(data_y^2)*(k_fold-1)
    
  # Find the model with least errors
  min_coef=which.min(errors)
  if(min_coef==1){
    coef_final=rep(0,n_col)
  }
  else{
  coef_final[coef_order>(n_col-min_coef+1)]=
    lm(data_y~data_X[,coef_order>(n_col-min_coef+1)]-1)%>%coef()
  }
  
  return(list("coef_order"=coef_order,"cv_error"=errors,
              "min#"=min_coef-1,"coef"=coef_final))
}


#Generate the data upto nth order

#' @param X (numeric) data of X
#' @param n_order (numeric) order of expansion
#' @return (numeric) expansion matrix

gene_poly_data=function(X,n_order){
  order_X=function(i){X^i}
  return(sapply((0:n_order), order_X))
}
