####### Competitors Functions ######

Hiernet_fct = function(y, X, X_test = X, y_test = y, strong = T){
  p = ncol(X)
  n = nrow(X)
  
  fit=hierNet(X,y,lam=50, strong=strong,
              trace = 0)
  
  #extract main effects and interactions
  beta = fit$bp-fit$bn
  Omega = fit$th
  
  #prediction
  y_hat = predict(fit,X)
  err = y-y_hat
  y_pred = predict(fit,X_test)
  err_pred = y_test-y_pred
  return(list(beta = beta,Omega = Omega,
              err = err,err_pred = err_pred))
}

BKMR_fct  = function(y,X,X_test = X, y_test = y, Z = NULL){
  
  # fit
  if (is.null(Z)){
    fit = kmbayes(y = y, Z = X, iter = 1000, 
                  verbose = T, varsel = TRUE)
  }else{
    fit = kmbayes(y = y, Z = X, iter = 1000, 
                  verbose = T, varsel = TRUE,
                  X = Z)
  }
  
  
  # prediction
  y_hat = ComputePostmeanHnew(fit, Znew = X, method = "exact")
  y_hat = y_hat$postmean
  y_pred = ComputePostmeanHnew(fit, Znew = X_test, method = "exact")
  y_pred = y_pred$postmean
  
  # error 
  err = y-y_hat
  err_pred = y_test-y_pred
  
  # posterior inclusion prob 
  post_inc_prob = apply(fit$r != 0,2,mean)
  return(list(post_inc_prob = post_inc_prob,
              err = err,err_pred = err_pred,
              fit = fit))
  
  
  
}



FAMILY_fct = function(y, X, X_test = X, y_test = y,
                      alphas = c(0.01,0.5,0.99),
                      lambdas = seq(0.1,1,length = 50)){
  fit = FAMILY(X, X, y, lambdas ,alphas,
               quad = TRUE,iter=500, 
               verbose = F)
  #Find optimal model with respect to alpha and lambda
  p = ncol(X)
  n = nrow(X)
  
  yhat =  predict(fit, X,X, XequalZ = T)
  mse = apply(yhat,c(2,3), "-" ,yhat)
  mse = apply(mse^2,c(2,3),sum)
  im = which(mse==min(mse),TRUE)
  
  #extract coeffs ana main effects
  coef_FAMILY = coef(fit, XequalZ = T)[[im[2]]][[im[1]]]
  beta = numeric(p)
  beta[coef_FAMILY$mains[,1]] = coef_FAMILY$mains[,2]
  int_Family = coef_FAMILY$interacts
  Omega = matrix(0,p,p)
  for(i in 1:nrow(int_Family)){
    Omega[int_Family[i,1],int_Family[i,2]] = int_Family[i,3]
  }
  Omega = (Omega + t(Omega))/2
  
  
  #prediction
  y_hat =  predict(fit, X,X, XequalZ = T)
  err = y-y_hat
  y_pred = predict(fit, X_test,X_test, XequalZ = T)
  err_pred = y_test-y_pred
  
  return(list(beta = beta,Omega = Omega,
              err = err,err_pred = err_pred))
}



PIE_fct = function(y, X, X_test = X, y_test = y){
  
  p = ncol(X)
  n = nrow(X)
  
  #estimation
  beta = as.vector(coef(cv.glmnet(X,y,nfolds = 5),
                        s="lambda.min"))[-1];  
  Omega = PIE(X,y-X%*%beta)
  
  #prediction
  err = y-X%*%beta - as.vector(diag(X%*%Omega%*%t(X)))
  err_pred = y_test-X_test%*%beta - 
    as.vector(diag(X_test%*%Omega%*%t(X_test)))
  
  return(list(beta = beta,Omega = Omega,
              err = err,err_pred = err_pred))
}


RAMP_fct = function(y, X, X_test = X, y_test = y, hier = "Strong",
                    max.iter = 100, inter = T){
  
  #estimate model
  fit = RAMP(X,y,max.iter = max.iter,
             hier = hier, inter = inter)
  
  p = ncol(X)
  n = nrow(X)
  
  #extract coef
  #beta
  ind_main_eff = fit$mainInd
  beta = numeric(p)
  beta[ind_main_eff] = fit$beta.m
  
  #Omega
  imp_int = fit$interInd
  #int_list = sort(fit$interInd.list[[max.iter]])
  
  int_list = character(p + p*(p-1)/2)
  count = 1
  for(i in 1:p){
    for(j in i:p){
      int_list[count] = paste("X",i,"X",j,sep="")
      count = count+1
    }
  }
  
  ind = match(imp_int,int_list)
  int = numeric(p+p*(p-1)/2)
  int[ind] = fit$beta.i
  Omega = matrix(0,p,p)
  Omega[lower.tri(Omega,diag = T)] = int/2
  Omega = Omega + t(Omega)
  
  #prediction
  y_hat = predict(fit,X)
  err = y-y_hat
  y_pred = predict(fit, X_test)
  err_pred = y_test-y_pred
  
  return(list(beta = beta,Omega = Omega,
              err = err,err_pred = err_pred))
  
}