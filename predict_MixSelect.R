predict_MixSelect = function(y, X, X_test = X, y_test = y, strong = T){
  # prediction
  X_sweep_pred = sweep(X_test, 2, l, `*`)
  mm = nrow(X_sweep_pred)
  N = mm + n
  K_big = rho*fields::Exp.cov(rbind(X_sweep_pred,X_sweep), p = exp_C)
  K_big_train = rho*fields::Exp.cov(rbind(X_sweep,X_sweep), p = exp_C)
  PKPt_big = P_big%*%K_big%*%t(P_big) + sigmasq*diag(N)
  PKPt_big_train = P_big_train%*%K_big_train%*%t(P_big_train) + sigmasq*diag(2*n)
  PKPt_12 = PKPt_big[1:mm,(mm+1):N] 
  PKPt_12_train = PKPt_big_train[1:n,(n+1):(2*n)] 
  Sig_22_inv = solve(PKPt_big[(mm+1):N,(mm+1):N])
  Sig_22_inv_train = solve(PKPt_big_train[(n+1):(2*n),(n+1):(2*n)])
  
  
  mu_pred = Z_test%*%beta_z + X_test%*%beta + X_int_test%*%lambda
  mu = Z%*%beta_z + X%*%beta + X_int%*%lambda
  
  y_pred[count,] = mu_pred + PKPt_12%*%Sig_22_inv%*%(y-mu)
  #y_hat[count,] = mu + PCPt%*%CI_inv%*%(y-mu)
  y_hat[count,] = mu + PKPt_12_train%*%Sig_22_inv_train%*%(y-mu)
  mu_pred_st[count,] = mu_pred
  mu_st[count,] = mu
  
  count = count + 1
}