# Gibbs sampler for MixSelect regression, based on a decomposition of the regression
#           surface in:
#           - main effects
#           - interactions with heredity constraint
#           - nonlinear deviation distributed as projected GP
# ArXiv: https://arxiv.org/abs/1911.01910
# Email: ff31@duke.edu for questions

# ARGUMENTS: y: response vector;
#            X: predictor matrix (n x p) to be included in projected GP;
#            Z: covariate adjustments matrix (n x q);
#            nrun: number of iterations;
#            burn: burn-in period;
#            thin: thinning interval; lambda_prior = "cont_spike",
#            tau: variance hyperparameter for continuous spike and slab on main effects and interactions
#            c: shrinkage parameter for continuous spike and slab on main effects and interactions
#            exp_C: exponent of exponential covariance, default is  (2.2) in paper
#            lambda_prior: "cont_spike" or "discrete_spike", the  latter is extremely slow and  not  suggested
#            na: T if analysis with missing data
#            na: T if analysis with missing data
#            k_na: number of latent factors for analysis with NAs, refer to Section 5.3 in paper
#            list_na: vector of dimension p+q indicating which variables in X and Z need to be included in  Factor Model,
#                   this is because maybe we don't want to use a factor model for all the covariates, or maybe there are some variables
#                   that are not missing
#            sd_lambda: hyperparamenter for element of \Lambda matrix
#            lod: T if analysis with data below the Limit of Detection (LOD)
#            X_lod: Matrix of dimension nxp, 1 if observation i,j has been observed  below LOD and 1  otherwise
#            lod_vec: vector or dimension  px1 containing  the LOD for each  variable  in  X

MixSelect = function(y, X , Z = NULL,
                     nrun = 2000, burn = 1000,thin = 1,
                     tau = 1, c = 0.1, exp_C = 2,
                     heredity = "strong", lambda_prior = "cont_spike",
                     na = F, k_na = NULL, list_na = NULL,sd_lambda = 0.1,
                     lod = F, X_lod = NULL, lod_vec = NULL,
                     verbose = T){
  
  
  n = nrow(X); p = ncol(X); n_test =  nrow(X_test)
  
  # --- Confounders ---
  if(is.null(Z)){
    Z = matrix(0,nrow = nrow(X), ncol = 1);
  } 

  # --- NA --- #
  if(na){
    
    XZ = cbind(X,Z); XZ_test = cbind(X_test,Z_test)
    
    if(is.null(list_na)) list_na = rep(T,ncol(XZ))
    W = XZ[,list_na]; W_test = XZ_test[,list_na]
    p_na = ncol(W)
    if(is.null(k_na)) k_na = floor(2*log(p_na))
    W_na = W %>% is.na()
    W_test_na = W_test %>% is.na()
    
    # set na to zero just to compute some of the quantities later
    X_na = X %>% is.na(); X_test_na = X_test %>% is.na()
    X[X_na] = 0; X_test[X_test_na] = 0;
    Z_na = Z %>% is.na(); Z_test_na = Z_test %>% is.na()
    Z[Z_na] = 0; Z_test[Z_test_na] = 0;
    
    # parameters of factor model
    Lambda_na = matrix(0, nrow = p_na, ncol = k_na)
    Plam_na =  matrix(sd_lambda, nrow = p_na, ncol = k_na)   
    eta_na = matrix(0, nrow = n, ncol = k_na)
    eta_test_na = matrix(0,n_test,k_na)
    ps_na = rep(1,p_na)
    
  }
  
  # --- LOD --- #
  if(lod){
    
    X[X_lod] = 0
    X_lod[X_na] = F
    
  }
  
  # --- interactions --- #
  X_int = model.matrix(y~.^2 - 1,as.data.frame(X))
  X_int = X_int[,(p+1):ncol(X_int)]
  Xt = t(X)

  # --- projection matrix --- #
  P = diag(n) - X%*%solve(t(X)%*%X)%*%t(X)
  Pt = t(P)
  X_big_proj_train = rbind(X,X)
  P_big_train = diag(2*n) - X_big_proj_train%*%solve(t(X_big_proj_train)%*%X_big_proj_train)%*%t(X_big_proj_train)
  
  # --- storage --- #
  S = floor((nrun - burn)/thin)
  beta_st = gamma_st = matrix(NA,nrow = S,ncol = p)
  beta_z_st = gamma_z_st = matrix(NA,nrow = S,ncol = ncol(Z))
  lambda_st = gamma_int_st = matrix(NA,nrow = S,ncol = ncol(X_int))
  l_st = gamma_l_st = matrix(NA,nrow = S,ncol = p)
  rho_st = sigma_sq_st = numeric(S)
  y_hat = y_hat2 = y_hat3 =  matrix(NA,nrow = S, ncol = n)
  mu_st = matrix(NA,nrow = S, ncol = nrow(X))
  Omega_st = array(0,c(S,p,p))
  Omega_01 = array(0,c(S,p,p))
  Lambda_st = array(0,c(S,p,k_na))
  
  # --- initial values --- #
  pi_prob = 0.5; pi_int = 0.1
  gamma = rbinom(p,1,pi_prob)
  gamma_int = rbinom(ncol(X_int),1,pi_int)
  G = matrix(0,p,p)
  beta = rep(0,p)
  lambda = rep(0,ncol(X_int))
  beta_z = rep(0,ncol(Z))
  gamma_z = rep(0,ncol(Z))
  sigmasq = 1
  l = rep(1,p)
  gamma_l = rep(1,p)
  X_sweep = sweep(X, 2, l, `*`)
  pi_prob = 0.5
  
  # --- initial value for rho --- #
  if(exp_C == 1){
    rho = 4
  }else if(exp_C == 2){
    rho = 1
  }else{
    rho = 1
  }

  # --- covariance function --- #
  C = rho*fields::Exp.cov(X_sweep, x2=NULL, p = exp_C)
  PCPt = P%*%C%*%Pt
  CI = PCPt + sigmasq*diag(n)
  chol = chol(CI)
  CI_inv = chol2inv(chol)
  logdetCI = 2*as.numeric(sum(log((diag(chol)))))
  
  at = ceiling(nrun/100)
  if(verbose) {
    pb = txtProgressBar(style = 3)
  }  
  
  for (s in 1:nrun){
    
    # --- Update missing data --- #
    if(na){
      
      # update missing data
      W_pred = eta_na %*% t(Lambda_na) + mvtnorm::rmvnorm(n, sigma = diag(1/ps_na))
      W[W_na] = W_pred[W_na]
      XZ[,list_na] = W
      X = XZ[,1:p]
      Z = XZ[,(p+1):ncol(XZ)]
      
      # update Sigma 
      Wtil = W - eta_na %*% t(Lambda_na) 
      ps_na = rgamma(p_na, 1 + 0.5*n, 1+0.5*colSums(Wtil^2))
      Sigma_na = diag(1/ps_na)
      
      # update Lambda
      eta2 = t(eta_na) %*% eta_na
      for (j in 1:p_na) {
        Qlam = diag(Plam_na[j, ]) + ps_na[j] * eta2
        blam = ps_na[j] * (t(eta_na) %*% W[, j])
        Llam = t(chol(Qlam))
        zlam = rnorm(k_na)
        vlam = solve(Llam, blam)
        mlam = solve(t(Llam), vlam)
        ylam = solve(t(Llam), zlam)
        Lambda_na[j, ] = t(ylam + mlam)
      }
      
      # update Eta
      Lmsg = Lambda_na * ps_na
      Veta1 = diag(k_na) + t(Lmsg) %*% Lambda_na
      eigs = eigen(Veta1)
      if(all(eigs$values > 1e-6)) {
        Tmat = sqrt(eigs$values) * t(eigs$vectors)
      } else {
        Tmat = chol(Veta1)
      }
      R = qr.R(qr(Tmat))
      S = solve(R)
      Veta = S %*% t(S)                                              
      Meta = W %*% Lmsg %*% Veta                                    
      eta_na = Meta + matrix(rnorm(n*k_na), nrow = n, ncol = k_na) %*% t(S)  
      
    }
    
    if(lod){
      
      for (i in 1:nrow(X_lod)) {
        for(j in 1:ncol(X_lod)){
          if(X_lod[i,j]){
            X[i, j] = truncnorm::rtruncnorm(n = 1, a = -Inf, 
                                            b = lod_vec[j],
                                            mean = Lambda_na[j, ] %*% eta_na[i, ], 
                                            sd = sqrt(1/ps_na[j]) )
          }
        }
      }
    }
    
    if(na){
      
      # update interactions matrix
      X_int = model.matrix(y~.^2 - 1,as.data.frame(X))
      X_int = X_int[,(p+1):ncol(X_int)]
      X_int_test = model.matrix(y_test~.^2 - 1,as.data.frame(X_test))
      X_int_test = X_int_test[,(p+1):ncol(X_int_test)]
      
      # update projection matrix
      P = diag(n) - X%*%solve(t(X)%*%X)%*%t(X)
      Pt = t(P)
      X_big_proj_train = rbind(X,X)
      P_big_train = diag(2*n) - X_big_proj_train%*%solve(t(X_big_proj_train)%*%X_big_proj_train)%*%t(X_big_proj_train)
      
    }
    
    # --- update gamma --- #
    
    set = sample(1:p,p)
    for (j in set){
      csi = y - (X_int%*%lambda + Z%*%beta_z)
      
      # j-th variable NOT in the model
      gamma[j] = 0
      ind = which(gamma != 0)
      if(length(ind) == 0){
        logpi_0 = 0
      }else{
        X_0 = X[,ind]
        A_delta_0_inv = t(X_0)%*%CI_inv%*%X_0 + diag(length(ind))
        m_0 = t(X_0)%*%CI_inv%*%csi
        S_0 = t(m_0)%*%solve(A_delta_0_inv)%*%m_0
        det0 = determinant(A_delta_0_inv,log = T)
        logpi_0 = -0.5*det0$modulus[1] +
          0.5*S_0 -
          0.5*log(det(diag(length(ind))))
        
        logpi_0 = as.numeric(logpi_0)
        
      }
      
      # j-th variable in the model
      gamma[j] = 1
      ind = which(gamma != 0)
      X_1 = X[,ind]
      A_delta_1_inv = t(X_1)%*%CI_inv%*%X_1 + diag(length(ind))
      m_1 = t(X_1)%*%CI_inv%*%csi
      S_1 = t(m_1)%*%solve(A_delta_1_inv)%*%m_1
      det1 = determinant(A_delta_1_inv,log = T)
      
      logpi_1 = -0.5*det1$modulus[1] + 
        0.5*S_1 -
        0.5*log(det(diag(length(ind))))
      logpi_1 = as.numeric(logpi_1)
      
      R_j = exp(logpi_0 - logpi_1)
      prob_inv = 1 + R_j*(1-pi_prob)/pi_prob
      gamma[j] = rbinom(1,1,1/prob_inv)
      
      if(gamma[j] == 0){
        beta[j] = 0
      }
      
      # set to zero the lambdas according  to heredity  condition
      G = matrix(0,p,p)
      if(heredity == "strong"){
        
        G = gamma%*%t(gamma)
        
      }else if(heredity == "weak"){
        
        for(h in 1:p){
          if(gamma[h] != 0){
            G[h,] = 1
            G[,h] = 1
          }
        }
      }else{
        stop("provide proper heredity condition")
      }
      H = G[lower.tri(G)]
      ind = which(H == 0)
      lambda[ind] = 0
    }
    
    # --- sample beta_j diff from zero ---
    csi = y - (X_int%*%lambda + Z%*%beta_z)
    ind = which(gamma != 0)
    if(length(ind) > 0){
      X_1 = X[,ind]
      D_inv = diag(length(ind))
      V = solve(t(X_1)%*%CI_inv%*%X_1 + D_inv)
      m = V%*%(t(X_1)%*%CI_inv%*%csi)
      beta[ind] = as.vector(bayesSurv::rMVNorm(n = 1, mean = m, Sigma = V))
      beta[-ind] = 0
    }
    
    
    # --- update pi --- #
    d1 = sum(gamma)
    pi_prob = rbeta(1,1+d1,20+p-d1)
    
    
    # --- update beta_z --- #
    D_inv = diag(1/(gamma_z*(tau)^2 + (1-gamma_z)*(tau*c)^2),nrow = length(beta_z))
    V = solve(t(Z)%*%CI_inv%*%Z + D_inv)
    csi = y - (X_int%*%lambda + X%*%beta)
    m = V%*%(t(Z)%*%CI_inv%*%csi)
    beta_z = as.vector(bayesSurv::rMVNorm(n = 1, mean = m, Sigma = V))
    
    
    # --- update gamma_z --- with continuos spike#
    for (j in 1:length(beta_z)){
      b = dnorm(beta_z[j],0,c*tau)
      a = dnorm(beta_z[j],0,tau)
      gamma_z[j] = rbinom(1,1,a/(a+b))
    }
    
    # --- update lambda --- #
    G = matrix(0,p,p)
    if(heredity == "strong"){
      
      G = gamma%*%t(gamma)
      
    }else if(heredity == "weak"){
      
      for(h in 1:p){
        if(gamma[h] != 0){
          G[h,] = 1
          G[,h] = 1
        }
      }
    }else{
      stop("provide proper heredity condition ")
    }
  
    
    if(lambda_prior == "cont_spike"){
      # continuous spike
      H = G[lower.tri(G)]
      ind = which(H != 0)
      
      if (length(ind)>0){
        lambda[-ind] = 0 
        X_int_curr = X_int[,ind]    
        E_inv = diag(1/(gamma_int[ind]*(tau)^2 + (1-gamma_int[ind])*(c*tau)^2),
                     nrow = length(ind))
        V = solve(t(X_int_curr)%*%CI_inv%*%X_int_curr + E_inv)
        csi = y - (Z%*%beta_z + X%*%beta)
        m = V%*%(t(X_int_curr)%*%CI_inv%*%csi)
        lambda[ind] = as.vector(bayesSurv::rMVNorm(n = 1, mean = m, Sigma = V))
      }else{
        lambda = rep(0,length(gamma_int))
      }
      
      # --- update gamma_int ---
      for (j in 1:length(lambda)){
        if (H[j] == 0){
          gamma_int[j] = 0
          lambda[j] = 0
        }else{
          b = dnorm(lambda[j],0,(c*tau))
          a = dnorm(lambda[j],0,tau)
          gamma_int[j] = rbinom(1,1,a/(a+b))
        }
      }
      
      
    }else{
      # discrete spike
      # very slow especially with weak heredity
      H = G[lower.tri(G)]
      ind_int = which(H != 0)
      lambda[-ind_int] = 0
      gamma_int[-ind_int] = 0
      csi = y - (X%*%beta + Z%*%beta_z)
      
      #ind_int = sample(ind_int, length(ind_int))
      if (length(ind_int)>0){
        for (j in ind_int){
          
          # j-th interaction variable not in the model
          gamma_int[j] = 0
          ind = which(gamma_int != 0)
          if(length(ind)==0){
            logpi_0 = 0
          }else{
            X_0 = X_int[,ind]
            A_delta_0_inv = t(X_0)%*%CI_inv%*%X_0 + diag(length(ind))
            m_0 = t(X_0)%*%CI_inv%*%csi
            S_0 = t(m_0)%*%solve(A_delta_0_inv)%*%m_0
            det0 = determinant(A_delta_0_inv,log = T)
            logpi_0 = -0.5*det0$modulus[1] +
              0.5*S_0-
              0.5*log(det(diag(length(ind))))
            logpi_0 = as.numeric(logpi_0)
          }
          
          
          # j-th variable in the model
          gamma_int[j] = 1
          ind = which(gamma_int != 0)
          X_1 = X_int[,ind]
          A_delta_1_inv = t(X_1)%*%CI_inv%*%X_1 + diag(length(ind))
          m_1 = t(X_1)%*%CI_inv%*%csi
          S_1 = t(m_1)%*%solve(A_delta_1_inv)%*%m_1
          det1 = determinant(A_delta_1_inv,log = T)
          logpi_1 = -0.5*det1$modulus[1] +
            0.5*S_1-
            0.5*log(det(diag(length(ind))))
          logpi_1 = as.numeric(logpi_1)
          
          
          R_j = exp(logpi_0 - logpi_1)
          prob_inv = 1+ R_j*(1-pi_int)/pi_int
          gamma_int[j] = rbinom(1,1,1/prob_inv)
          
          if(gamma_int[j] == 0){
            lambda[j] = 0
          }
        }
      }else{
        lambda = rep(0,length(gamma_int))
        gamma_int = rep(0,length(gamma_int))
      }
      
      # --- update lambda_j diff from zero ---
      
      ind = which(gamma_int != 0)
      if(length(ind) > 0){
        X_1 = X_int[,ind]
        D_inv = diag(length(ind))
        V = solve(t(X_1)%*%CI_inv%*%X_1 + D_inv)
        m = V%*%(t(X_1)%*%CI_inv%*%csi)
        lambda[ind] = as.vector(bayesSurv::rMVNorm(n = 1, mean = m, Sigma = V))
      }
    }
    
    
    
    # --- update pi_int ---
    d2 = sum(gamma_int)
    pi_int = rbeta(1,1+d2,5+ncol(X_int)-d2)
    
    
    # mean vector for all the updates later
    mu = y - (Z%*%beta_z + X%*%beta + X_int%*%lambda)
    
    # --- updateL's: one add or remove move ---
    
    if(rho == 0){
      
      l = rep(0,p)
      gamma_l = rep(0,p)
      X_sweep = sweep(X, 2, l, `*`)
      C = matrix(0,n,n)
      PCPt = matrix(0,n,n)
      CI = diag(n)*sigmasq
      chol = chol(CI)
      logdetCI = 2*as.numeric(sum(log((diag(chol)))))
      CI_inv = diag(n)*(1/sigmasq)
      
    }else{
      set = sample(1:p,1)
      for (j in set){
        
        if (gamma_l[j] == 0){
          #propose from the prior so that log prop cancels out with the 
          l_starj = rgamma(1,1)
          gamma_star = gamma_l
          gamma_star[j] = 1
          l_star = l; l_star[j] = l_starj
          X_sweep_star = sweep(X, 2, l_star, `*`) 
          
        }else{
          l_starj = 0
          gamma_star = gamma_l
          gamma_star[j] = 0
          l_star = l; l_star[j] = l_starj
          X_sweep_star = sweep(X, 2, l_star, `*`)
          
          logdet_prior = - log(9/10) + log(1/10)
        }
        
        C_star = rho*fields::Exp.cov(X_sweep_star, x2=NULL, p = exp_C)
        PCPt_star = P%*%C_star%*%Pt
        CI_star = PCPt_star + diag(n)*sigmasq
        
        chol = chol(CI_star)
        logdetCI_star = 2*as.numeric(sum(log((diag(chol)))))
        CI_inv_star = chol2inv(chol)
        
        logr = logdetCI_star - logdetCI +
          t(mu)%*%(CI_inv_star - CI_inv)%*%mu
        logr = -0.5*logr 
        
        logu = log(runif(1))
        
        if (logr > logu){
          l = l_star
          X_sweep = X_sweep_star
          gamma_l = gamma_star
          C = (C_star + t(C_star))/2
          CI = (CI_star + t(CI_star))/2
          CI_inv = (CI_inv_star + t(CI_inv_star))/2
          logdetCI = logdetCI_star
          PCPt = PCPt_star
        }
      }
    }
    
    
    # ---- update L's s.t. gamma_l is != 0 ---
    ind = which(l != 0)
    
    if (length(ind) > 0){
      
      if (length(ind) == 1){
        set = ind
      }else{
        set = sample(c(ind),1)
      }
      
      
      for (j in set){
        
        l_starj = rgamma(1,l[j])
        l_star = l; l_star[j] = l_starj
        X_sweep_star = sweep(X, 2, l_star, `*`)
        
        C_star = rho*fields::Exp.cov(X_sweep_star, x2=NULL, p = exp_C)
        PCPt_star = P%*%C_star%*%Pt
        CI_star = PCPt_star + diag(n)*sigmasq
        
        chol = chol(CI_star)
        logdetCI_star = 2*as.numeric(sum(log((diag(chol)))))
        CI_inv_star = chol2inv(chol)
        
        loglik = logdetCI_star - logdetCI +
          t(mu)%*%(CI_inv_star - CI_inv)%*%mu
        loglik = -0.5*loglik
        
        logprop = dgamma(l[j],l_starj,log = T) - dgamma(l_starj,l[j],log = T)
        logprior =  0.5*(dgamma(l_starj,1,log = T) - dgamma(l[j],1,log = T))
        
        logr = loglik + logprop + logprior
        logu = log(runif(1))
        
        if (logr > logu){
          l = l_star
          C = (C_star + t(C_star))/2
          X_sweep = X_sweep_star
          CI = (CI_star + t(CI_star))/2
          CI_inv = (t(CI_inv_star) + CI_inv_star)/2
          logdetCI = logdetCI_star
          PCPt = PCPt_star
        }
        
      }
      
    }
    
    
    
    
    
    
    # --- Update rho --- #

    if (rho == 0){
      
      rho_star = rgamma(1,1)
      jstar = sample(1:p,1)
      lj_star = rgamma(1,1)
      l_star = rep(0,p)
      l_star[jstar] = lj_star
      X_sweep_star = sweep(X, 2, l_star, `*`)
      
      # covariance function
      C_star = rho_star*fields::Exp.cov(X_sweep_star, x2=NULL, p = exp_C)
      PCPt_star = P%*%C_star%*%Pt
      CI_star = PCPt_star + sigmasq*diag(n)
      chol = chol(CI_star)
      CI_inv_star = chol2inv(chol)
      logdetCI_star = 2*as.numeric(sum(log((diag(chol)))))
      
      
      loglik = logdetCI_star - logdetCI +
        t(mu)%*%(CI_inv_star - CI_inv)%*%mu
      loglik = -0.5*loglik
      
      logr = loglik
      logu = log(runif(1))
      
      if (logr > logu){
        
        rho = rho_star
        CI = (t(CI_star) + CI_star)/2
        CI_inv = (CI_inv_star + t(CI_inv_star))/2
        logdetCI = logdetCI_star
        PCPt = PCPt_star
        
        l = l_star
        X_sweep = X_sweep_star
        gamma_l = rep(0,p)
        gamma_l[jstar] = 1
        
      }else{
        
        l = rep(0,p)
        gamma_l = rep(0,p)
        X_sweep = sweep(X, 2, l, `*`)
        
      }
      
      
    }else{
      
      rho_star = 0
      
      PCPt_star = matrix(0,n,n)
      CI_star = diag(n)*sigmasq
      chol = chol(CI_star)
      logdetCI_star = 2*as.numeric(sum(log((diag(chol)))))
      CI_inv_star = chol2inv(chol)
      
      loglik = logdetCI_star - logdetCI +
        t(mu)%*%(CI_inv_star - CI_inv)%*%mu
      loglik = -0.5*loglik
      
      logr = loglik
      logu = log(runif(1))
      
      if (logr > logu){
        
        rho = rho_star
        CI = CI_star
        CI_inv = (CI_inv_star + t(CI_inv_star))/2
        logdetCI = logdetCI_star
        PCPt = PCPt_star
        l = rep(0,p)
        gamma_l = rep(0,p)
        X_sweep = sweep(X, 2, l, `*`)
        
      }
    }
    
    
    # sampling step for rho when is != 0 
    
    if (rho != 0){
      rho_jump = .1
      rho_star = rgamma(1, shape = rho^2/rho_jump^2,
                        rate = rho/rho_jump^2)
      CI_star = (rho_star/rho)*PCPt + diag(n)*sigmasq
      chol = chol(CI_star)
      logdetCI_star = 2*as.numeric(sum(log((diag(chol)))))
      CI_inv_star = chol2inv(chol)
      
      
      loglik = logdetCI_star - logdetCI +
        t(mu)%*%(CI_inv_star - CI_inv)%*%mu
      loglik = -0.5*loglik
      
      logprop = dgamma(rho,shape = rho_star^2/rho_jump^2,
                       rate = rho_star/rho_jump^2,log = T) -
        dgamma(rho_star,shape = rho^2/rho_jump^2,
               rate = rho/rho_jump^2,log = T)
      logprior = dgamma(rho_star,1,log = T) -
        dgamma(rho,1,log = T)
      
      logr = loglik + logprop + logprior
      logu = log(runif(1))
      
      if (logr > logu){
        rho = rho_star
        CI = (t(CI_star) + CI_star)/2
        CI_inv = (CI_inv_star + t(CI_inv_star))/2
        logdetCI = logdetCI_star
        PCPt = P%*%CI%*%Pt
      }
      
    }
    
    
    
    # --- update sigmasq --- #
    
    sigmasq_jump = 0.1
    sigmasq_star = rgamma(1, shape = sigmasq^2/sigmasq_jump^2, 
                          rate = sigmasq/sigmasq_jump^2)
    CI_star = PCPt + diag(n)*sigmasq_star
    chol = chol(CI_star)
    logdetCI_star = 2*as.numeric(sum(log((diag(chol)))))
    CI_inv_star = chol2inv(chol)
    
    
    loglik = logdetCI_star - logdetCI +
      t(mu)%*%(CI_inv_star - CI_inv)%*%mu
    loglik = -0.5*loglik
    
    logprop = dgamma(sigmasq,shape = sigmasq_star^2/sigmasq_jump^2, 
                     rate = sigmasq_star/sigmasq_jump^2,log = T) - 
      dgamma(sigmasq_star,shape = sigmasq^2/sigmasq_jump^2, 
             rate = sigmasq/sigmasq_jump^2,log = T)
    logprior = dgamma(sigmasq_star,1,log = T) - 
      dgamma(sigmasq,1,log = T)
    
    logr = loglik + logprop + logprior
    logu = log(runif(1))
    
    if (logr > logu){
      sigmasq = sigmasq_star
      CI = (t(CI_star) + CI_star)/2
      CI_inv = (CI_inv_star + t(CI_inv_star))/2
      logdetCI = logdetCI_star
    }
    
    
    if(verbose & (s %% at == 0)) setTxtProgressBar(pb, s / nrun)
    
    if(s > burn & s%%thin == 0){
      
      # --- storage of coeffs ---
      beta_st[count,] = beta
      gamma_st[count,] = gamma
      beta_z_st[count,] = beta_z
      gamma_z_st[count,] = gamma_z
      lambda_st[count,] = lambda
      gamma_int_st[count,] = gamma_int
      rho_st[count] = rho
      sigma_sq_st[count] = sigmasq
      l_st[count,] = l
      
      Omega_curr = matrix(0,p,p)
      Omega_curr[lower.tri(Omega_curr)] = lambda/2
      Omega_st[count,,] = Omega_curr + t(Omega_curr)
      
      Omega_curr = matrix(0,p,p)
      Omega_curr[lower.tri(Omega_curr)] = gamma_int
      Omega_01[count,,] = Omega_curr + t(Omega_curr)
      
      
      # in sample prediction
      K_big_train = rho*fields::Exp.cov(rbind(X_sweep,X_sweep), p = exp_C)
      PKPt_big_train = P_big_train%*%K_big_train%*%t(P_big_train) + sigmasq*diag(2*n)
      Sig_22_inv_train = solve(PKPt_big_train[(n+1):(2*n),(n+1):(2*n)])
      
          
      mu = Z%*%beta_z + X%*%beta + X_int%*%lambda
      y_hat[count,] = mu + PKPt_12_train%*%Sig_22_inv_train%*%(y-mu)

      count = count + 1
      
    }

    
  }
  
  ret_list = list(beta = beta_st,
                  gamma_beta = gamma_st,
                  beta_z = beta_z_st,
                  gamma_z = gamma_z_st,
                  lambda = lambda_st,
                  gamma_int = gamma_int_st,
                  sigmasq = sigma_sq_st,
                  rho = rho_st,
                  exp_C = exp_C,
                  l = l_st,
                  Omega = Omega_st,
                  Omega_01 = Omega_01,
                  y_hat = y_hat,
                  X = X)
  
  if(na){
    ret_list[["Lambda"]] = Lambda_st
  }
  
  return(ret_list)
  
  
}

