update_sigmasq = function(params, constants) {
  
  sigmasq = params$sigmasq
  mu = params$mu; CI_inv = params$CI_inv; CI = params$CI
  PCPt = params$PCPt; logdetCI = params$logdetCI
  n = constants$n
  
  sigmasq_jump = 0.1
  sigmasq_star = rgamma(1, shape = sigmasq^2/sigmasq_jump^2, 
                        rate = sigmasq/sigmasq_jump^2)
  CI_star = PCPt + diag(n)*sigmasq_star
  chol = chol(CI_star + 0.01*diag(1, nrow = nrow(CI_star)))
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
  
  params[["sigmasq"]] = sigmasq
  params[["CI"]] = CI
  params[["CI_inv"]] = CI_inv
  params[["logdetCI"]] = logdetCI
  
  return(params)
}