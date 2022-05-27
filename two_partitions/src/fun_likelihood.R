likelihood_cl <- function(part, k, A_or_B, a1, a2, b1, b2, rho, Y,X,A_block){
  t <- dim(Y)[2]
  cluster <- part[[k]]
  nk <- length(cluster)
  E = diag(nk)
  Omega_k = matrix(0, ncol = nk, nrow= nk)
  Y_k = Y[cluster, , drop=FALSE]
  X_k = X[cluster, , drop=FALSE]
  if(A_or_B==1){
    c1 = a1
    c2 = a2
    tXX = t*E
  } else {
    c1 = b1
    c2 = b2
    tXX = diag(as.numeric(rowSums(X_k*X_k)), nrow = nk, ncol = nk)
  }
  if(nk == 1){
    Omega_k[1,1] = 1/(c1/(1-rho)+c2)
    if(A_or_B == 1){
      XtY = sum(Y_k)
    } else {
      XtY = sum(Y_k * X_k)
    }
  } else {
    O_k = matrix(1, ncol =nk, nrow = nk)
    A_block_k = A_block[cluster, cluster]
    D = diag(rowSums(A_block_k))
    M = rho * (D-A_block_k) + (1-rho)*diag(nk)
    Omega_k = (M/c1) - ((1-rho)/c1) * ((1-rho)/c1) * O_k / ( nk*(1-rho)/c1 + 1/c2)
    if(A_or_B == 1){
      XtY = rowSums(Y_k)
    } else {
      XtY = rowSums(Y_k * X_k)
    }
  }
  P_k = solve(tXX + Omega_k)
  Sigma_det = E - P_k * tXX
  
  quad_form = - as.numeric(t(XtY) %*% P_k %*% XtY)
  log_det = as.numeric(determinant(Sigma_det)$modulus)
  
  return(list(quad_form = quad_form, log_det = log_det))
}
likelihood_part <- function(part, A_or_B, a1, a2, b1, b2, rho, Y,X,A_block){
  log_det = 0
  quad_form = 0
  for(k in 1:length(part)){
    tmp <- likelihood_cl(part, k, A_or_B, a1,a2,b1,b2, rho, Y,X,A_block)
    log_det <- log_det + tmp$log_det
    quad_form <- quad_form + tmp$quad_form
  }
  return(list(quad_form = quad_form, log_det = log_det))
}
likelihood <- function(par, partA, partB, a2, b2, rho, Y,X,A_block, alpha_sigma, beta_sigma){
    a1 <- par[1]
    b1 <- par[2]
    if(a1 < 0 || a2 < 0 || b1 < 0 || b2 < 0){
      return(Inf)
    }
    N <- dim(Y)[1]
    t <- dim(Y)[2]
    log_det = 0
    quad_form = sum(Y*Y)
    tmp <- likelihood_part(partA, 1, a1, a2, b1, b2, rho, Y,X,A_block)
    log_det <- log_det + tmp$log_det
    quad_form <- quad_form + tmp$quad_form
    tmp <- likelihood_part(partB, 0, a1, a2, b1, b2, rho, Y,X,A_block)
    log_det <- log_det + tmp$log_det
    quad_form <- quad_form + tmp$quad_form
    log_like = lgamma(alpha_sigma +  ( N * t/2)) + alpha_sigma * log(beta_sigma)
    if(beta_sigma + quad_form/2 < 0){
      cat("negative log\n")
    }
    log_like = log_like - lgamma(alpha_sigma) - (alpha_sigma + ( N * t/2)) * log(beta_sigma + quad_form/2)
    log_like = log_like + 0.5 * log_det
    return(-log_like)
  }

log_pyp_prior <- function(part, eta, sigma){
  # for each cluster we need (1-sigma)_{n_i-1} = Gamma(n_i - sigma)/Gamma(1-sigma)
  x <- sapply(part, length)
  n <- sum(x)
  l <- sum(lgamma(x-sigma)-lgamma(1-sigma))
  l <- l + sum(log(eta + (0:(length(x)-1))*sigma))
  l <- l - sum(log(eta + 0:(n-1)))
  return(l)
}

log_post <- function(par, priorA_num, priorB_num, partA, partB, 
                     a2, b2, rho, Y,X,A_block, alpha_sigma, beta_sigma,
                     eta_py, sigma_py){
  out <- likelihood(par, partA, partB, a2, b2, rho, Y,X,A_block, alpha_sigma, beta_sigma) #likelihood returns -loglike.
  if(priorA_num == 0){
    out <- out - log_pyp_prior(partA, eta_py, sigma_py)
  }
  if(priorB_num == 0){
    out <- out - log_pyp_prior(partB, eta_py, sigma_py)
  }
  return(out)
}
