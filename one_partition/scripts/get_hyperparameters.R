get_hyperparameters <- function(Y, K, rho){
  N <- nrow(Y)
  T <- ncol(Y)
  
  alpha_mle <- rowMeans(Y)
  sigma2_mle <- (T-1)/T * apply(Y, MAR = 1, FUN = var)
  
  mean_sigma2 <- mean(sigma2_mle)
  var_sigma2 <- var(sigma2_mle)
  
  nu_sigma <- 2* ((mean_sigma2 * mean_sigma2)/var_sigma2 + 2)
  lambda_sigma <- mean_sigma2 * (nu_sigma-2)/nu_sigma
  
  sigma2_est <- mean(sigma2_mle)
  
  a1 <-  (max(alpha_mle) - min(alpha_mle))^2/(4 * (K+1)^2 * mean_sigma2 /(1 - rho))
  a2 <- max(alpha_mle)^2/(4 * mean_sigma2) - a1/(1 - rho)
  if( (a2 < 0)) hyper_parameters <- NULL
  else hyper_params <- list(a1 = a1, a2 = a2, nu_sigma = nu_sigma, lambda_sigma = lambda_sigma)
  return(hyper_params)
}