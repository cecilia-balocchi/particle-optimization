library(optimization)
error_safe <- function(expr){
  tryCatch(expr,
           error = function(e){
             message("An error occurred:\n", e)
             NA
           })
}

summary_hyperpar <- function(Y, X, A_block,  
                          gamma_init_A, gamma_init_B,
                          eta_input, rho_input = 0.9,
                          priorA = "ep", priorB = "ep"){
  if(priorA == "ep"){
    priorA_num = 0
  } else if(priorA == "unif"){
    priorA_num = 3
  } else {
    warning("priorA should be a string `ep` for Ewens-Pitman or `unif` for the uniform.")
  }
  if(priorB == "ep"){
    priorB_num = 0
  } else if(priorB == "unif"){
    priorB_num = 3
  } else {
    warning("priorB should be a string `ep` for Ewens-Pitman or `unif` for the uniform.")
  }
  Rcpp::sourceCpp("src/particle_summary.cpp")
  source("src/fun_likelihood.R")
  
  Xorig <- X
  Yorig <- Y
  Xmeans <- rowMeans(X)
  X <- X - Xmeans
  betas_mle <- numeric(N)
  for(i in 1:N)
  betas_mle[i] <- cov(Y[i,],X[i,])/var(X[i,])
  Y <- Y - betas_mle*Xmeans
  eta_py = eta_input
  sigma_py = 0
  rho = rho_input
  
  N <- dim(Y)[1]
  t <- dim(Y)[2]
  
  n_tr <- dim(X)[1]
  betas_mle <- numeric(n_tr)
  for(i in 1:n_tr)
    betas_mle[i] <- cov(Y[i,],X[i,])/var(X[i,])
  alphas_mle <- rowMeans(Y) - betas_mle * rowMeans(X) 
  sigmas <- numeric(n_tr)
  for(i in 1:n_tr)
    sigmas[i] <- sd(lm(Y[i,]~X[i,])$residuals)
  sigma2 <- mean(sigmas^2)
  
  mu <- mean(sigmas^2)
  v <- var(sigmas^2)
  alpha_sigma <- mu^2/v + 2
  beta_sigma <- mu*(alpha_sigma-1)
  
  K <- round(log(n_tr))
  
  tmp <- (max(alphas_mle)-min(alphas_mle))/(K+1)/2
  a1 <- tmp^2/sigma2*(1-0.8)
  a2 <- (max(abs(alphas_mle))/2)^2/sigma2 - (a1/(1-rho))
  
  tmp <- (max(betas_mle)-min(betas_mle))/(K+1)/2
  b1 <- tmp^2/sigma2*(1-0.8)
  b2 <- (max(abs(betas_mle))/2)^2/sigma2 - (b1/(1-rho))
  
  partA <- gamma_init_A
  partB <- gamma_init_B
  
  log_post2 <- function(par){
    log_post(par, priorA_num, priorB_num, partA, partB, 
             a2, b2, rho, Y,X,A_block, alpha_sigma, beta_sigma,
             eta_py, sigma_py)
  }
  
  tmp_nm <- error_safe(optim(par = c(a1,b1), fn = log_post2, method = "Nelder-Mead"))
  if(!any(is.na(tmp_nm))){
    a1_new <- tmp_nm$par[1]
    b1_new <- tmp_nm$par[2]
  } else {
    a1_new <- a1
    b1_new <- b1
  }
  
  tmp_new <- particle_summary(Y, X, A_block, 
                              gamma_init_A = partA, gamma_init_B = partB,
                              a1_input = a1_new, b1_input = b1_new, a2_input = a2, b2_input = b2, 
                              alpha_sigma_input = alpha_sigma, beta_sigma_input = beta_sigma,
                              priorA_input = priorA_num, priorB_input = priorB_num, 
                              eta_input = eta_py, rho_input = rho)

  final_hyperpar <- c(a1_new, a2, b1_new, b2, alpha_sigma, beta_sigma)
  
  return(list(adjusted = tmp_new, optim = tmp_nm, hyperpar = final_hyperpar))
}