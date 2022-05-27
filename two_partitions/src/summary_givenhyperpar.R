library(optimization)
error_safe <- function(expr){
  tryCatch(expr,
           error = function(e){
             message("An error occurred:\n", e)
             NA
           })
}

summary_givenhyperpar <- function(Y, X, A_block,  
                          hyperpar, 
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
  
  K <- round(log(N))
  
  partA <- gamma_init_A
  partB <- gamma_init_B

  a1 <- hyperpar[1]
  a2 <- hyperpar[2]
  b1 <- hyperpar[3]
  b2 <- hyperpar[4]
  alpha_sigma <- hyperpar[5]
  beta_sigma <- hyperpar[6]
  
  tmp <- particle_summary(Y, X, A_block, 
                              gamma_init_A = partA, gamma_init_B = partB,
                              a1_input = a1, b1_input = b1, a2_input = a2, b2_input = b2, 
                              alpha_sigma_input = alpha_sigma, beta_sigma_input = beta_sigma,
                              priorA_input = priorA_num, priorB_input = priorB_num, 
                              eta_input = eta_py, rho_input = rho)
  
  return(list(output = tmp, hyperpar = hyperpar))
}