MCMCfunc <- function(Y, E, W, n.rep, num.C, num.D, time, rho, tau, theta.C, 
                     lambda, sigma, theta.D, normal.prior.var, 
                     gamma.prior.scale, gamma.prior.shape){

n <- ncol(W)
n.time <- ncol(Y)
  
## Set initial parameter values
slope <- rep(NA, nrow(Y))
intercept <- rep(NA, nrow(Y))     
 for(i in 1:n)
 {
   mod <- glm(Y[i, ]~offset(log(E[i ,])) + time, family="poisson")
   intercept[i] <- mod$coefficients[1]
   slope[i] <- mod$coefficients[2]
 }

kmean <- kmeans(x=intercept, centers=num.C, nstart=1000)     
alpha.temp <- kmean$centers 
alpha.order <- order(alpha.temp)
cluster.temp <- kmean$cluster
alpha <- sort(alpha.temp)
C <- rep(NA, nrow(Y))
     for(i in 1:nrow(Y))
     {
      C[i] <- which(alpha.temp[cluster.temp[i]]==alpha)    
     }
     
kmean <- kmeans(x=slope, centers=num.D, nstart=1000)     
beta.temp <- kmean$centers 
beta.order <- order(beta.temp)
cluster.temp <- kmean$cluster
beta <- sort(beta.temp)
D <- rep(NA, nrow(Y))
     for(i in 1:nrow(Y))
     {
      D[i] <- which(beta.temp[cluster.temp[i]]==beta)    
     }

phi <- rep(0, nrow(Y))
delta <- rep(0, nrow(Y))
     

     
###Set initial parameters
CC <- C
alpha.list <- alpha[CC]
C.bar <- floor((num.C+1)/2)
centres.C <- ((1:num.C)-C.bar)^2


DD <- D
beta.list <- beta[DD]
D.bar <- floor((num.D+1)/2)
centres.D <- ((1:num.D)-D.bar)^2

logE <- log(E)


###Setting W as a double
n.neighbours <- as.numeric(apply(W, 1, sum))
W.duplet <- c(NA, NA)
     for(i in 1:n)
     {
          for(j in 1:n)
          {
               if(W[i,j]==1)
               {
               W.duplet <- rbind(W.duplet, c(i,j))     
               }else{}
          }
     }
W.duplet <- W.duplet[-1, ]     
n.duplet <- nrow(W.duplet) 


## Create the list object
Wlist <- as.list(rep(NA,n))     
     for(i in 1:n)
     {
     Wlist[[i]] <- which(W[i, ]==1)     
     }

     


###Create stores###
alpha.store <- matrix(nrow=n.rep, ncol=num.C)
phi.store <- matrix(nrow=n.rep,ncol=n)
tau.store <- rep(0,n.rep)
rho.store <- rep(0,n.rep)
C.store <- matrix(nrow=n.rep, ncol=n)
theta.C.store <- rep(0,n.rep)

beta.store <- matrix(nrow=n.rep, ncol=num.D)
delta.store <- matrix(nrow=n.rep,ncol=n)
sigma.store <- rep(0,n.rep)
lambda.store <- rep(0,n.rep)
D.store <- matrix(nrow=n.rep, ncol=n)
theta.D.store <- rep(0,n.rep)



###Create acceptance stores###
accept.alpha <- c(0,0)
accept.phi <- c(0,0)
accept.rho <- c(0,0)
accept.beta <- c(0,0)
accept.delta <- c(0,0)
accept.lambda <- c(0,0)
accept.C <- c(0,0)
accept.D <- c(0,0)
accept.theta.C <- c(0,0)
accept.theta.D <- c(0,0)

accept.alpha.all <- c(0,0)
accept.phi.all <- c(0,0)
accept.rho.all <- c(0,0)
accept.beta.all <- c(0,0)
accept.delta.all <- c(0,0)
accept.lambda.all <- c(0,0)
     
alphapropvar <- 0.1     
betapropvar <- 0.1   
deltapropvar <- 0.1   
phipropvar <- 0.1   
rhopropvar <- 0.1   
lambdapropvar <- 0.1  
thetapropvar <- 0.05
     
     
## Create the set of determinants     
Wstar <- diag(n.neighbours) - W
Wstar.eigen <- eigen(Wstar)
Wstar.val <- Wstar.eigen$values

Q.rho <- rho*Wstar + (1 - rho)*diag(1, n, n)
det.Q.rho <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))    

Q.lambda <- lambda*Wstar + (1 - lambda)*diag(1, n, n)
det.Q.lambda <-  0.5 * sum(log((lambda * Wstar.val + (1-lambda))))    

     
##########################################################################
##########################################################################

###Start the loop###
for(i in 1:n.rep){

if(floor(i/100)==i/100){print(i)}

###Update C###
if(num.C>1){
  
prop.C <- rep(0,n)
nums <- 1:num.C
for(j in 1:n){
  options <- nums[-CC[j]]
  prop.C[j] <- sample(x=options,size=1)
}
prop.alpha <- alpha[prop.C]

prob1 <- theta.C * (CC - C.bar)^2 - theta.C * (prop.C - C.bar)^2
offset <- E*exp(matrix(rep(phi,n.time),nrow=n)+(matrix(rep(beta.list,n.time),nrow=n)+matrix(rep(delta,n.time),nrow=n))*matrix(rep(time,n),nrow=n, byrow=TRUE))

test=Cupdate(Y, offset, prob1, CC, prop.C, alpha.list, prop.alpha, n.time, n)
CC <- test[[1]]

accept.C[1] <- accept.C[1] + test[[2]]
accept.C[2] <- accept.C[2] + n

}else{
  accept.C[1] <- accept.C[1] + n
  accept.C[2] <- accept.C[2] + n
}

alpha.list <- alpha[CC]
C.store[i,] <- CC


###Update D###
if(num.D > 1){
prop.D <- rep(0,n)
nums <- 1:num.D
for(j in 1:n){
  options <- nums[-DD[j]]
  prop.D[j] <- sample(x=options,size=1)
}
prop.beta <- beta[prop.D]

prob1 <- theta.D * (DD - D.bar)^2 - theta.D * (prop.D - D.bar)^2
offset <- E*exp(matrix(rep(alpha.list,n.time),nrow=n)+matrix(rep(phi,n.time),nrow=n)+(matrix(rep(delta,n.time),nrow=n))*matrix(rep(time,n),nrow=n, byrow=TRUE))

test=Dupdate(Y, offset, prob1, DD, prop.D, beta.list, prop.beta, n.time, n, time)

DD <- test[[1]]

accept.D[1] <- accept.D[1] + test[[2]]
accept.D[2] <- accept.D[2] + n

}else{
  accept.D[1] <- accept.D[1] + n
  accept.D[2] <- accept.D[2] + n
}

beta.list <- beta[DD]
D.store[i,] <- DD




###Update alpha###
proposal <- c(-1000, alpha, 1000)
for(j in 1:num.C)
{
  proposal[(j+1)] <- rtrunc(n=1, spec="norm", a=proposal[j], b=proposal[(j+2)], mean=proposal[(j+1)], sd=alphapropvar)    
}
prop.alpha <- proposal[2:(num.C+1)]
prop.alpha.list <- prop.alpha[CC]

offset <- E*exp(matrix(rep(phi,n.time),nrow=n)+(matrix(rep(beta.list,n.time),nrow=n)+matrix(rep(delta,n.time),nrow=n))*matrix(rep(time,n),nrow=n, byrow=TRUE))

test=clustalphaupdate(Y, offset, normal.prior.var, alpha, prop.alpha, alpha.list, prop.alpha.list, n.time, n)
alpha <- test[[1]]
accept.alpha[1] <- accept.alpha[1]+test[[2]]     
accept.alpha[2] <- accept.alpha[2]+1
alpha.store[i,] <- alpha




###Update beta###
proposal <- c(-1000, beta, 1000)
for(j in 1:num.D)
{
  proposal[(j+1)] <- rtrunc(n=1, spec="norm", a=proposal[j], b=proposal[(j+2)], mean=proposal[(j+1)], sd=betapropvar)    
}
prop.beta <- proposal[2:(num.D+1)]
prop.beta.list <- prop.beta[DD]

offset <- E*exp(matrix(rep(alpha.list,n.time),nrow=n)+matrix(rep(phi,n.time),nrow=n)+(matrix(rep(delta,n.time),nrow=n))*matrix(rep(time,n),nrow=n, byrow=TRUE))

test=clustbetaupdate(Y, offset, normal.prior.var, beta, prop.beta, beta.list, prop.beta.list, n.time, n, time)
beta <- test[[1]]
accept.beta[1] <- accept.beta[1]+test[[2]]     
accept.beta[2] <- accept.beta[2]+1
beta.store[i,] <- beta



###Update phi###
offset <- E*exp(matrix(rep(alpha.list,n.time),nrow=n)+(matrix(rep(beta.list,n.time),nrow=n)+matrix(rep(delta,n.time),nrow=n))*matrix(rep(time,n),nrow=n, byrow=TRUE))

test = clustpoissoncarupdate(Y, offset, Wlist, n.neighbours, phi, rho, tau, phipropvar, n)

phi <- test[[1]]
accept.phi[1] <- accept.phi[1]+test[[2]]     
accept.phi[2] <- accept.phi[2]+n

for(j in 1:num.C){
  phi[which(CC==j)] <- phi[which(CC==j)] - mean(phi[which(CC==j)])
}

phi.store[i,] <- phi




###Update delta###
offset <- E*exp(matrix(rep(alpha.list,n.time),nrow=n)+matrix(rep(phi,n.time),nrow=n) + matrix(rep(beta.list,n.time),nrow=n)*matrix(rep(time,n),nrow=n, byrow=TRUE))

test = clustpoissoncarupdate2(Y, offset, Wlist, n.neighbours, time,
                      delta, lambda, sigma, deltapropvar, n)

delta <- test[[1]]
accept.delta[1] <- accept.delta[1]+test[[2]]     
accept.delta[2] <- accept.delta[2]+n

for(j in 1:num.D){
  delta[which(DD==j)] <- delta[which(DD==j)] - mean(delta[which(DD==j)])
}

delta.store[i,] <- delta





##Update tau2###
tau.shape <- n/2 + gamma.prior.shape
tau.scale <- 0.5*t(phi)%*%Q.rho%*%phi + gamma.prior.scale
tau <- rinvgamma(1,tau.shape,tau.scale)

tau.store[i] <- tau




##Update sigma###
sigma.shape <- n/2 + gamma.prior.shape
sigma.scale <- 0.5*t(delta)%*%Q.lambda%*%delta + gamma.prior.scale
sigma <- rinvgamma(1,sigma.shape,sigma.scale)

sigma.store[i] <- sigma





###Update rho###
prop.rho <-  rtrunc(n=1, spec="norm", a=0, b=1, mean=rho, sd=rhopropvar)

prop.Q.rho <- prop.rho*Wstar + (1 - prop.rho)*diag(1, n, n)
prop.det.Q.rho <-  0.5 * sum(log((prop.rho * Wstar.val + (1-prop.rho))))    
full.rho <- det.Q.rho - 0.5*(t(phi)%*%Q.rho%*%phi)/tau
full.prop.rho <- prop.det.Q.rho - 0.5*(t(phi)%*%prop.Q.rho%*%phi)/tau
ratio.rho <- exp(full.prop.rho - full.rho)


if(runif(1,0,1) < ratio.rho)
{
  rho <- prop.rho
  Q.rho <- prop.Q.rho
  det.Q.rho <- prop.det.Q.rho
  accept.rho[1] <- accept.rho[1]+1
}

accept.rho[2] <- accept.rho[2]+1
rho.store[i] <- rho



###Update lambda###
prop.lambda <-  rtrunc(n=1, spec="norm", a=0, b=1, mean=lambda, sd=lambdapropvar)

prop.Q.lambda <- prop.lambda*Wstar + (1 - prop.lambda)*diag(1, n, n)
prop.det.Q.lambda <-  0.5 * sum(log((prop.lambda * Wstar.val + (1-prop.lambda))))    
full.lambda <- det.Q.lambda - 0.5*(t(delta)%*%Q.lambda%*%delta)/tau
full.prop.lambda <- prop.det.Q.lambda - 0.5*(t(delta)%*%prop.Q.lambda%*%delta)/tau
ratio.lambda <- exp(full.prop.lambda - full.lambda)


if(runif(1,0,1) < ratio.lambda)
{
  lambda <- prop.lambda
  Q.lambda <- prop.Q.lambda
  det.Q.lambda <- prop.det.Q.lambda
  accept.lambda[1] <- accept.lambda[1]+1
}

accept.lambda[2] <- accept.lambda[2]+1
lambda.store[i] <- lambda



  
###Update theta.C###
prop.theta.C <- rtrunc(n=1, spec="norm", a=1, b=100, mean=theta.C, sd=thetapropvar)   

prob1 <- sum((CC-C.bar)^2) * (theta.C - prop.theta.C)     
prob2 <- n*log(sum(exp(-theta.C * centres.C))) - n*log(sum(exp(-prop.theta.C * centres.C)))     
ratio.theta.C <- exp(prob1 + prob2)
 

if(runif(1,0,1) < ratio.theta.C)
{
  theta.C <- prop.theta.C
  accept.theta.C[1] <- accept.theta.C[1]+1
} 

accept.theta.C[2] <- accept.theta.C[2]+1  

theta.C.store[i] <- theta.C



###Update theta.D###
prop.theta.D <- rtrunc(n=1, spec="norm", a=1, b=100, mean=theta.D, sd=thetapropvar)   

prob1 <- sum((DD-D.bar)^2) * (theta.D - prop.theta.D)     
prob2 <- n*log(sum(exp(-theta.D * centres.D))) - n*log(sum(exp(-prop.theta.D * centres.D)))     
ratio.theta.D <- exp(prob1 + prob2)

if(runif(1,0,1) < ratio.theta.D)
{
  theta.D <- prop.theta.D
  accept.theta.D[1] <- accept.theta.D[1]+1
} 

accept.theta.D[2] <- accept.theta.D[2]+1  

theta.D.store[i] <- theta.D



     
#### Tune the acceptance probabiltiies
     k <- i/100
          if(ceiling(k)==floor(k))
          {
          #### Determine the acceptance probabilities
          betaACC <- 100 * accept.beta[1] / accept.beta[2]
          alphaACC <- 100 * accept.alpha[1] / accept.alpha[2]
          phiACC <- 100 * accept.phi[1] / accept.phi[2]
          deltaACC <- 100 * accept.delta[1] / accept.delta[2]
          rhoACC <- 100 * accept.rho[1] / accept.rho[2]
          lambdaACC <- 100 * accept.lambda[1] / accept.lambda[2]
          
          accept.beta.all <- accept.beta.all + accept.beta
          accept.alpha.all <- accept.alpha.all + accept.alpha
          accept.phi.all <- accept.phi.all + accept.phi
          accept.delta.all <- accept.delta.all + accept.delta
          accept.rho.all <- accept.rho.all + accept.rho
          accept.lambda.all <- accept.lambda.all + accept.lambda
           
          accept.alpha <- c(0,0)
          accept.phi <- c(0,0)
          accept.rho <- c(0,0)
          accept.beta <- c(0,0)
          accept.delta <- c(0,0)
          accept.lambda <- c(0,0)
          
                    
          #### beta tuning parameter
               if(betaACC > 70)
               {
               betapropvar <- 2 * betapropvar
               }else if(betaACC < 50)              
               {
               betapropvar <- 0.5 * betapropvar
               }else
               {
               }

               #### alpha tuning parameter
               if(alphaACC > 70)
               {
               alphapropvar <- 2 * alphapropvar
               }else if(alphaACC < 50)              
               {
               alphapropvar <- 0.5 * alphapropvar
               }else
               {
               }
          #### phi tuning parameter
               if(phiACC > 70)
               {
               phipropvar <- 2 * phipropvar
               }else if(phiACC < 50)              
               {
               phipropvar <- 0.5 * phipropvar
               }else
               {
               }
           #### delta tuning parameter
               if(deltaACC > 70)
               {
               deltapropvar <- 2 * deltapropvar
               }else if(deltaACC < 50)              
               {
               deltapropvar <- 0.5 * deltapropvar
               }else
               {
               }              
          
           #### rho tuning parameter
               if(rhoACC > 70)
               {
               rhopropvar <- min(2 * rhopropvar,0.5)
               }else if(rhoACC < 50)              
               {
               rhopropvar <- 0.5 * rhopropvar
               }else
               {
               }                         
               
           #### lambda tuning parameter
               if(lambdaACC > 70)
               {
               lambdapropvar <- min(2 * lambdapropvar,0.5)
               }else if(lambdaACC < 50)              
               {
               lambdapropvar <- 0.5 * lambdapropvar
               }else
               {
               }                         
               
               
               }else
          {   
          }
     
}

###Return the results and acceptance rates###
result <- list(alpha.store=alpha.store, phi.store=phi.store, tau.store=tau.store,
               rho.store=rho.store, C.store=C.store, theta.C.store=theta.C.store, beta.store=beta.store, 
               delta.store=delta.store, sigma.store=sigma.store, lambda.store=lambda.store, D.store=D.store, 
               theta.D.store=theta.D.store, accept.alpha=accept.alpha.all[1]/accept.alpha.all[2], 
               accept.phi=accept.phi.all[1]/accept.phi.all[2], accept.rho=accept.rho.all[1]/accept.rho.all[2], 
               accept.beta=accept.beta.all[1]/accept.beta.all[2], accept.delta=accept.delta.all[1]/accept.delta.all[2],
               accept.lambda=accept.lambda.all[1]/accept.lambda.all[2], accept.C=accept.C[1]/accept.C[2],accept.D=accept.D[1]/accept.D[2], 
               accept.theta.C=accept.theta.C[1]/accept.theta.C[2],accept.theta.D=accept.theta.D[1]/accept.theta.D[2])

return(result)
}