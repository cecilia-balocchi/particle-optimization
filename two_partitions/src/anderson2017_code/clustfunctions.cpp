  #include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Cupdate(NumericMatrix Y, NumericMatrix offset, NumericVector prob1, 
             NumericVector C, NumericVector propC, NumericVector alpha, 
             NumericVector propalpha, const int ntime, const int n)
{

double logaccept=0, accept, accepted=0;
double prob2, prob3;

//Update each C value in turn
for(int j = 0; j < n; j++){
   
   //Compute the acceptance probability
   prob2 = sum(Y( j, _) * (propalpha[j]-alpha[j]));
   prob3 = sum(offset( j, _) * (exp(alpha[j])-exp(propalpha[j])));
   logaccept = prob1[j]+prob2+prob3;
    
   //Accept or not
   accept = exp(logaccept);
   if(runif(1)[0] <= accept) 
          {
          C[j] = propC[j];
          accepted = accepted + 1;
          }
          else
          { 
          }  
}





List out(3);
out[0] = C;
out[1] = accepted;
out[2] = logaccept;
return out;
}




// [[Rcpp::export]]
List Dupdate(NumericMatrix Y, NumericMatrix offset, NumericVector prob1,
             NumericVector D, NumericVector propD, NumericVector beta, 
             NumericVector propbeta, const int ntime, const int n, 
             NumericVector time)
{

double logaccept=0, accept, accepted=0;
double prob2, prob3;

//Update each D value in turn
for(int j = 0; j < n; j++){
   
   //Compute the acceptance probability
   prob2 = sum(Y( j, _) * time * (propbeta[j]-beta[j]));
   prob3 = sum(offset( j, _) * (exp(beta[j]*time)-exp(propbeta[j]*time)));
   logaccept = prob1[j]+prob2+prob3;

   //Accept or not
   accept = exp(logaccept);
   if(runif(1)[0] <= accept) 
          {
          D[j] = propD[j];
          accepted = accepted + 1;
          }
          else
          { 
          }  
}





List out(2);
out[0] = D;
out[1] = accepted;
return out;
}




// [[Rcpp::export]]
List clustalphaupdate(NumericMatrix Y, NumericMatrix offset, double normalpriorvar,
                 NumericVector alpha, NumericVector propalpha, 
                 NumericVector alphalist, NumericVector propalphalist,
                 const int ntime, const int n)
{
     

//Compute the acceptance probability
//NumericVector prob1(n), prob2(n);
double logaccept=0, accept, accepted=0;
double prob1, prob2;

for(int t = 0; t < ntime; t++){
     prob1 = sum(Y( _, t) * (propalphalist-alphalist));
     prob2 = sum(offset( _, t) * (exp(alphalist)-exp(propalphalist)));
     logaccept = logaccept + prob1+prob2;
}

accept = exp(logaccept);

//Accept or not
   if(runif(1)[0] <= accept) 
          {
          alpha = propalpha;
          accepted = 1;
          }
          else
          { 
          }


List out(2);
out[0] = alpha;
out[1] = accepted;
return out;
}









// [[Rcpp::export]]
List clustbetaupdate(NumericMatrix Y, NumericMatrix offset, double normalpriorvar,
                 NumericVector beta, NumericVector propbeta, 
                 NumericVector betalist, NumericVector propbetalist,
                 const int ntime, const int n, NumericVector time)
{
     


//Compute the acceptance probability
//NumericVector prob1(n), prob2(n);
double logaccept=0, accept, accepted=0;
double prob1, prob2;

for(int t = 0; t < ntime; t++){
     prob1 = sum(Y( _, t) * time[t] * (propbetalist-betalist));
     prob2 = sum(offset( _, t) * (exp(betalist*time[t])-exp(propbetalist*time[t])));
     logaccept = logaccept + prob1+prob2;
}

accept = exp(logaccept);

//Accept or not
   if(runif(1)[0] <= accept) 
          {
          beta = propbeta;
          accepted = 1;
          }
          else
          { 
          }


List out(2);
out[0] = beta;
out[1] = accepted;
return out;
}





// [[Rcpp::export]]
List clustpoissoncarupdate(NumericMatrix Y, NumericMatrix offset, List Wlist,
                      NumericVector nneighbours,
                      NumericVector phi, double rho, double tau2,
                      double phipropvar, const int n)
{
// Update the spatially correlated random effects 
//Create new objects
double logaccept=0, accepted=0;
double acceptance, sumphi;
double oldpriorbit, newpriorbit;
double priordenom, priormean, priorvar;
double propphi;
double lik1, lik2;
NumericVector phinew(n);
 
   
//  Update each random effect in turn
phinew = phi;
     for(int j = 0; j < n; j++)
     {
     // calculate prior mean and variance
     IntegerVector neighbourvec = Wlist[j];
     int m = neighbourvec.size();
     sumphi = 0;
          for(int l = 0; l < m; l++) sumphi += phinew[(neighbourvec[l]-1)];
      priordenom = (nneighbours[j] * rho + (1-rho));
      priorvar = tau2 / priordenom;
      priormean = rho * sumphi / priordenom; 
      
      // propose a value  
      propphi = rnorm(1, phinew[j], sqrt(priorvar*phipropvar))[0];
      
      // Accept or reject it
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
     
      lik1 = sum(offset(j,_) * (exp(phinew[j]) - exp(propphi)));
      lik2 = sum(Y(j,_) * (propphi - phinew[j]));
      
        
      logaccept = lik1 + lik2 + oldpriorbit - newpriorbit;
      acceptance = exp(logaccept);
      
      
          if(runif(1)[0] <= acceptance) 
          {
          phinew[j] = propphi;
          accepted = accepted + 1;
          }
          else
          { 
          }
    }


List out(2);
out[0] = phinew;
out[1] = accepted;
return out;
}



// [[Rcpp::export]]
List clustpoissoncarupdate2(NumericMatrix Y, NumericMatrix offset, List Wlist,
                      NumericVector nneighbours, NumericVector time,
                      NumericVector delta, double lambda, double sigma,
                      double deltapropvar, const int n)
{
// Update the spatially correlated random effects 
//Create new objects
double logaccept=0, accepted=0;
double acceptance, sumdelta;
double oldpriorbit, newpriorbit;
double priordenom, priormean, priorvar;
double propdelta;
double lik1, lik2;
NumericVector deltanew(n);
 
   
//  Update each random effect in turn
deltanew = delta;
     for(int j = 0; j < n; j++)
     {
     // calculate prior mean and variance
     IntegerVector neighbourvec = Wlist[j];
     int m = neighbourvec.size();
     sumdelta = 0;
          for(int l = 0; l < m; l++) sumdelta += deltanew[(neighbourvec[l]-1)];
      priordenom = (nneighbours[j] * lambda + (1-lambda));
      priorvar = sigma / priordenom;
      priormean = lambda * sumdelta / priordenom; 
      
      // propose a value  
      propdelta = rnorm(1, deltanew[j], sqrt(priorvar*deltapropvar))[0];
      
      // Accept or reject it
      newpriorbit = (0.5/priorvar) * pow((propdelta - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((deltanew[j] - priormean), 2);
     
      lik1 = sum(offset(j,_) * (exp(deltanew[j] * time) - exp(propdelta * time)));
      lik2 = sum(Y(j,_) * time * (propdelta - deltanew[j]));
      
        
      logaccept = lik1 + lik2 + oldpriorbit - newpriorbit;
      acceptance = exp(logaccept);
      
      
          if(runif(1)[0] <= acceptance) 
          {
          deltanew[j] = propdelta;
          accepted = accepted + 1;
          }
          else
          { 
          }
    }


List out(2);
out[0] = deltanew;
out[1] = accepted;
return out;
}

