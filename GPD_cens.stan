 functions {
real gpdcens_lpdf(vector x, vector rightcens, vector slow, real sigma, real xi, real xmax)
  {
   // Declare variables first
  real loglik;
  loglik = 0.0;
  if(xi < 0){
    if(xmax > -sigma/xi){
	loglik = log(0.0);
  	return loglik;
  	}
  }		
  for(i in 1:num_elements(x)){
  	if(rightcens[i] < 0.5){
  	 loglik += - log(sigma)-(1.0+1.0/xi)*log(1.0+xi*x[i]/sigma);
  	} else{
 	 loglik += - 1.0/xi*log(1.0+xi*x[i]/sigma) ;
  	}
  	loglik += 1.0/xi*log(1.0+xi*slow[i]/sigma);
  }
  return loglik;
   }
 }


 data {
	int<lower=0> n; 	// Sample size
	vector[n] x; 		// Data vector
	vector[n] slow; 	// Lower truncation
	vector[n] rightcens; 	// Right-censoring indicator
	real xmax;
 }

 parameters {
	real<lower=0> sigma;
	real xi;
 }

 model{
 // Define prior for hyperparameters
	xi  ~ normal(0, 1);
	sigma ~ cauchy(0,5); //half-normal prior, avoids normalization
	 x ~ gpdcens(rightcens, slow, sigma, xi, xmax);
	 //target += gpdcens_lpdf(x | rightcens, slow, sigma, xi, xmax);
 }
