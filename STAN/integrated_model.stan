data{
  int nsite;
  int npixel;
  int npar_det;
  int npar_state;
  int npar_thin; 
  int m; // total number of data points
  int CONSTANT;// 
  int W;
  matrix[npixel, npar_det] v; // detection covariates for the entire region B
  matrix[npixel, npar_state] x_s; // covariates for latent state
  matrix[npixel, npar_thin] h_s; //covariates for thinning probability
  real cell_area;//log area of grid cell 
  int pa_pixel[nsite]; //The pixels the presence only data occurred
  int po_pixel[m]; //The pixels the presence only data occurred
  int ones[m];
  int y[nsite];
  
  
  
}
parameters{
  vector[npar_det] a;  //det/non-det data model regression coefficients
  vector[npar_state] beta; //latent state model regression coefficients
  vector[npar_thin] cc; //presence-only data model regression coefficients
}

transformed parameters{
  vector<lower=0>[npixel] lambda;
  vector<lower=0,upper=1>[npixel] psi; // probability of Species presence in a gridcell 
  vector[npixel] b; //  presence only thinning prob linear predictor
  real po_denominator; //The presence_only data model denominator 
  vector[nsite] rho; // detection/non-detection data model linear predictor
  vector[m] po_p;
  vector[nsite] log_lik_occ;
  
  //   # Bayesian version of the Koshkina (2017) model.
  // #
  // # The latent-state model
  // for(pixel in 1:npixel){
  //   lambda[pixel] = exp(dot_product(x_s[pixel,], beta) + cell_area);
  //   psi[pixel] = 1 - exp(-lambda[pixel]);
  //   // # presence only thinning prob linear predictor
  //   b[pixel] = inv_logit(dot_product(h_s[pixel,] , cc));
  // }
    lambda = exp(x_s * beta + cell_area);
    psi = 1 - exp(-lambda);
    // # presence only thinning prob linear predictor
    b = inv_logit(h_s * cc);
  //  # The presence_only data model denominator, which
  // #  is the thinned poisson process across the
  // #  whole region (divided by the total number of 
  // #  data points because it has to be 
  // #  evaluated for each data point).
  // po_denominator = dot_product(to_row_vector(lambda) , to_row_vector(b)) / m;
  po_denominator = dot_product(lambda , b) / m;
  // po_denominator = to_row_vector(lambda) * (b) / m;

  // for(po in 1:m){
  //   po_p[po] = exp( log(lambda[po_pixel[po]] * b[po_pixel[po]] ) -
  //         log(po_denominator) 
  //         ) / CONSTANT ;
  // }
  po_p = exp( log(lambda[po_pixel] .* b[po_pixel] ) -
          log(po_denominator) 
          ) / CONSTANT ;
      // # detection/non-detection data model linear predictor
  // for(site in 1:nsite)  rho[site] = inv_logit(dot_product(v[pa_pixel[site], ], a));
   rho = inv_logit(v[pa_pixel, ]* a);
  
  // # The number of detections for site is a binomial
// #  process with Pr(rho[site]) | z = 1 with
// #  W sampling occasions per site.
for(site in 1:nsite){
// https://mbjoseph.github.io/posts/2020-04-28-a-step-by-step-guide-to-marginalizing-over-discrete-parameters-for-ecologists-using-stan/
  if (y[site] > 0) {
      log_lik_occ[site] = log(psi[pa_pixel[site]]) + binomial_lpmf(y[site] | W, rho[site]);
    } else {
      log_lik_occ[site] = log_sum_exp(
        log(psi[pa_pixel[site]]) + binomial_lpmf(y[site] | W, rho[site]), 
        log1m(psi[pa_pixel[site]])
      );
    }
} 
}


model{
  // # Priors for latent state model
  to_vector(beta) ~ normal(0, 1);
  // # Priors for presence-only data model
  to_vector(cc) ~ logistic(0, 1);
  // # Priors for det/non-det data model
  to_vector(a) ~ logistic(0, 1);

  // # Loop through each presence-only datapoint
  // #  using Bernoulli one's trick. The numerator
  // #  is just the thinned poisson process for
  // #  the ith data point.
  target += bernoulli_lpmf(ones | po_p);
  target += log_lik_occ;
}


generated quantities{
  // # Derived parameter, the number of cells occupied
  real zsum;
  vector[npixel] z;
  for (p in 1:npixel) z[p] = bernoulli_rng(psi[p]);
  zsum = sum(z);

}


