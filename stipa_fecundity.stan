//
// This Stan program defines a model for stipa fecundity based on tiller counts
//

// Input data
data {
  int<lower=1> N; // number of rows or data points
  vector[N] stipa_tillers; // response variable
  vector[N] avena_density; // predictor variable
  vector[N] stipa_bg_count; // predictor variable
  // what about the random effects for year and block? 
}

// Parameters to be estimated by the model
parameters {
  real intercept;
  real<lower=0> sigma;
  real b_avena;
  real b_stipa;  # should these each have separate mu and sigma?
}

// Optional section to transform parameters as needed for convergence issues
transformed parameters{
  
}


// The model to be estimated. 
model {
  // Declare objects necessary for the rest of the model
  
  // Set priors
  
  
  // Implement the model
  
}

