// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N; // Number of plots
  int<lower = 1> S; // Number of species
  int fecundity[N]; // Fecundity of the focal species in each plot
  int size_big[N];   // Indicator variable for what size the species is
  matrix[N,S] sp_matrix; // Matrix of abundances for each species (including abundances of non-focal individuals of the focal species)
  vector[N] drought;   // Environmental values for each plot (ambient or drought)
 // int<lower = 0> Intra[S]; // Indicator boolean variable to identify the focal species (0 for non-focal and 1 for focal). Included for easier calculations
  // The below values define the regularized horseshoe priors used for species-specific parameters
  real tau0; 		// determines the scale of the global shrinkage parameter (tau)
  real slab_scale;	// scale for significant alpha_sp values
  real slab_df;		// effective degrees of freedom for significant alpha_sp values
}

transformed data{
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5*slab_df;
}

parameters{
  real lambda_0;
  real lambda_size;
  real lambda_drought;
  real alpha_generic_tilde_0;
  real alpha_generic_tilde_size;
  real alpha_generic_tilde_drought;
 // vector[2] alpha_generic_tilde; // first element is intercept, second is drought effect
 // vector[2] alpha_intra_tilde;
  vector[S] alpha_hat_ij_tilde; // non-generic intercepts
  vector[S] alpha_hat_eij_tilde; // non-generic drought effect
  vector<lower = 0>[S] local_shrinkage_ij;
  vector<lower = 0>[S] local_shrinkage_eij;
  real<lower = 0> c2_tilde;
  real<lower = 0> tau_tilde;
}

transformed parameters{
  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the normalized (and thus easier to sample)
  // 	counterparts declared in the parameters block
  real c2;
  real tau;
  vector[S] alpha_hat_ij;
  vector[S] local_shrinkage_ij_tilde;
  vector[S] alpha_hat_eij;
  vector[S] local_shrinkage_eij_tilde;
  real alpha_generic_0;
  real alpha_generic_size;
  real alpha_generic_drought;
  // vector[2] alpha_generic;
  // vector[2] alpha_intra;
  
  tau = tau0*tau_tilde; 	// tau ~ cauchy(0, tau0)
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)
  
  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  for(s in 1:S){
    local_shrinkage_ij_tilde[s] = sqrt( c2 * square(local_shrinkage_ij[s]) / (c2 + square(tau) * square(local_shrinkage_ij[s])) );
    alpha_hat_ij[s] = tau * local_shrinkage_ij_tilde[s] * alpha_hat_ij_tilde[s];
    
    local_shrinkage_eij_tilde[s] = sqrt( c2 * square(local_shrinkage_eij[s]) / (c2 + square(tau) * square(local_shrinkage_eij[s])) );
    alpha_hat_eij[s] = tau * local_shrinkage_eij_tilde[s] * alpha_hat_eij_tilde[s];
  }
  
  // scale the lambdas and alphas values
  alpha_generic_0 = 3 * alpha_generic_tilde_0 - 6;
  //  alpha_intra[1] = 3 * alpha_intra_tilde[1] - 6;
  alpha_generic_size = 0.5 * alpha_generic_tilde_size;
  alpha_generic_drought = 0.5 * alpha_generic_tilde_drought;
  //  alpha_intra[2] = 0.5 * alpha_intra_tilde[2];
}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] F_hat;
  vector[N] interaction_effects;
  matrix[N,S] alpha_eij;
  vector[N] lambda_ei;

  // set regular priors
  alpha_generic_tilde_0 ~ normal(0,1);
  alpha_generic_tilde_size ~ normal(0,1);
  alpha_generic_tilde_drought ~ normal(0,1);
//  alpha_intra_tilde ~ normal(0,1);
  lambda_0 ~ normal(0, 1);
  lambda_size ~ normal(0, 1);
  lambda_drought ~ normal(0, 1);


  // set the hierarchical priors for the Finnish horseshoe (regularized horseshoe) (Piironen and Vehtari 2017)
  // Following the stan implementation from https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html
  alpha_hat_ij_tilde ~ normal(0,1);
  local_shrinkage_ij ~ cauchy(0,1);

  alpha_hat_eij_tilde ~ normal(0,1);
  local_shrinkage_eij ~ normal(0,1);
  tau_tilde ~ cauchy(0,1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

  // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = exp(lambda_0 + lambda_size*size_big[i] + lambda_drought*drought[i]);
    for(s in 1:S){
        alpha_eij[i,s] = exp(alpha_generic_0 + alpha_generic_size*size_big[i] + alpha_hat_ij[s] + (alpha_generic_drought + alpha_hat_eij[s]) * drought[i]);
    }
    interaction_effects[i] = sum(alpha_eij[i,] .* sp_matrix[i,]);
    
    F_hat[i] = lambda_ei[i] / (1 + interaction_effects[i]);
  }
  fecundity ~ poisson(F_hat);
}
