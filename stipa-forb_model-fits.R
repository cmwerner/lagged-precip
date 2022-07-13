# This script will run the empirical model fits for each focal species and each
#       environmental covariate. A separate script will then make the empirical
#       figures for the manuscript

#rm(list = ls())
library(rstan)
library(here)
library(HDInterval)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load in the appropriate data
#load(here("BH_simulations/test_multiple_simulations.RData")) # load in stipa_dat file

# number of interacting species -- not treating the intraspecific competition differently here
species_names <- stipa_dat %>% 
  select(-c(block, size, treatment, seed)) %>%
  names()
S <- species_names %>% length() 

# extract data needed to run the model 
N <- as.integer(nrow(stipa_dat)) # number of rows in the data
fecundity <- as.integer(stipa_dat$seed)
size_big <- ifelse(stipa_dat$size == "big", 1, 0)
drought <- ifelse(stipa_dat$treatment == "drought", 1, 0)
sp_matrix <- stipa_dat %>% 
  select(-c(block, size, treatment, seed)) %>%
  as.matrix()

# Set the parameters defining the regularized horseshoe prior
tau0 <- 1
slab_df <- 4 
slab_scale <- sqrt(2) 

# Create the data vectors to be passed to rstan for subsequent model fits
data_vec <- c("N", "S", "fecundity", "size_big", "drought", "sp_matrix", "tau0", "slab_scale", "slab_df")

prelim_fit <- stan(file = "stipa-forb_preliminary.stan", data = data_vec, iter = 3000, 
                  chains = 3, control = list(adapt_delta = 0.99, max_treedepth = 15))
prelim_posteriors <- extract(prelim_fit)

####
save(prelim_fit, prelim_posteriors, file = "stipa-forb_prelim-fit_3.rdata")
#load("stipa-forb_prelim-fit_3.rdata")
####

PrelimFit <- prelim_fit # just changing between my and Topher's notation style
PrelimPosteriors <- prelim_posteriors # just changing between my and Topher's notation style

# Examine diagnostic plots and determine if the model fit is adequate to move
#       forward with the final fit
# First examine the distribution of Rhat and effective sample size values
hist(summary(PrelimFit)$summary[,"Rhat"])
hist(summary(PrelimFit)$summary[,"n_eff"])
# Visually examine the traceplots for key parameters
traceplot(PrelimFit, pars = "lambda_0")
traceplot(PrelimFit, pars = "lambda_size")
traceplot(PrelimFit, pars = "lambda_drought")
traceplot(PrelimFit, pars = "alpha_generic_0")
traceplot(PrelimFit, pars = "alpha_generic_size")
# traceplot(PrelimFit, pars = "alpha_hat_ij") # takes a long time to draw with a lot of species
# traceplot(PrelimFit, pars = "alpha_hat_eij")
# Check for parameter correlations and divergent transitions
pairs(PrelimFit, pars = c("lambda_0","lambda_size","lambda_drought", "alpha_generic_0"))
# Check for autocorrelation in key parameters
acf(PrelimPosteriors$alpha_generic[,1])
acf(PrelimPosteriors$alpha_generic[,2])
acf(PrelimPosteriors$alpha_intra[,1])
acf(PrelimPosteriors$alpha_intra[,2])


#### If the diagnostic plots don't reveal any problems wiht the model fit, now
#       move on to determining which parameters warrant inclusion in the final
#       model (i.e. the data pulled their posteriors away from 0). The final model
#       will then be run with only these species-specific parameters, but without
#       the regularized horseshoe priors.
Inclusion_ij <- rep(0, S)
Inclusion_eij <- rep(0, S)
IntLevel <- 0.95 # previously used 0.5
for(s in 1:S){
  Ints_ij <- hdi(PrelimPosteriors$alpha_hat_ij[,s], credMass = IntLevel)
  Ints_eij <- hdi(PrelimPosteriors$alpha_hat_eij[,s], credMass = IntLevel)
  if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
    Inclusion_ij[s] <- 1
  }
  if(Ints_eij[1] > 0 | Ints_eij[2] < 0){
    Inclusion_eij[s] <- 1
  }
}
Inclusion_ij
Inclusion_eij

plot(PrelimFit, pars = "alpha_hat_ij")
plot(PrelimFit, pars = "alpha_hat_eij")

include.df <- data.frame(species_names, Inclusion_ij, Inclusion_eij)

# want to dis-include species only present in one plot in either drought or not-drought treatments
# as these single values can skew inclusion. This is a first pass, but might want to fine-tune it
sp.abundance <- stipa %>%
  group_by(species, treatment) %>%
  dplyr::summarise(
    count.sum = sum(count)
  ) %>% filter(!(species == 'none'))

sp.low.presence <- unique(sp.abundance$species[sp.abundance$count.sum < 2])
include.df <- include.df %>% mutate(
  sufficient = case_when(
    species_names %in% sp.low.presence ~ 0,
    TRUE ~ 1
  ),
  Inclusion_ij = Inclusion_ij * sufficient,
  Inclusion_eij = Inclusion_eij * sufficient
)

Inclusion_ij <- include.df$Inclusion_ij
Inclusion_eij <- include.df$Inclusion_eij

# Reset initial conditions with values from the preliminary fit
# WORK ON THIS LATER
# ChainInitials <- list(lambda_0 = colMeans(PrelimPosteriors$lambda_0), 
#                       lambda_size = colMeans(PrelimPosteriors$lambda_size), 
#                       lambda_drought = colMeans(PrelimPosteriors$lambda_drought), 
#                       alpha_generic_tilde = colMeans(PrelimPosteriors$alpha_generic_tilde), 
#                       alpha_hat_ij_tilde = colMeans(PrelimPosteriors$alpha_hat_ij_tilde), 
#           #            local_shrinkage_ij = colMeans(PrelimPosteriors$local_shrinkage_ij), 
#           #           c2_tilde = mean(PrelimPosteriors$c2_tilde), tau_tilde = mean(PrelimPosteriors$tau_tilde), 
#           #            local_shrinkage_eij = colMeans(PrelimPosteriors$local_shrinkage_eij)
#                       alpha_hat_eij_tilde = colMeans(PrelimPosteriors$alpha_hat_eij_tilde))
# InitVals <- list(ChainInitials, ChainInitials, ChainInitials)

# Run the final fit of the model
FinalDataVec <- c("N", "S", "fecundity", "size_big", "drought", "sp_matrix",
             "Inclusion_ij", "Inclusion_eij")
FinalFit <- stan(file = "stipa-forb_final.stan", data = FinalDataVec, iter = 3000,
                 chains = 3, control = list(adapt_delta = 0.95))
Posteriors <- extract(FinalFit)

# If the fit looks good, safe the final output here
save(FinalFit, Posteriors, Inclusion_ij, Inclusion_eij, file = "stipa-forb_final-fit_3.rdata")
#load("stipa-forb_final-fit.rdata")

# Examine the same diagnostic plots as before to check for problems with the 
#       final model fit
pairs(FinalFit, pars = c("lambda_0","lambda_size","lambda_drought", "alpha_generic_0"))
hist(summary(FinalFit)$summary[,"Rhat"])
hist(summary(FinalFit)$summary[,"n_eff"])
traceplot(FinalFit, pars = "lambda_0")
traceplot(FinalFit, pars = "lambda_size")
traceplot(FinalFit, pars = "lambda_drought")
traceplot(FinalFit, pars = "alpha_generic_0")
which(Inclusion_ij == 1)
#traceplot(FinalFit, pars = "alpha_hat_ij")
which(Inclusion_eij == 1)
#traceplot(FinalFit, pars = "alpha_hat_eij")

# Double check the autocorrelation
acf(Posteriors$lambda_0)
acf(Posteriors$lambda_drought)
acf(Posteriors$alpha_generic[,1])
acf(Posteriors$alpha_generic[,2])
acf(Posteriors$alpha_intra[,1])
acf(Posteriors$alpha_intra[,2])
for(s in 1:S){
  if(Inclusion_ij[s] == 1){
    quartz()
    acf(Posteriors$alpha_hat_ij[,s])
  }
  if(Inclusion_eij[s] == 1){
    quartz()
    acf(Posteriors$alpha_hat_eij[,s])
  }
}



