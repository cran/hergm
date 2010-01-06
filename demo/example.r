library(hergm) # Load library: see help(package = "hergm") 
data(example) # Load data: see ?example
n <- 15 # Number of nodes
k <- 3 # Number of blocks
simulate <- FALSE # If TRUE, simulation of networks, otherwise Bayesian inference
parallel <- 1 # Parallel computing: number of computing nodes; if more than one, computing is parallel
samplesize <- 1.2e+5 # If simulate = TRUE, number of networks to be sampled, otherwise number of draws from posterior; if computing is parallel, number of draws per computing node
returned_samplesize <- min(6e+3, samplesize) # If simulate = FALSE, number of recorded draws from posterior; if compting is parallel, number of recorded draws per computing node
name <- "example"
mcmc <- hergm( # See ?hergm
     d ~ edges_i(k), # Formula: network ~ terms: see ?network, ?hergm.terms, ?edges_i
     alpha_shape = 1.0, # Shape parameter of Gamma prior of clustering parameter
     alpha_rate = 1.0, # Rate (inverse scale) parameter of Gamma prior of clustering parameter 
     eta_mean_mean = -2.0,  # Means of Gaussian prior of means of Gaussian baseline distribution of truncated Dirichlet process / stick-breaking prior
     eta_mean_sd = 1e+1, # Standard deviations of the Gaussian baseline distribution of truncated Dirichlet process / stick-breaking prior
     eta_precision_shape = 1.0, # Shape parameter of Gamma prior of precision parameter of baseline distribution of truncated Dirichlet process / stick-breaking prior 
     eta_precision_rate = 1e-1, # Rate (inverse scale) parameter of Gamma prior of precision parameter of baseline distribution of truncated Dirichlet process / stick-breaking prior  
     parallel = parallel, # Parallel computing: number of computing nodes; if more than one, computing is parallel
     simulate = simulate, # If TRUE, simulation of networks, otherwise Bayesian inference
     seeds = c(0.0842877348260602, # Seed of pseudo-random number generator; if computing is parallel, number of seeds must equal number of computing nodes
               0.5533580354178491, 
               0.7950986713125660, 
               0.1252978907576672, 
               0.1832790205073627, 
               0.8640834893401168, 
               0.2271214423797233, 
               0.6870248310503468, 
               0.9876543708220379, 
               0.2423582594866208),   
     samplesize = samplesize, # If simulate = TRUE, number of networks to be sampled, otherwise number of draws from posterior; if computing is parallel, number of draws per computing node
     burnin = 1e+4, # If simulate = TRUE, number of burn-in iterations 
     interval = 1e+2, # If simulate = TRUE, number of proposals between sampled networks
     mh_scale = 1.0, # If simulate = FALSE, scale factof of candidate-generating distribution of Metropolis-Hastings algorithm
     output = TRUE, # Output: if TRUE, full output, otherwise limited output
     name = name, # If output = TRUE, name of project is used to name output files
     verbose = 1 # Console output: -1 none, 0 short, +1 long
     )
if (simulate == FALSE) processed_mcmc <- hergm.postprocess( # See ?hergm.postprocess
     n = n, # Number of nodes
     k = k, # Number of blocks
     d1 = 0, # Number of ergm terms
     d2 = 1, # Number of hergm terms
     parallel = parallel, # Parallel computing: number of computing nodes used to generate MCMC sample
     burnin = 1e+3, # Number of burn-in iterations; if parallel > 1, 
     samplesize = returned_samplesize, # MCMC sample size, including number of burn-in iterations; if parallel > 1, MCMC sample size per computing node
     mcmc = mcmc, # MCMC sample in the form of vector
     relabel = TRUE, # Relabel MCMC sample
     output = TRUE, # If TRUE, full output, including relabeled MCMC sample; otherwise limited output
     name = name # If output = TRUE, name of project is used to name output files
     )

