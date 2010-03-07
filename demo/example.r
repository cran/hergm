library(hergm) # Load library: see help(package = "hergm") 
data(example) # Load data: see ?example
n <- 15 # Number of nodes
k <- 3 # Number of blocks
simulate <- FALSE # If TRUE, simulation of networks, otherwise Bayesian inference
parallel <- 1 # Parallel computing: number of computing nodes; if more than one, computing is parallel
samplesize <- 1.2e+5 # If simulate = TRUE, number of networks to be sampled, otherwise number of draws from posterior; if computing is parallel, number of draws per computing node
returned_samplesize <- min(1.2e+4, samplesize) # If simulate = FALSE, number of recorded draws from posterior; if compting is parallel, number of recorded draws per computing node
name <- "example"
mcmc <- hergm( # See ?hergm
     d ~ edges_i(k), # Formula: network ~ terms: see ?network, ?hergm.terms, ?edges_i
     alpha_shape = 1.0, # Hyperprior: shape parameter of Gamma prior of scaling parameter
     alpha_rate = 1.0, # Hyperprior: rate (inverse scale) parameter of Gamma prior of scaling parameter
     eta_mean_mean = 0.0, # Hyperprior: mean of mean of Gaussian base distribution
     eta_mean_sd = 10.0, # Hyperprior: standard deviation of mean of Gaussian base distribution
     eta_precision_shape = 1.0, # Hyperprior: shape of Gamma prior of precision of Gaussian base distribution
     eta_precision_rate = 0.1, # Hyperprior: rate (inverse scale) of Gamma prior of precision of Gaussian base distribution
     eta_mean = c(0.0), # Prior: means of Gaussian base distribution and of Gaussian prior of covariate parameters
     eta_sd = c(10.0), # Prior: standard deviation of Gaussian base distribution and of Gaussian prior of covariate parameters
     parallel = parallel, # Parallel computing: number of computing nodes; if more than one, computing is parallel
     simulate = simulate, # If TRUE, simulation of networks, otherwise Bayesian inference
     seeds = c(3395,1664482,4095377,1473557,9170474,6344540,3736735,692447,3227407,3237189), # Seed of pseudo-random number generator; if computing is parallel, number of seeds must equal number of processors
     samplesize = samplesize, # If simulate = TRUE, number of networks to be sampled, otherwise number of draws from posterior; if computing is parallel, number of draws per computing node
     burnin = 1e+4, # If simulate = TRUE, number of burn-in iterations 
     interval = 1e+2, # If simulate = TRUE, number of proposals between sampled networks
     mh_scale = 1.0, # If simulate = FALSE, scale factor of candidate-generating distribution of Metropolis-Hastings algorithm
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
     burnin = trunc(returned_samplesize / 10), # Number of burn-in iterations; if parallel > 1, 
     samplesize = returned_samplesize, # MCMC sample size, including number of burn-in iterations; if parallel > 1, MCMC sample size per computing node
     mcmc = mcmc, # MCMC sample in the form of vector
     relabel = TRUE, # Relabel MCMC sample
     output = TRUE, # If TRUE, full output
     name = name # If output = TRUE, name of project is used to name output files
     )

