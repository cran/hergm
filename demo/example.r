.libPaths("~/r/runs") # Specify library path
setwd("~/r/runs") # Specify work directory
library(hergm) # Load library: see help(package = "hergm") 
data(example) # Load data: see ?example
k <- 3 # Number of blocks
simulate <- FALSE # If TRUE, simulation of networks, otherwise Bayesian inference
parallel <- 1 # Parallel computing: number of computing nodes; if more than one, computing is parallel
samplesize <- 1e+4 # If simulate = TRUE, number of networks to be sampled, otherwise number of draws from posterior; if computing is parallel, number of draws per computing node
total_samplesize <- min(parallel * 1200, parallel * samplesize)
name <- "example"
mcmc <- hergm( # See ?hergm
     d ~ edges_i(k), # Formula: network ~ terms: see ?network, ?hergm.terms, ?edges_i
     alpha = 1, # Clustering parameter of truncated Dirichlet process / stick-breaking prior
     alpha_shape = 1, # Shape parameter of Gamma prior of clustering parameter
     alpha_rate = 1, # Rate (inverse scale) parameter of Gamma prior of clustering parameter 
     eta = c(-1), # Natural parameters of exponential-family model
     eta_mean = c(-1),  # Means of Gaussian baseline distribution of truncated Dirichlet process / stick-breaking prior
     eta_sd = c(2), # Standard deviations of the Gaussian baseline distribution of truncated Dirichlet process / stick-breaking prior
     parallel = parallel, # Parallel computing: number of computing nodes; if more than one, computing is parallel
     simulate = simulate, # If TRUE, simulation of networks, otherwise Bayesian inference
     seeds = c(-2010, 2010), # Seed of pseudo-random number generator; if computing is parallel, number of seeds must equal number of computing nodes
     samplesize = samplesize, # If simulate = TRUE, number of networks to be sampled, otherwise number of draws from posterior; if computing is parallel, number of draws per computing node
     burnin = 1e+4, # If simulate = TRUE, number of burn-in iterations 
     interval = 1e+2, # If simulate = TRUE, number of proposals between sampled networks
     output = TRUE, # Output: if TRUE, full output, otherwise limited output
     name = name, # If output = TRUE, name of project is used to name output files
     verbose = 1 # Console output: -1 none, 0 short, +1 long
     )
if (simulate == FALSE) processed_mcmc <- hergm.postprocess( # See ?hergm.postprocess
     n = 15, # Number of nodes
     k = k, # Number of blocks
     d1 = 0, # Number of ergm terms
     d2 = 1, # Number of hergm terms
     burnin = 0, # Number of burn-in iterations
     samplesize = total_samplesize, # Total sample size, including number of burn-in iterations and, if parallel > 1, summed across computing nodes
     mcmc = mcmc, # MCMC sample in the form of vector
     output = TRUE, # If TRUE, full output, including relabeled MCMC sample; otherwise limited output
     name = name, # If output = TRUE, name of project is used to name output files
     )

