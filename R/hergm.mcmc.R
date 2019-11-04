###########################################################################
# Copyright 2009 Michael Schweinberger                                    #
#                                                                         #
# This file is part of hergm.                                             #
#                                                                         # 
#    hergm is free software: you can redistribute it and/or modify        #
#    it under the terms of the GNU General Public License as published by #
#    the Free Software Foundation, either version 3 of the License, or    #
#    (at your option) any later version.                                  #
#                                                                         # 
#    hergm is distributed in the hope that it will be useful,             #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
#    GNU General Public License for more details.                         #
#                                                                         #
#    You should have received a copy of the GNU General Public License    #
#    along with hergm.  If not, see <http://www.gnu.org/licenses/>.       #
#                                                                         # 
###########################################################################

hergm.mcmc <- function(parameterization, method, sample_size_multiplier_blocks, formula, max_number, initialize, initialization_method, network, model, hyper_prior, parametric, MHproposal, MCMCparams, verbose, scaling, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, mean_between, eta, indicator, parallel, simulate, seeds, mh_scale, variational, temperature, predictions, perturb, estimate_parameters, n_em_step_max, max_iter, initial_estimate, NR_step_len, NR_step_len_multiplier, NR_max_iter) 
{

  # print("hergm.mcmc.R: eta")
  # print(eta)

  # Prepare
  Clist <- ergm.Cprepare(network, model)
  if (Clist$dir == FALSE) maxedges <- Clist$n * (Clist$n - 1) / 2 # Undirected
  else maxedges <- Clist$n * (Clist$n - 1) # Directed
  if ((simulate == FALSE) && (method == "bayes")) scalefactor <- hergm.set.mcmc(method, max_number, initialize, network, model, hyper_prior, parametric, MHproposal, MCMCparams, verbose, indicator, scaling, alpha_shape, alpha_rate, alpha,  eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, mean_between, eta, simulate, parallel, seeds, predictions, variational, temperature, mh_scale, perturb) # The last argument is the initial value of the scale factor
  if ((simulate == TRUE) || (method == "bayes")) hergm_list <- hergm.preprocess(method, max_number, initialize, network, model, hyper_prior, parametric, Clist, MHproposal, MCMCparams, maxedges, scaling, alpha_shape, alpha_rate, alpha,  eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, mean_between, eta, indicator, simulate, parallel, variational, temperature, predictions, verbose, perturb)
  if ((simulate == FALSE) && (method == "bayes")) 
    {
    hergm_list$scalefactor <- scalefactor # Set scale factor
    if (hergm_list$hyper_prior == 1)
      {
      k <- number_clusters(hergm_list$alpha, hergm_list$Clist$n)
      if (verbose >= 0) 
        { 
          cat("\nDirichlet process:")
          cat("\n- scaling parameter: ", hergm_list$alpha, sep = "")
          cat("\n- mean of number of non-empty blocks: ", 
              formatC(k$mean, digits = 2, width = 4, format = "f", mode = "real"), sep = "") 
          cat("\n- variance of number of non-empty blocks: ", 
              formatC(k$variance, digits = 2, width = 4, format = "f", mode = "real"), sep = "") 
          cat("\n")
        }   
      }
    }
  # Run 
  if (simulate == TRUE) parallel <- 1
  if ((simulate == TRUE) || (method == "bayes"))
    {
    if (parallel > 1) # Parallel computing
      {
      number <- parallel # Specify number of computing nodes on cluster
      if (length(seeds) < number) 
        {
        maximum <- as.integer(.Machine$integer.max / 10)
        cluster.seeds <- sample(0:maximum, size=number, replace=FALSE)
        }
      else 
        {
        cluster.seeds <- seeds
        }
      cluster <- makePSOCKcluster(number)
      clusterEvalQ(cluster, library(hergm))
      s <- clusterApplyLB(cluster, cluster.seeds[1:number], hergm.wrapper, hergm_list)
      stopCluster(cluster)
      mcmc <- append(s[[1]]$mcmc, s[[2]]$mcmc)
      if (number > 2) 
        {
        for (i in 3:number) mcmc <- append(mcmc, s[[i]]$mcmc)
        }
      sample <- list()
      sample$mcmc <- mcmc
      }
    else sample <- hergm.wrapper(seeds[1], hergm_list) # Non-parallel computing
    }
  else
    {
    if (is.null(max_number) == TRUE)
      {
      if (is.null(indicator) == TRUE) max_number <- 1
      else max_number <- length(unique(indicator))
    }
    if (max_number == 1) {all_indicators_fixed = TRUE; indicator = rep(1,network$gal$n)}
    if ((!is.null(indicator)) && (length(indicator) == Clist$n) && (sum(is.na(indicator)) == 0) && (min(indicator) >= 0) && (max(indicator) <= max_number)) all_indicators_fixed <- TRUE # ...use given indicator if indicator is not NULL and the length is correct and there are no NA
    else all_indicators_fixed <- FALSE # ...otherwise estimate all of them
    if (all_indicators_fixed == FALSE) indicator <- NULL
    sample <- hergm.large(network = network, formula = formula, parameterization = parameterization, 
                          indicator = indicator, same_between_blocks = TRUE, max_number = max_number,
                          number_cores = parallel, initialization_method = initialization_method, estimate_parameters = estimate_parameters, 
                          verbose = verbose, n_em_step_max = n_em_step_max, max_iter = max_iter, 
                          initial_estimate = initial_estimate, MCMCparams = MCMCparams, sample_size_multiplier_blocks = sample_size_multiplier_blocks, seeds = seeds, NR_step_len = NR_step_len, NR_step_len_multiplier = NR_step_len_multiplier, NR_max_iter = NR_max_iter)
    if (sample$estimation_status == "failed") 
      { 
      estimate_parameters <- FALSE
      }
    }
  # Store
  object <- list() 
  object$method <- method
  object$n <- Clist$n
  object$network <- network
  object$formula <- formula 
  object$parameterization <- parameterization
  object$estimate_parameters <- estimate_parameters
  object$max_number <- max_number
  if (method == "bayes") 
    {
    object$number_fixed <- hergm_list$number_fixed
    object$d1 <- hergm_list$d1
    object$d2 <- hergm_list$d2
    object$hyper_prior <- hergm_list$hyper_prior
    object$parallel <- hergm_list$parallel
    object$simulate <- hergm_list$simulate
    object$sample_size <- min(10000, hergm_list$MCMCparams$samplesize)
    }
  else 
    {
    object$simulate <- simulate
    object$relabel <- 0
    }
  if (simulate == TRUE) 
    {
    number_edges <- sample$sample_heads[1]
    if (number_edges == 0)
      {
      sample_heads <- NULL
      sample_tails <- NULL
      }
    else
      { 
      sample$sample_heads <- sample$sample_heads[-1]
      sample$sample_tails <- sample$sample_tails[-1]
      object$heads <- sample$sample_heads[1:number_edges]
      object$tails <- sample$sample_tails[1:number_edges]
      object$edgelist <- cbind(object$heads, object$tails)
      }
    }
  object$predictions <- predictions
  if (method == "bayes") 
    {
    object$sample <- sample$mcmc
    }
  else 
    {
    object$results <- sample
    object$mlergm_out <- sample$mlergm_out 
    }
  object$verbose <- verbose
  object
}

