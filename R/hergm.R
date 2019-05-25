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

hergm <- function(formula, 
                  max_number = 2,
                  hierarchical = TRUE,
                  parametric = FALSE,
                  parameterization = "offset",
                  initialize = FALSE,
                  initialization_method = 1,
                  estimate_parameters = TRUE,
                  initial_estimate = NULL,
                  n_em_step_max = 100,
                  max_iter = 4,
                  perturb = FALSE,
                  scaling = NULL,
                  alpha = NULL,
                  alpha_shape = NULL, 
                  alpha_rate = NULL, 
                  eta = NULL,
                  eta_mean = NULL, 
                  eta_sd = NULL,
                  eta_mean_mean = NULL,
                  eta_mean_sd = NULL,
                  eta_precision_shape = NULL,
                  eta_precision_rate = NULL,
                  mean_between = NULL,
                  indicator = NULL,
                  parallel = 1, 
                  simulate = FALSE, 
                  method = "ml",
                  seeds = NULL, 
                  sample_size = NULL, 
                  sample_size_multiplier_blocks = 20,
                  NR_max_iter = 200,
                  NR_step_len = 1,
                  NR_step_len_multiplier = 0.2, 
                  interval = 1024,
                  burnin = 16*interval, 
                  mh.scale = 0.25,
                  variational = FALSE,
                  temperature = c(1,100),
                  predictions = FALSE,
                  posterior.burnin = 2000,
                  posterior.thinning = 1,
                  relabel = 1,
                  number_runs = 1,
                  verbose = 0, 
                  ...) 
{
  original.formula <- formula
  oopts <- options(warn = -1)
  control <- control.ergm()
  options()
  network <- hergm.getnetwork(formula, max_number)
  # if (sum(network[,] == 1) == 0) stop("\nNetwork is extreme: terminating...\n\n") # Simplistic check
  control$drop <- FALSE
  model <- ergm_model(formula, network, drop=control$drop, expanded=TRUE)
  Clist <- ergm.Cprepare(network, model)
  for (i in 1:Clist$nterms) 
    {
    if (model$terms[[i]]$name %in% c("edges_i", "arcs_i", "arcs_j", "ctriple_ijk", "ttriple_ijk")) method <- "bayes"
    }
  if ((simulate == FALSE) && (is.null(indicator) == FALSE)) method <- "bayes" 
  if (is.null(sample_size)) 
    {
    if (method == "ml") sample_size <- 2500
    else sample_size = 1e+5 
    }
  if ((method == "bayes") && (simulate == FALSE) && (sample_size < 100)) sample_size <- 100
  MCMCsamplesize <- sample_size
  ## Commenting old line and putting changed line underneath 
  ## ergm.design changed in ergm-master / changing for compatibility 
  #Clist.miss <- ergm.design(network, model, verbose=FALSE)
  Clist.miss <- ergm.design(network, verbose=FALSE)
  constraints <- ~.
  MHproposal <- ergm_proposal(constraints, weights=control$MCMC.prop.weights, control$MCMC.prop.args, network, class="c")
  MCMCparams <- c(control, list(samplesize=MCMCsamplesize, burnin=burnin, interval=interval,
                                maxit=1, Clist.miss=Clist.miss, mcmc.precision=control$MCMLE.mcmc.precision))
  if (((min(10000, sample_size) - posterior.burnin) / posterior.thinning) < 1000) 
    {
    posterior.burnin <- 0
    posterior.thinning <- 1
    }
  s <- min(MCMCparams$samplesize, 10000)
  MCMCparams$stats <- matrix(0,ncol=Clist$nstats,nrow=s)
  MCMCparams$target.stats <- Clist$target.stats
  object <- hergm.mcmc(parameterization, method, sample_size_multiplier_blocks, original.formula, max_number, initialize, initialization_method, network, model, hyper_prior=hierarchical, parametric, MHproposal, MCMCparams, verbose, scaling, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, mean_between, eta, indicator, parallel, simulate, seeds, mh.scale, variational, temperature, predictions, perturb, estimate_parameters, n_em_step_max, max_iter, initial_estimate = initial_estimate, NR_step_len = NR_step_len, NR_step_len_multiplier = NR_step_len_multiplier, NR_max_iter = NR_max_iter)
  if ((max_number >= 10) && (relabel == 1)) relabel <- 2 
  if ((simulate == FALSE) && (method == "bayes")) 
    {
    object <- hergm.postprocess(object=object, burnin=posterior.burnin, thinning=posterior.thinning, relabel=relabel, number_runs=number_runs)
    object$mcmc.diagnostics <- mcmc.diagnostics.hergm(object)
    }
  on.exit(options(oopts))

  return(structure(object, class="hergm"))
}

