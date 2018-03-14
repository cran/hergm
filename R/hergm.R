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
                  max_number = NULL,
                  hierarchical = TRUE,
                  parametric = FALSE,
                  initialize = FALSE,
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
                  seeds = NULL, 
                  sample_size = 1e+5, 
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
  if (sum(network[,] == 1) == 0) stop("\nNetwork is extreme: terminating...\n\n") # Simplistic check
  control$drop <- FALSE
  model <- ergm.getmodel(formula, network, drop=control$drop, expanded=TRUE)
  MCMCsamplesize <- sample_size
  Clist <- ergm.Cprepare(network, model)
  Clist.miss <- ergm.design(network, model, verbose=verbose)
  constraints <- ~.
  MHproposal <- MHproposal(constraints, weights=control$MCMC.prop.weights, control$MCMC.prop.args, network, class="c")
  MCMCparams=c(control,list(samplesize=MCMCsamplesize,burnin=burnin,interval=interval,maxit=1,Clist.miss=Clist.miss,mcmc.precision=control$MCMLE.mcmc.precision))
  if (((min(10000, sample_size) - posterior.burnin) / posterior.thinning) < 1000) 
    {
    posterior.burnin <- 0
    posterior.thinning <- 1
    }
  s <- min(MCMCparams$samplesize, 10000)
  MCMCparams$stats <- matrix(0,ncol=Clist$nstats,nrow=s)
  MCMCparams$target.stats <- Clist$target.stats
  object <- hergm.mcmc(original.formula, max_number, initialize, network, model, hyper_prior=hierarchical, parametric, MHproposal, MCMCparams, verbose, scaling, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, mean_between, eta, indicator, parallel, simulate, seeds, mh.scale, variational, temperature, predictions, perturb)
  if (is.null(max_number)) max_number <- 1
  if ((max_number >= 10) && (relabel == 1)) relabel <- 2 
  object.hergm <- hergm.postprocess(object=object, burnin=posterior.burnin, thinning=posterior.thinning, relabel=relabel, number_runs=number_runs)
  if (simulate == FALSE) object.hergm$mcmc.diagnostics <- mcmc.diagnostics.hergm(object.hergm)
  on.exit(options(oopts))

  return(structure(object.hergm, class="hergm"))
}

