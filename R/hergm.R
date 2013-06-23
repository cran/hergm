###########################################################################
# Copyright 2009 Nobody                                                   #
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
                 indicator = NULL,
                 parallel = 1, 
                 simulate = FALSE, 
                 seeds = NULL, 
                 samplesize = 1e+5, 
                 burnin = 1e+4, 
                 interval = 1e+2,
                 mh_scale = NULL,
                 temperature = c(1,10),
                 output = TRUE,
                 verbose = -1, 
                 name = NULL,
                 ...) 
{
  options(warn = -1)
  control <- control.ergm()
  options()
  nw <- ergm.getnetwork(formula)
  control$drop <- FALSE
  model <- ergm.getmodel(formula, nw, drop=control$drop, expanded=TRUE)
  MCMCsamplesize <- samplesize
  Clist <- ergm.Cprepare(nw, model)
  Clist.miss <- ergm.design(nw, model, verbose=verbose)
  if (verbose >= 2)
    {
    for (i in 1:Clist$nterms)
      {
      if (model$terms[[i]]$name == "edges_i") 
        { 
        cat("\n\n")
        print(summary(nw ~ degree(0:(Clist$n-1)), drop = TRUE))
        }
      else if (model$terms[[i]]$name == "arcs_i") 
        {
        cat("\n\n")
        print(summary(nw ~ odegree(0:(Clist$n-1)), drop = TRUE))
        }
      else if (model$terms[[i]]$name == "arcs_j") 
        {
        cat("\n\n")
        print(summary(nw ~ idegree(0:(Clist$n-1)), drop = TRUE))
        }
      }
    }
  d <- Clist$nstats
  constraints <- ~.
  MHproposal <- MHproposal(constraints, weights=control$MCMC.prop.weights, control$MCMC.prop.args, nw, class="c")
# MHproposal.miss <- MHproposal("randomtoggleNonObserved",arguments=control$MCMC.prop.args, nw=nw,reference=~Bernoulli)
  MCMCparams=c(control,list(samplesize=MCMCsamplesize,burnin=burnin,interval=interval,maxit=1,Clist.miss=Clist.miss,mcmc.precision=control$MCMLE.mcmc.precision))
  MCMCparams$stats <- matrix(0,ncol=Clist$nstats,nrow=MCMCparams$samplesize)
  MCMCparams$target.stats <- Clist$target.stats
  print(
    system.time(
      sample <- hergm.mcmc(nw, model, MHproposal, MCMCparams, verbose, name, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, indicator, parallel, simulate, seeds, mh_scale, temperature, output)
    )
  )
  cat("\n")
  sample
}


