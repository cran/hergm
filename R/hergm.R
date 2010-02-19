hergm <- function(formula, 
                 alpha = NULL,
                 alpha_shape = NULL, 
                 alpha_rate = NULL, 
                 eta = NULL,
                 eta_mean = NULL, 
                 eta_sd = NULL,
                 parallel = 1, 
                 simulate = FALSE, 
                 seeds = NULL, 
                 samplesize = 1e+5, 
                 burnin = 1e+4, 
                 interval = 1e+2,
                 output = FALSE,
                 verbose = -1, 
                 name = NULL,
                 ...) 
{
  options(warn = -1)
  control <- control.ergm()
  options()
  nw <- ergm.getnetwork(formula)
  model <- ergm.getmodel(formula, nw, drop=control$drop, expanded=TRUE)
  MCMCsamplesize <- samplesize
  Clist <- ergm.Cprepare(nw, model)
  Clist.miss <- ergm.design(nw, model, verbose=verbose)
  d <- Clist$nstats
  constraints <- ~.
  MHproposal <- MHproposal(constraints,weights=control$prop.weights,control$prop.args,nw,model,class="c")
  MHproposal.miss <- MHproposal("randomtoggleNonObserved",control$prop.args, nw, model)
  MCMCparams=c(control,list(samplesize=MCMCsamplesize,burnin=burnin,interval=interval,maxit=1,Clist.miss=Clist.miss,mcmc.precision=control$mcmc.precision))
  MCMCparams$stats <- matrix(0,ncol=Clist$nstats,nrow=MCMCparams$samplesize)
  MCMCparams$meanstats <- Clist$meanstats
  print(
    system.time(
      sample <- hergm.mcmc(nw, model, MHproposal, MCMCparams, verbose, name, alpha_shape, alpha_rate, alpha, eta_mean, eta_sd, eta, parallel, simulate, seeds, output)
    )
  )
  cat("\n")
  sample
}


