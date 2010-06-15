hergm.set.mcmc <- function(nw, model, MHproposal, MCMCparams, verbose, name, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, parallel, simulate, seeds, output, scalefactor)
{

  test_verbose <- -1
  verbose <- test_verbose

  # Prepare I
  if (verbose >= 0) cat("\nMetropolis-Hastings algorithm: scale factor and acceptance rate:")
  cp_samplesize <- MCMCparams$samplesize # Store
  MCMCparams$samplesize <- round(parallel * cp_samplesize / 100)
  if (MCMCparams$samplesize <= 10) MCMCparams$samplesize <- 10
  cp_parallel <- parallel # Store
  parallel <- 1

  # Prepare II
  Clist <- ergm.Cprepare(nw, model)
  maxedges <- max(50000, Clist$nedges)
  hergm_list <- hergm.preprocess(nw, model, Clist, MHproposal, MCMCparams, maxedges, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, simulate = FALSE, parallel, output = FALSE, name = "", verbose = test_verbose)

  # Metropolis-Hastings: finding scale factor
  if ((is.null(scalefactor)) || (scalefactor <= 0.0)) scalefactor <- 1.0
  hergm_list$scalefactor <- scalefactor
  if (verbose >= 0) cat("\n")
  s <- hergm.wrapper(seeds[1], hergm_list)
  iteration <- 1
  if (verbose >= 0) cat("\n(", iteration, ")", " ", 
                      formatC(scalefactor, digits = 4, width = 6, format = "f", mode = "real"), 
                      " ",
                      formatC(s$mh_accept, digits = 4, width = 6, format = "f", mode = "real"), 
                      sep = "")
  while ((s$mh_accept < 0.25) && (iteration <= 10))
    {  
    iteration <- iteration + 1
    hergm_list <- hergm.preprocess(nw, model, Clist, MHproposal, MCMCparams, maxedges, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, simulate = FALSE, parallel, output = FALSE, name = "", verbose = test_verbose)
    scalefactor <- scalefactor / 2
    hergm_list$scalefactor <- scalefactor
    if (verbose >= 0) cat("\n")
    s <- hergm.wrapper(seeds[1], hergm_list)
    if (verbose >= 0) cat("\n(", iteration, ")", " ", 
                        formatC(scalefactor, digits = 4, width = 6, format = "f", mode = "real"), 
                        " ",
                        formatC(s$mh_accept, digits = 4, width = 6, format = "f", mode = "real"), 
                          sep = "")
    }

  parallel <- cp_parallel # Reset
  MCMCparams$samplesize <- cp_samplesize # Reset

  scalefactor

}

