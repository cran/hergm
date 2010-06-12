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

hergm.set.mcmc <- function(nw, model, MHproposal, MCMCparams, verbose, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, parallel, seeds, output, scalefactor)
{

  # Prepare I
  if (verbose >= 0) cat("\nMetropolis-Hastings algorithm: scale factor and acceptance rate:")
  cp_samplesize <- MCMCparams$samplesize # Store
  MCMCparams$samplesize <- round(parallel * cp_samplesize / 100)
  if (MCMCparams$samplesize < 100) MCMCparams$samplesize <- 100
  else if (MCMCparams$samplesize > 10000) MCMCparams$samplesize <- 10000

  # Prepare II
  Clist <- ergm.Cprepare(nw, model)
  maxedges <- max(50000, Clist$nedges)
  hergm_list <- hergm.preprocess(nw, model, Clist, MHproposal, MCMCparams, maxedges, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, indicator = NULL, simulate = FALSE, parallel = 1, output = FALSE, name = "", verbose = -1)

  # Metropolis-Hastings: finding scale factor
  min_accept <- 0.30
  hergm_list$scalefactor <- scalefactor
  s <- hergm.wrapper(seeds[1], hergm_list)
  iteration <- 1
  if (verbose >= 0) cat("\n(", iteration, ")", " ", 
                      formatC(scalefactor[1], digits = 4, width = 6, format = "f", mode = "real"), 
                      " ",
                      formatC(scalefactor[2], digits = 4, width = 6, format = "f", mode = "real"), 
                      " ",
                      formatC(s$mh_accept[1], digits = 4, width = 6, format = "f", mode = "real"), 
                      " ",
                      formatC(s$mh_accept[2], digits = 4, width = 6, format = "f", mode = "real"), 
                      sep = "")
  while ((min(s$mh_accept) < min_accept) && (iteration <= 10))
    {  
    iteration <- iteration + 1
    hergm_list <- hergm.preprocess(nw, model, Clist, MHproposal, MCMCparams, maxedges, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, indicator = NULL, simulate = FALSE, parallel = 1, output = FALSE, name = "", verbose = -1)
    if (s$mh_accept[1] < min_accept) scalefactor[1] <- scalefactor[1] / 2
    if (s$mh_accept[2] < min_accept) scalefactor[2] <- scalefactor[2] / 2
    hergm_list$scalefactor <- scalefactor
    s <- hergm.wrapper(seeds[1], hergm_list)
    if (verbose >= 0) cat("\n(", iteration, ")", " ", 
                          formatC(scalefactor[1], digits = 4, width = 6, format = "f", mode = "real"), 
                          " ",
                          formatC(scalefactor[2], digits = 4, width = 6, format = "f", mode = "real"), 
                          " ",
                          formatC(s$mh_accept[1], digits = 4, width = 6, format = "f", mode = "real"), 
                          " ",
                          formatC(s$mh_accept[2], digits = 4, width = 6, format = "f", mode = "real"), 
                          sep = "")
    }

  MCMCparams$samplesize <- cp_samplesize # Reset

  scalefactor

}

