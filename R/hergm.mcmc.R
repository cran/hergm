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

hergm.mcmc <- function(nw, model, MHproposal, MCMCparams, verbose, name, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, indicator, parallel, simulate, seeds, mh_scale, output) 
{

  # Prepare
  if (simulate == FALSE) scalefactor <- hergm.set.mcmc(nw, model, MHproposal, MCMCparams, verbose, alpha_shape, alpha_rate, alpha,  eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, parallel, seeds, output, mh_scale) # The last argument is the initial value of the scale factor
  Clist <- ergm.Cprepare(nw, model)
  if (Clist$dir == FALSE) maxedges <- Clist$n * (Clist$n - 1) / 2 # Undirected
  else maxedges <- Clist$n * (Clist$n - 1) # Directed
  hergm_list <- hergm.preprocess(nw, model, Clist, MHproposal, MCMCparams, maxedges, alpha_shape, alpha_rate, alpha,  eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, indicator, simulate, parallel, output, name, verbose)
  if (simulate == FALSE) hergm_list$scalefactor <- scalefactor # Set scale factor
  else if (hergm_list$hyper_prior == 1)
    {
    k <- number_clusters(hergm_list$alpha, hergm_list$Clist$n)
    cat("\nDirichlet process:")
    cat("\n- scaling parameter: ", hergm_list$alpha, sep = "")
    cat("\n- mean of number of non-empty blocks: ", formatC(k$mean, digits = 2, width = 4, format = "f", mode = "real"), sep = "") 
    cat("\n- variance of number of non-empty blocks: ", formatC(k$variance, digits = 2, width = 4, format = "f", mode = "real"), sep = "") 
    cat("\n")   
    }
  flush.console() # Windows

  # Run
  if (simulate == TRUE) parallel <- 1
  if (parallel > 1) # Parallel computing
    {
    number <- parallel # Specify number of computer nodes on cluster
    cluster <- makeSOCKcluster(number)
    clusterEvalQ(cluster, library(hergm))
    s <- clusterApplyLB(cluster, seeds[1:number], hergm.wrapper, hergm_list)
    stopCluster(cluster)
    sample <- append(s[[1]]$mcmc, s[[2]]$mcmc)
    if (number > 2) 
      {
      for (i in 3:number) sample <- append(sample, s[[i]]$mcmc)
      }
    }
  else sample <- hergm.wrapper(seeds[1], hergm_list) # Non-parallel computing

  # Postprocess
  output_list <- list() 
  output_list$n <- Clist$n
  output_list$max_number <- hergm_list$max_number
  output_list$d1 <- hergm_list$d1
  output_list$d2 <- hergm_list$d2
  output_list$parallel <- hergm_list$parallel
  output_list$sample_size <- min(12000, hergm_list$MCMCparams$samplesize)
  if (simulate == TRUE) 
    {
    output_list$sample <- sample$sample
    number_edges <- sample$sample_heads[1]
    sample$sample_heads <- sample$sample_heads[-1]
    sample$sample_tails <- sample$sample_tails[-1]
    output_list$heads <- sample$sample_heads[1:number_edges]
    output_list$tails <- sample$sample_tails[1:number_edges]
    }
  output_list$sample <- sample$mcmc 

  output_list
}

