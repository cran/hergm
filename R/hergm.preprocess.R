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

hergm.preprocess <- function(nw, model, Clist, MHproposal, MCMCparams, maxedges, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, indicator, simulate, parallel, temperature, output, name, verbose) # Michael
{
  if (is.null(verbose)) verbose <- -1 
  terms <- Clist$nterms # Number of hergm terms
  hierarchical <- vector(mode = "integer", length = terms) # Indicator: hierarchical hergm term
  max_number <- Clist$n # Default: (maximum) number of categories
  min_size <- Clist$n # Structural parameters corresponding to categories with min_size..n nodes show up in hergm pmf
  dependence <- 0 # Default: no dyad-dependence
  q_i <- vector(mode = "numeric", length = Clist$n) # Proposal distribution of nodes: sample nodes, then sample category indicator of sampled node
  for (i in 1:Clist$n) 
    {
    q_i[i] <- 1 / Clist$n # Default: discrete uniform; depending on the model specification, default values may be overwritten
    }
  edges <- 0
  edges_i <- 0
  edges_ij <- 0
  mutual <- 0
  mutual_i <- 0
  mutual_ij <- 0
  twostar_ijk <- 0
  triangle_ijk <- 0
  ttriple_ijk <- 0
  for (i in 1:terms) # For given hergm term...
    {
    if (model$terms[[i]]$name == "edges")
      {
      hierarchical[i] <- 0
      edges <- 1
      }
    else if (model$terms[[i]]$name == "mutual")
      {
      mutual <- 1
      hierarchical[i] <- 0
      dependence <- 1 # See ergm package: if dyad-independence term, non-null and FALSE, otherwise null
      }
    else if (model$terms[[i]]$name %in% c("altkstar", "balance", "ctriple", "cycle", "cyclicalties", "cyclicalweights", "dsp", "esp", "gwdegree", "gwdsp", "gwdsp", "gwesp", "gwidegree", "gwnsp", "gwodegree", "intransitive", "istar", "kstar", "localtriangle", "mstar", "nsp", "opentriad", "ostar", "simmelian", "simmelianties", "threepath", "transitive", "transitiveties", "transitiveweights", "triadcensus", "triangle", "tripercent", "ttriple", "twopath")) # ergm term
      {
      hierarchical[i] <- 0
      dependence <- 1
      }
    else if (model$terms[[i]]$name == "edges_i") # hergm term
      {
      edges_i <- 1
      hierarchical[i] <- 1
      min_size_i <- 1
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i   
      if (simulate == FALSE)
        {
        degree <- vector(mode = "numeric", length = Clist$n)    
        for (i in 1:Clist$n) 
          { 
          degree[i] <- sum(nw[i,])
          if (verbose >= 2) cat("\ndegree[",i,"] = ",degree[i])
          }
        range <- max(degree) - min(degree) 
        max_odds <- 6
        T <- range / log(max_odds) 
        sum <- 0
        for (i in 1:Clist$n) 
          { 
          q_i[i] <- exp(degree[i] / T) # Implication: the log-odds of the probability of selecting nodes with maximum degree cannot exceed 2, or the odds cannot exceed exp(range / T) = exp(log(max_odds)) = max_odds; note that T can be interpreted as the tempature (inverse canonical parameter of discrete exponential family distribution 
          sum <- sum + q_i[i]
          }
        if (verbose >= 2) cat("\nProbability of selection of nodes:")
        for (i in 1:Clist$n) 
          { 
          q_i[i] <- q_i[i] / sum
          if (verbose >= 2) cat(" ",q_i[i]) 
          }
        if (verbose >= 2) cat("\n")
        }     
      }
    else if (model$terms[[i]]$name == "arcs_i") # hergm term
      {
      hierarchical[i] <- 1 
      min_size_i <- 1
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "arcs_j") # hergm term
      {
      hierarchical[i] <- 1 
      min_size_i <- 1
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "edges_ij") # hergm term
      {
      edges_ij <- 1
      hierarchical[i] <- 1 
      min_size_i <- 2
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "mutual_i") # hergm term
      {
      mutual_i <- 1
      hierarchical[i] <- 1 
      dependence <- 1
      min_size_i <- 2
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "mutual_ij") # hergm term
      {
      mutual_ij <- 1
      hierarchical[i] <- 1 
      dependence <- 1
      min_size_i <- 2
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "twostar_ijk") # hergm term
      {
      twostar_ijk <- 1
      hierarchical[i] <- 1 
      dependence <- 1
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "triangle_ijk") # hergm term
      {
      triangle_ijk <- 1
      hierarchical[i] <- 1 
      dependence <- 1
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "ttriple_ijk") # hergm term
      {
      ttriple_ijk <- 1
      hierarchical[i] <- 1 
      dependence <- 1
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "ctriple_ijk") # hergm term
      {
      hierarchical[i] <- 1 
      dependence <- 1
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else 
      {
      hierarchical[i] <- 0 # Indicator: non-hierarchical hergm term
      if (is.null(model$terms[[i]]$dependence) || (model$terms[[i]]$dependence == 1)) dependence <- 1
      }
    }
  decomposable <- 1 # Indicator of whether the decomposition of "inputs" into dyadic "inputs" is trivial; 
  # Note 1: the implementation is "lazy": "inputs" can be decomposed, but the current implementation is too "lazy" to try unless the decomposition is trivial
  # Note 2: if "inputs" can be decomposed into dyadic "inputs", subgraph sampling and computations can be facilitated in hergm
  for (i in 1:terms)
    {
    if (hierarchical[i] == 0) # non-hierarchical term
      {
      if (model$terms[[i]]$inputs[3] > 0) decomposable <- 0 # If the total number of input parameters and covariates is positive, decomposition is non-trivial; one such model term is sufficient to set decomposable to 0, meaning non-trivial
      }
    }
  d <- Clist$nstats # Number of parameters
  structural <- vector(mode = "integer", length = d) # Indicator: structural parameters 
  theta <- vector(mode = "numeric", length = d) 
  d1 <- 0
  d2 <- 0
  l <- 0
  for (i in 1:terms) # For given hergm term... 
    {
    number <- model$terms[[i]]$inputs[2] # Number of change statistics = number of parameters 
    if (hierarchical[i] == 0) # ergm term
      {
      d1 <- d1 + number # Increment number of non-structural parameters
      for (k in 1:number) # Set indicator: structural parameter
        {
        l <- l + 1
        structural[l] <- 0 # Non-structural parameter
        }
      if (is.null(eta)) theta[i] <- 0
      else theta[i] <- eta[i]
      }     
    else 
      {
      d2 <- d2 + number # Increment number of structural parameters
      for (k in 1:number) # Set indicator: structural parameter
        {
        l <- l + 1
        structural[l] <- 1 # Structural parameter
        }
      theta[i] <- 1 
      }
    }
  if (d2 == 0) max_number <- 1
  between <- vector(mode = "numeric", length = d2) # Default: no between-block terms
  for (i in 1:terms) # For given hergm term... 
    {
    if (hierarchical[i] == 1)
      {
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i != max_number) 
        {
        cat("\n\n")
        error_message <- paste("The number of blocks should either be left unspecified or the same number of blocks must be specified.")      
        stop(error_message, call. = FALSE)
        } 
      if (model$terms[[i]]$name == "edges_ij") between[i] = 1
      else if (model$terms[[i]]$name == "mutual_ij") between[i] = 1
      }     
    }
  if (edges + edges_i + edges_ij > 1) nonidentifiable <- FALSE
  else if (mutual + mutual_i + mutual_ij > 1) nonidentifiable <- FALSE
  else identifiable <- TRUE
  if (identifiable == FALSE) 
    {
    cat("\n\n")
    error_message <- paste("The model is non-identifiable: check model terms and drop some of them.")      
    stop(error_message, call. = FALSE)
    } 
  # Prior
  null <- list()
  null$alpha_shape <- is.null(alpha_shape)
  null$alpha_rate <- is.null(alpha_rate)
  null$alpha <- is.null(alpha)
  null$eta_mean_mean <- is.null(eta_mean_mean)
  null$eta_mean_precision <- is.null(eta_mean_sd)
  null$eta_precision_shape <- is.null(eta_precision_shape)
  null$eta_precision_rate <- is.null(eta_precision_rate)
  null$eta_mean <- is.null(eta_mean)
  null$eta_precision <- is.null(eta_sd)
  null$eta <- is.null(eta)
  null$indicator <- is.null(indicator)
  hyper_prior <- ((null$alpha_shape==FALSE) && (null$alpha_rate==FALSE) && (null$eta_mean_mean==FALSE) && (null$eta_mean_precision==FALSE) && (null$eta_precision_shape==FALSE) && (null$eta_precision_rate==FALSE))
  if (hyper_prior == TRUE) cat("\nHyper-prior is fully specified and will be used.\n")
  else cat("\nHyper-prior is not fully specified and will not be used.\n")
  if (null$eta) eta <- matrix(data = 0, nrow = d2, ncol = max_number)
  if (null$alpha_shape) alpha_shape <- 1
  if (null$alpha_rate) alpha_rate <- 1
  if (null$eta_mean_mean) eta_mean_mean <- vector(mode = "numeric", length = d2)
  eta_mean_precision <- vector(mode = "numeric", length = d2)
  if (null$eta_mean_precision) 
    {
    for (i in 1:d2) eta_mean_precision[i] <- 1 
    }
  else
    {
    for (i in 1:d2) eta_mean_precision[i] <- 1 / (eta_mean_sd[i] * eta_mean_sd[i])
    }
  if (null$eta_precision_shape) eta_precision_shape <- 10
  if (null$eta_precision_rate) eta_precision_rate <- 10
  if (null$eta_mean) eta_mean <- vector(mode = "numeric", length = d)
  eta_sigma <- matrix(data = 0, nrow = d, ncol = d)
  for (i in 1:d) 
    {
    if (null$eta_precision) eta_sigma[i,i] <- 1 
    else eta_sigma[i,i] <- eta_sd[i] * eta_sd[i]
    }
  if (null$indicator) indicator <- rep.int(0, Clist$n)
  # Marginal Gaussian priors:
  eta_mean1 <- vector(mode = "numeric", length = d1)
  eta_mean2 <- vector(mode = "numeric", length = d2) 
  eta_sigma11 <- matrix(data = 0, nrow = d1, ncol = d1)
  eta_sigma12 <- matrix(data = 0, nrow = d1, ncol = d2) 
  eta_sigma21 <- matrix(data = 0, nrow = d2, ncol = d1) 
  eta_sigma22 <- matrix(data = 0, nrow = d2, ncol = d2)
  i1 <- 0
  i2 <- 0
  for (i in 1:d) 
    {
    if (structural[i] == 0) # ergm term
      {
      i1 <- i1 + 1
      eta_mean1[i1] <- eta_mean[i]     
      j1 <- 0
      j2 <- 0
      for (j in 1:d)
        { 
        if (structural[j] == 0) # ergm term
          {
          j1 <- j1 + 1
          eta_sigma11[i1,j1] <- eta_sigma[i,j]
          }
        else # Hierarchical, non-hierarchical hergm term
          {
          j2 <- j2 + 1
          eta_sigma12[i1,j2] <- eta_sigma[i,j]
          }
        }
      }
    else # hergm term
      {
      i2 <- i2 + 1
      eta_mean2[i2] <- eta_mean[i]
      j1 <- 0
      j2 <- 0
      for (j in 1:d)
        {
        if (structural[j] == 0) # Hierarchical, non-hierarchical hergm term
          {
          j1 <- j1 + 1
          eta_sigma21[i2,j1] <- eta_sigma[i,j]
          }
        else # Hierarchical, hierarchical hergm term
          {
          j2 <- j2 + 1
          eta_sigma22[i2,j2] <- eta_sigma[i,j]
          }        
        }
      }
    }
  # Marginal Gaussian prior of non-structural parameters...
  if (d1 == 0) 
    {
    cf1 <- 1
    p1 <- 1
    }
  else 
    {
    cf1 <- t(chol(eta_sigma11)) # ...Cholesky factor of covariance matrix satisfying eta_sigma11 = cf1 * t(cf1)
    p1 <- solve(eta_sigma11) # ...precision (inverse covariance) matrix
    }
  #print(cf1 %*% t(cf1))
  #print(p1 %*% eta_sigma11)
  # Conditional Gaussian prior of structural parameters given non-structural parameters...
  if (d2 == 0)
    {
    b <- 1
    cf2 <- 1
    p2 <- 1
    }
  else 
    {
    b <- eta_sigma21 %*% p1
    eta_sigma2 <- eta_sigma22 - (b %*% eta_sigma12) # ...covariance matrix
    cf2 <- t(chol(eta_sigma2)) #...Cholesky factor of covariance matrix satisfying eta_sigma2 = cf2 * t(cf2)
    p2  <- solve(eta_sigma2) # ...precision (inverse covariance) matrix
    }
  #print(cf2 %*% t(cf2))
  #print(p2 %*% eta_sigma2)
  if (null$eta) eta <- matrix(data = 0, nrow = d2, ncol = max_number)
  if (null$alpha) alpha <- 1
  if (null$eta) eta <- vector(mode = "numeric", length = d)
  if (null$indicator) indicator <- vector(mode = "numeric", length = Clist$n)
  if (min(indicator) < 0) indicator <- indicator - min(indicator) # Shift to 0
  if (max(indicator) > (max_number - 1)) indicator <- indicator - (max(indicator) - (max_number - 1)) # Shift to max_number-1
  max_iteration <- MCMCparams$samplesize
  number_terms <- length_mcmc(d1, d2, max_number, Clist$n)
  if (simulate == TRUE) dimension <- MCMCparams$samplesize
  else dimension <- min(MCMCparams$samplesize, 10000)
  mcmc <- vector(mode = "numeric", length = (dimension * number_terms))
  if (Clist$dir == FALSE) max_edges <- Clist$n * (Clist$n - 1) / 2 # Undirected
  else max_edges <- Clist$n * (Clist$n - 1) # Directed
  if ((simulate == TRUE) && (output == TRUE))
    {
    sample_heads <- vector(mode = "numeric", length = (max_iteration * (max_edges + 1))) # max_edges + 1: the number of edges is added on every iteration
    sample_tails <- vector(mode = "numeric", length = (max_iteration * (max_edges + 1))) # max_edges + 1: the number of edges is added on every iteration
    }
  else
    {
    sample_heads <- 0
    sample_tails <- 0
    }
  mh_accept <- vector(mode = "numeric", length = 2 + (Clist$n - min_size))
  if ((terms == 2) && (edges + edges_ij > 0) && (twostar_ijk + triangle_ijk + ttriple_ijk > 0)) model_type <- 1 # Some algorithms, e.g., variational algorithms, are restricted to a subset of model specifications; the specificiations it works: hergms with two model terms, one ergm term (either edges or edges_ij) and one hierarchical term (either twostar_ijk or triangle_ijk or ttriple_ijk)
  else model_type <- 0
  if ((simulate == FALSE) && (dependence > 0))
    {
    default <- c(0,10)
    if (is.null(temperature))
      {
      temperature <- default
      cat("\nWarning: minimum and maximum temperature are NULL and are replaced by ",temperature[1]," and ",temperature[2],", respectively.\n",sep="")
      }
    else if (length(temperature) < 2)
      {
      temperature <- default
      cat("\nWarning: either minimum or maximum temperature are unspecified and are replaced by ",temperature[1]," and ",temperature[2],", respectively.\n",sep="")
      }
    temperature[1] <- abs(temperature[1])
    temperature[2] <- abs(temperature[2])
    if (temperature[1] > temperature[2])
      {
      min <- min(temperature)
      max <- max(temperature)
      temperature[1] <- min
      temperature[2] <- max
      }
    cat("\nMinimum or maximum temperature: ",temperature[1]," and ",temperature[2],", respectively.\n",sep="")
    }
 
  # Build object hergm_list
  hergm_list <- list()
  hergm_list$model_type <- model_type
  hergm_list$hyper_prior <- hyper_prior
  hergm_list$dependence <- dependence
  hergm_list$hierarchical <- hierarchical
  hergm_list$decomposable <- decomposable
  hergm_list$d <- d
  hergm_list$d1 <- d1
  hergm_list$d2 <- d2
  hergm_list$structural <- structural  
  hergm_list$min_size <- min_size
  hergm_list$max_number <- max_number
  hergm_list$null <- null
  hergm_list$alpha_shape <- alpha_shape
  hergm_list$alpha_rate <- alpha_rate
  hergm_list$alpha <- alpha
  hergm_list$eta_mean_mean <- eta_mean_mean
  hergm_list$eta_mean_precision <- eta_mean_precision
  hergm_list$eta_precision_shape <- eta_precision_shape
  hergm_list$eta_precision_rate <- eta_precision_rate
  hergm_list$eta_mean1 <- eta_mean1
  hergm_list$eta_mean2 <- eta_mean2
  hergm_list$theta <- theta
  hergm_list$b <- as.vector(b)
  hergm_list$cf1 <- as.vector(cf1)
  hergm_list$cf2 <- as.vector(cf2)
  hergm_list$p1 <- as.vector(p1)
  hergm_list$p2 <- as.vector(p2)
  hergm_list$eta <- eta
  hergm_list$indicator <- indicator
  hergm_list$max_iteration <- max_iteration
  hergm_list$terms <- terms
  hergm_list$between <- between
  hergm_list$output <- output
  hergm_list$name <- name
  hergm_list$verbose <- verbose
  hergm_list$MHproposal <- MHproposal
  hergm_list$maxedges <- maxedges
  hergm_list$Clist <- Clist
  hergm_list$simulate <- simulate
  hergm_list$mh_accept <- mh_accept
  hergm_list$q_i <- q_i
  hergm_list$parallel <- parallel
  hergm_list$temperature <- temperature
  hergm_list$model <- model

  #print(hergm_list) # If you want to print hergm_list, print here, before the long vectors are added to hergm_list

  hergm_list$MCMCparams <- MCMCparams
  hergm_list$sample_heads <- as.vector(sample_heads)
  hergm_list$sample_tails <- as.vector(sample_tails)
  hergm_list$mcmc <- as.vector(mcmc)

  hergm_list
}

