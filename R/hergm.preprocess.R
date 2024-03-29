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

hergm.preprocess <- function(method, max_number, initialize, network, model, hyper_prior, parametric, Clist, MHproposal, MCMCparams, maxedges, scaling, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, mean_between, eta, indicator, simulate, parallel, variational, temperature, predictions, verbose, perturb)
{

  # print("hergm.preprocess.R I: eta")
  # print(eta)

  if (is.null(verbose)) verbose <- -1 
  terms <- Clist$nterms # Number of hergm terms
  hierarchical <- vector(mode = "integer", length = terms) # Indicator: hierarchical hergm term
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
  transitiveties_ijk <- 0
  triangle_ijk <- 0
  ttriple_ijk <- 0
  if (is.null(scaling)) size.dependent <- FALSE
  else size.dependent <- TRUE
  scaling <- vector(length=terms) # Size-dependent parameterizations
  for (i in 1:terms) # For given hergm term...
    {
    if (model$terms[[i]]$name == "edges")
      {
      hierarchical[i] <- 0
      edges <- 1
      scaling[i] <- 2 # Order of subgraph count
      }
    else if (model$terms[[i]]$name == "mutual")
      {
      hierarchical[i] <- 0
      mutual <- 1
      scaling[i] <- 2 # Order of subgraph count
      }
    else if (model$terms[[i]]$name %in% c("altkstar", "balance", "ctriple", "cycle", "cyclicalties", "cyclicalweights", "dsp", "esp", "gwdegree", "gwdsp", "gwdsp", "gwesp", "gwidegree", "gwnsp", "gwodegree", "intransitive", "istar", "kstar", "localtriangle", "mstar", "mutual", "nsp", "opentriad", "ostar", "simmelian", "simmelianties", "threepath", "transitive", "transitiveties", "transitiveweights", "triadcensus", "triangle", "tripercent", "ttriple", "twopath")) # ergm term
      {
      hierarchical[i] <- 0
      dependence <- 1
      scaling[i] <- 3 # In general, the resulting rescaling is not the best possible, but the best possible rescaling is too burdensome to implement
      }
    else if (model$terms[[i]]$name == "edges_i") # hergm term
      {
      method <- "bayes"
      edges_i <- 1
      hierarchical[i] <- 1
      scaling[i] <- 2 # Order of subgraph count
      min_size_i <- 1
      if (min_size_i < min_size) min_size <- min_size_i
      if (simulate == FALSE)
        {
        degree <- vector(mode = "numeric", length = Clist$n)    
        for (i in 1:Clist$n) 
          { 
          degree[i] <- sum(network[i,])
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
      method <- "bayes"
      hierarchical[i] <- 1 
      scaling[i] <- 2 # Order of subgraph count
      min_size_i <- 1
      if (min_size_i < min_size) min_size <- min_size_i
      }     
    else if (model$terms[[i]]$name == "arcs_j") # hergm term
      {
      method <- "bayes"
      hierarchical[i] <- 1 
      scaling[i] <- 2 # Order of subgraph count
      min_size_i <- 1
      if (min_size_i < min_size) min_size <- min_size_i
      }     
    else if (model$terms[[i]]$name == "edges_ij") # hergm term
      {
      edges_ij <- 1
      hierarchical[i] <- 1 
      scaling[i] <- 2 # Order of subgraph count
      min_size_i <- 2
      if (min_size_i < min_size) min_size <- min_size_i
      }     
    else if (model$terms[[i]]$name == "mutual_i") # hergm term
      {
      mutual_i <- 1
      hierarchical[i] <- 1 
      scaling[i] <- 2 # Order of subgraph count
      dependence <- 1
      min_size_i <- 2
      if (min_size_i < min_size) min_size <- min_size_i
      }     
    else if (model$terms[[i]]$name == "mutual_ij") # hergm term
      {
      mutual_ij <- 1
      hierarchical[i] <- 1 
      scaling[i] <- 2 # Order of subgraph count
      dependence <- 1
      min_size_i <- 2
      if (min_size_i < min_size) min_size <- min_size_i
      }     
    else if (model$terms[[i]]$name == "twostar_ijk") # hergm term
      {
      twostar_ijk <- 1
      hierarchical[i] <- 1 
      scaling[i] <- 3 # Order of subgraph count
      dependence <- 1
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      }
    else if (model$terms[[i]]$name == "transitiveties_ijk") # hergm term
      {
      transitiveties_ijk <- 1
      hierarchical[i] <- 1
      scaling[i] <- 3 # Order of subgraph count
      dependence <- 1
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      }     
    else if (model$terms[[i]]$name == "triangle_ijk") # hergm term
      {
      triangle_ijk <- 1
      hierarchical[i] <- 1 
      scaling[i] <- 3 # Order of subgraph count
      dependence <- 1
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      }     
    else if (model$terms[[i]]$name == "ttriple_ijk") # hergm term
      {
      method <- "bayes"
      ttriple_ijk <- 1
      hierarchical[i] <- 1 
      scaling[i] <- 3 # Order of subgraph count
      dependence <- 1
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      }     
    else if (model$terms[[i]]$name == "ctriple_ijk") # hergm term
      {
      method <- "bayes"
      hierarchical[i] <- 1 
      scaling[i] <- 3 # Order of subgraph count
      dependence <- 1
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      }     
    else 
      {
      hierarchical[i] <- 0 # Indicator: non-hierarchical hergm term
      if (is.null(model$terms[[i]]$dependence) || (model$terms[[i]]$dependence == 1)) dependence <- 1
      }
    }
  if ((size.dependent == TRUE) && (verbose >= 0)) cat("\nSize-dependent parameterizations.\n") 
  d <- Clist$nterms # Number of parameters
  structural <- vector(mode = "integer", length = d) # Indicator: structural parameters 
  theta <- vector(mode = "numeric", length = d) 
  d1 <- 0
  d2 <- 0
  l <- 0
  decomposable <- 1 # Indicator of whether the decomposition of "inputs" into dyadic "inputs" is trivial; 
  # Note 1: the implementation is "lazy": "inputs" can be decomposed, but the current implementation is too "lazy" to try unless the decomposition is trivial
  # Note 2: if "inputs" can be decomposed into dyadic "inputs", subgraph sampling and computations can be facilitated in hergm
  for (i in 1:terms) # For given hergm term... 
    {
    number <- model$terms[[i]]$inputs[2] # Number of change statistics = number of parameters 
    if (hierarchical[i] == 0) # ergm term
      {
      d1 <- d1 + number # Increment number of non-structural parameters
      if (model$terms[[i]]$inputs[3] > 0) decomposable <- 0 # If the total number of input parameters and covariates is positive, decomposition is non-trivial; one such model term is sufficient to set decomposable to 0, meaning non-trivial
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
  if (is.null(max_number)) max_number <- Clist$n # Unspecified by user: default
  else max_number <- min(max_number, Clist$n) # Specified by user by using max_number; first option 
  if (!is.null(indicator)) indicator <- as.integer(indicator) # If non-null, make sure all elements of indicator are integers
  # ...otherwise leave indicator NULL, because NULL means non-fixed, whereas non-NULL means fixed
  if (simulate == TRUE) all_indicators_fixed <- TRUE # To simulate data, use given indicators
  else # To estimate model... 
    {
    if ((!is.null(indicator)) && (length(indicator) == Clist$n) && (sum(is.na(indicator)) == 0) && (min(indicator) >= 0) && (max(indicator) <= max_number)) all_indicators_fixed <- TRUE # ...use given indicator if indicator is not NULL and the length is correct and there are no NA
    else all_indicators_fixed <- FALSE # ...otherwise estimate all of them
    }
  if (all_indicators_fixed == TRUE) initialize <- FALSE
  if (is.null(indicator)) number_fixed <- 0
  else number_fixed <- Clist$n - sum(is.na(indicator)) # Number of fixed indicators
  if (number_fixed > 0) # Make minimum 0 
    {
    indicator <- indicator - min(indicator, na.rm=TRUE)
    }
  indicator <- as.integer(indicator) # One more time: making sure that all elements of indicator are integers
  between <- vector(mode = "numeric", length = d2) # Default: no between-block terms
  for (i in 1:terms) # For given hergm term... 
    {
    if ((hierarchical[i] == 1) && (model$terms[[i]]$name %in% c("edges_ij", "mutual_ij"))) between[i] <- 1
    else between[i] <- 0 
    }
  if (sum(between) >= 1) # The following is used when simulate==TRUE and is.null(eta)==TRUE and unused otherwise
    {
    if (is.null(mean_between)) mean_between <- 3 / (Clist$n - 1) # If a homogeneous Bernoulli model without covariates governs between-block relations and node i is in its own block and all n-1 other nodes are in other blocks, then the expected degree of node i is 3/(n-1) * (n-1) = 3; otherwise, it is less  
    }
  if ((edges + edges_i + edges_ij > 1) || (mutual + mutual_i + mutual_ij > 1)) 
    {
    if (verbose >= 0) cat("\n\n")
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
  if (simulate == TRUE) hyper_prior <- 0
  else # simulate == FALSE
    {
    if (d2 == 0) hyper_prior <- 0 # No hierarchical model terms
    }
  #if ((simulate == FALSE) && (all_indicators_fixed == TRUE)) hyper_prior <- 0
  #else if ((simulate == TRUE) && (null$indicator == FALSE)) hyper_prior <- 0
  prior_assumptions <- vector(length=2)
  if (hyper_prior == 0) prior_assumptions[1] <- 0 # If TRUE, hierarchical prior, otherwise non-hierarchical prior
  else prior_assumptions[1] <- 1 
  if (parametric == FALSE) prior_assumptions[2] <- 0 # If TRUE, parametric prior (which may or may not be hierarchical)
  else prior_assumptions[2] <- 1 
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
  length.eta <- d1 + ((max_number+1)*d2)
  if (null$eta || (length(eta) < length.eta)) eta <- vector(mode = "numeric", length = length.eta)
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
      eta_sigma11[i1,i1] <- eta_sigma[i,i]
      }
    else # hergm term
      {
      i2 <- i2 + 1
      eta_mean2[i2] <- eta_mean[i]
      eta_sigma22[i2,i2] <- eta_sigma[i,i]
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
  if (null$alpha) alpha <- 1
  if (null$indicator) 
    {
    indicator <- vector(mode = "numeric", length = Clist$n)
    }
  if (initialize == TRUE) 
    {
    indicator <- hergm.initialize(network, max_number, perturb)
    }
  indicator <- indicator - min(indicator, na.rm=TRUE) # Make sure that the block memberhips start at 0
  for (i in 1:Clist$n) 
    {
    if (is.na(indicator[i]) == FALSE) 
      {
      if (number_fixed > 0) indicator[i] <- -abs(indicator[i]) - 1 # The additional subtraction of -1 is important, because then all fixed memberships are -max_number...-1, which is what the C code expects: see h_ergm_latent.c
      else indicator[i] <- indicator[i] # Specified and non-fixed indicator of block membership
      }
    else indicator[i] <- 0 
    } 
  for (i in 1:terms) # We need to reset the input vectors of model terms; otherwise, when the user does not specify the number of blocks as an argument of the model terms, then the input vectors are intialized as vectors of length C n, where C > 0
    {
    hergm.theta <- rep.int(1, max_number+1)
    if (hierarchical[i] == 1) model$terms[[i]]$inputs <- c(0, 1, 1+length(indicator)+length(hergm.theta), c(max_number, indicator, hergm.theta))
    }
  Clist <- ergm.Cprepare(network, model)
  # Important note: the following works as long as the statistics are either monotone subgraph counts or functions of the graph such that the function takes its minimum value at the empty graph and its maximum value at the complete graph
  #print(model$terms)
  network$terms <- terms
  network$max_number <- max_number
  g0 <- matrix(0, Clist$n, Clist$n)
  g1 <- matrix(1, Clist$n, Clist$n)
  g0 <- as.network(g0, type="adjacency", directed=is.directed(network))
  g1 <- as.network(g1, type="adjacency", directed=is.directed(network))
  g0$terms <- terms
  g0$max_number <- max_number
  g1$terms <- terms
  g1$max_number <- max_number
  for (i in 1:terms)
    {
    if (hierarchical[i] == 1)
      {
      #print(model$terms[[i]]$name)
      f0 <- paste("g0", "~", model$terms[[i]]$name)
      f0 <- as.formula(f0)
      model$minval <- summary(f0)
      f1 <- paste("g1", "~", model$terms[[i]]$name)
      f1 <- as.formula(f1)
      model$maxval <- summary(f1)
      }
    }
  if ((simulate == FALSE) && (prior_assumptions[1] == 1) && (d2 > 0) && (verbose >= 0) && (method == "bayes"))
    {
    cat("\nPrior:")
    cat("\n- alpha ~ Gamma(", alpha_shape, ",", alpha_rate, ")", sep="")
    cat("\n- eta_mean ~ N(", eta_mean_mean[1], ",", (1 / eta_mean_precision[1]), ")", sep="")
    cat("\n- eta_precision ~ Gamma(", eta_precision_shape[1], ",", eta_precision_rate[1], ")", sep="")
    cat("\n")
    }
  max_iteration <- MCMCparams$samplesize
  if (simulate == TRUE) predictions <- TRUE
  number_terms <- length_mcmc(d1, d2, max_number, Clist$n, predictions)
  if (simulate == TRUE) dimension <- MCMCparams$samplesize
  else dimension <- min(MCMCparams$samplesize, 10000)
  mcmc <- vector(mode = "numeric", length = (dimension * number_terms))
  if (Clist$dir == FALSE) max_edges <- Clist$n * (Clist$n - 1) / 2 # Undirected
  else max_edges <- Clist$n * (Clist$n - 1) # Directed
  if ((simulate == TRUE) && (predictions == TRUE))
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
  if ((variational == TRUE) && (terms == 2) && (edges + edges_ij > 0) && (twostar_ijk + transitiveties_ijk + triangle_ijk + ttriple_ijk > 0)) model_type <- 1 # Some algorithms, e.g., variational algorithms, are restricted to a subset of model specifications; the specificiations it works: hergms with two model terms, one ergm term (either edges or edges_ij) and one hierarchical term (either twostar_ijk or transitiveties + triangle_ijk or ttriple_ijk)
  else model_type <- 0
  if ((simulate == FALSE) && (dependence > 0))
    {
    default <- c(1,10)
    if (is.null(temperature))
      {
      temperature <- default
      warning("minimum and maximum temperature are NULL and are replaced by ",temperature[1]," and ",temperature[2],", respectively.\n",sep="")
      }
    else if (length(temperature) < 2)
      {
      temperature <- default
      warning("either minimum or maximum temperature are unspecified and are replaced by ",temperature[1]," and ",temperature[2],", respectively.\n",sep="")
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
    if ((verbose == 1) && (simulate == FALSE) && (method == "bayes"))
      {
      cat("\nMinimum or maximum temperature: ",temperature[1]," and ",temperature[2],", respectively.\n",sep="")
      }
    }

  # print("hergm.preprocess.R II: eta")
  # print(eta)
 
  # Build object hergm_list
  hergm_list <- list()
  hergm_list$model_type <- model_type
  hergm_list$prior_assumptions <- prior_assumptions
  hergm_list$hyper_prior <- prior_assumptions[1]
  hergm_list$dependence <- dependence
  hergm_list$hierarchical <- hierarchical
  hergm_list$scaling <- scaling
  hergm_list$d <- d
  hergm_list$d1 <- d1
  hergm_list$d2 <- d2
  hergm_list$structural <- structural  
  hergm_list$min_size <- min_size
  hergm_list$max_number <- max_number
  hergm_list$number_fixed <- number_fixed
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
  hergm_list$mean_between <- mean_between
  hergm_list$predictions <- predictions
  hergm_list$verbose <- verbose
  hergm_list$MHproposal <- MHproposal
  hergm_list$maxedges <- 5*maxedges # Please note: the multiplication is motivated by the fact that otherwise, when ERGM simulates complete graphs, it will return empty graphs; therefore, by adding 1, we are tricking ERGM into believing that complete graphs are not complete graphs; ERGM does the same in ergm.san.R: if(z$status==1) maxedges <- 5*maxedges 
  hergm_list$Clist <- Clist
  hergm_list$simulate <- simulate
  hergm_list$all_indicators_fixed <- all_indicators_fixed
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

  # print("hergm.preprocess.R III: hergm_list$eta")
  # print(hergm_list$eta)

  hergm_list
}

