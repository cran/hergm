hergm.preprocess <- function(nw, model, Clist, MHproposal, MCMCparams, maxedges, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, indicator = NULL, simulate, parallel, output, name, verbose) # Michael
{
  if (is.null(verbose)) verbose <- -1 
  terms <- Clist$nterms # Number of hergm terms	
  hierarchical <- vector(mode = "integer", length = terms) # Indicator: hierarchical hergm term
  max_number <- Clist$n # Default: (maximum) number of categories
  min_size <- Clist$n # Structural parameters corresponding to categories with min_size..n nodes show up in hergm pmf
  dependence <- 0 # Default: no dyad-dependence
  n_between <- 0 # Default: no between-block terms
  q_i <- vector(mode = "numeric", length = Clist$n) # Proposal distribution of nodes: sample nodes, then sample category indicator of sampled node
  for (i in 1:Clist$n) 
    {
    q_i[i] <- 1 / Clist$n # Default: discrete uniform; depending on the model specification, default values may be overwritten
    }
  for (i in 1:terms) # For given hergm term... 
    {
    if (model$terms[[i]]$name != "mutual")
      {
      if (is.null(model$terms[[i]]$dependence)) dependence <- 1 # See ergm package: if dyad-independence term, non-null and FALSE, otherwise null
      else if (model$terms[[i]]$dependence == TRUE) dependence <- 1 # Dyad-dependence
      }
    if (model$terms[[i]]$name == "edges_i") # hergm term
      {
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
      hierarchical[i] <- 1 
      min_size_i <- 2
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      n_between <- 1
      }     
    else if (model$terms[[i]]$name == "mutual_i") # hergm term
      {
      hierarchical[i] <- 1 
      min_size_i <- 2
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "mutual_ij") # hergm term
      {
      hierarchical[i] <- 1 
      min_size_i <- 2
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "triangle_ijk") # hergm term
      {
      hierarchical[i] <- 1 
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "ttriple_ijk") # hergm term
      {
      hierarchical[i] <- 1 
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "ctriple_ijk") # hergm term
      {
      hierarchical[i] <- 1 
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "twostar_i") # hergm term
      {
      hierarchical[i] <- 1 
      min_size_i <- 1
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else if (model$terms[[i]]$name == "twostar_ijk") # hergm term
      {
      hierarchical[i] <- 1 
      min_size_i <- 3
      if (min_size_i < min_size) min_size <- min_size_i
      max_number_i <- model$terms[[i]]$inputs[4] # (Maximum) number of categories: 1st model$terms[[i]]$inputs, thus 4th element of inputs; see InitErgm.R
      if (max_number_i < max_number) max_number <- max_number_i
      }     
    else hierarchical[i] <- 0 # Indicator: non-hierarchical hergm term
    }
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
  # Prior
  if (simulate == TRUE) 
    {
    if ((is.null(alpha)) || (is.null(eta)) || (is.null(indicator))) hyper_prior <- TRUE # Posterior predictive check
    else hyper_prior <- FALSE # Simulation
    }
  else hyper_prior <- TRUE # Dirichlet process prior with hyper prior
  if (is.null(eta)) eta <- matrix(data = 0, nrow = d2, ncol = max_number)
  if (is.null(alpha_shape)) alpha_shape <- 1
  if (is.null(alpha_rate)) alpha_rate <- 1
  if (is.null(eta_mean_mean)) eta_mean_mean <- vector(mode = "numeric", length = d2)
  eta_mean_precision <- vector(mode = "numeric", length = d2)
  if (is.null(eta_mean_sd)) 
    {
    for (i in 1:d2) eta_mean_precision[i] <- 0.25
    }
  else
    {
    for (i in 1:d2) eta_mean_precision[i] <- 1 / (eta_mean_sd[i] * eta_mean_sd[i])
    }
  if (is.null(eta_precision_shape)) eta_precision_shape <- 1
  if (is.null(eta_precision_rate)) eta_precision_rate <- 1
  if (is.null(eta_mean)) eta_mean <- vector(mode = "numeric", length = d)
  eta_sigma <- matrix(data = 0, nrow = d, ncol = d)
  for (i in 1:d) 
    {
    if (is.null(eta_sd)) eta_sigma[i,i] <- 4
    else eta_sigma[i,i] <- eta_sd[i] * eta_sd[i]
    }
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
  if (is.null(eta)) eta <- matrix(data = 0, nrow = d2, ncol = max_number)
  if (is.null(alpha)) alpha <- 1
  if (is.null(eta)) eta <- vector(mode = "numeric", length = d)
  if (is.null(indicator)) indicator <- vector(mode = "numeric", length = Clist$n)

  max_iteration <- MCMCparams$samplesize
  terms <- length_mcmc(d1, d2, max_number, Clist$n)
  if (simulate == TRUE) dimension <- MCMCparams$samplesize
  else dimension <- min(MCMCparams$samplesize, 12000)
  mcmc <- vector(mode = "numeric", length = (dimension * terms))

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
  mh_accept <- 0
  call_RNGstate <- 1

  # Build object hergm_list
  hergm_list <- list()
  hergm_list$hyper_prior <- hyper_prior
  hergm_list$dependence <- dependence
  hergm_list$hierarchical <- hierarchical
  hergm_list$d <- d
  hergm_list$d1 <- d1
  hergm_list$d2 <- d2
  hergm_list$structural <- structural  
  hergm_list$min_size <- min_size
  hergm_list$max_number <- max_number
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
  hergm_list$n_between <- n_between
  hergm_list$mcmc <- as.vector(mcmc)
  hergm_list$output <- output
  hergm_list$sample_heads <- as.vector(sample_heads)
  hergm_list$sample_tails <- as.vector(sample_tails)
  hergm_list$name <- name
  hergm_list$verbose <- verbose
  hergm_list$MCMCparams <- MCMCparams
  hergm_list$MHproposal <- MHproposal
  hergm_list$maxedges <- maxedges
  hergm_list$Clist <- Clist
  hergm_list$simulate <- simulate
  hergm_list$mh_accept <- mh_accept
  hergm_list$q_i <- q_i
  hergm_list$call_RNGstate <- call_RNGstate
  hergm_list$parallel <- parallel

  hergm_list
}

