find_phi_coefficient <- function (z_memb, true_memb) {
  n <- length(z_memb)
  if (n < 2000) {
    true_matrix <- matrix(0, n, n)
    est_matrix <- true_matrix
    for (i in unique(z_memb)) {
      est_matrix[z_memb == i, z_memb == i] <- 2
    }
    for (i in unique(true_memb)) {
      true_matrix[true_memb == i, true_memb == i] <- 1
    }
    diag(true_matrix) <- diag(est_matrix) <- 0
    res_matrix <- true_matrix - est_matrix
    n_00 <- as.numeric(sum(res_matrix == 0))
    n_01 <- as.numeric(sum(res_matrix == 1))
    n_10 <- as.numeric(sum(res_matrix == -2))
    n_11 <- as.numeric(sum(res_matrix == -1))
  } else {
    n_00 <- n_01 <- n_10 <- n_11 <- 0
    n_blocks <- (n%/%1000)
    if (n%%1000 == 0) {
      sizes <- c(rep(1000, n_blocks))
    } else {
      sizes <- c(rep(1000, n_blocks - 1), n - 1000 * n_blocks)
    }
    cum_sizes <- c(0, cumsum(sizes))
    for (i in 1:n_blocks) {
      for (j in i:n_blocks) {
        temp_matrix_full_true <- matrix(0, sizes[i], 
                                        sizes[j])
        temp_matrix_full_est <- matrix(0, sizes[i], sizes[j])
        for (k in unique(z_memb)) {
          temp_matrix_full_est[z_memb[(cum_sizes[i] + 
                                         1):(cum_sizes[i + 1])] == k, z_memb[(cum_sizes[j] + 
                                                                                1):(cum_sizes[j + 1])] == k] <- 2
        }
        for (k in unique(true_memb)) {
          temp_matrix_full_true[true_memb[(cum_sizes[i] + 
                                             1):(cum_sizes[i + 1])] == k, true_memb[(cum_sizes[j] + 
                                                                                       1):(cum_sizes[j + 1])] == k] <- 1
        }
        diag(temp_matrix_full_true) <- diag(temp_matrix_full_est) <- 0
        res_matrix <- temp_matrix_full_true - temp_matrix_full_est
        n_temp_01 <- as.numeric(sum(res_matrix == 1))
        n_temp_10 <- as.numeric(sum(res_matrix == -2))
        n_temp_11 <- as.numeric(sum(res_matrix == -1))
        n_temp_00 <- sizes[i] * sizes[j] - (n_temp_01 + 
                                              n_temp_10 + n_temp_11)
        n_00 <- n_00 + n_temp_00 * ((i != j) + 1)
        n_01 <- n_01 + n_temp_01 * ((i != j) + 1)
        n_10 <- n_10 + n_temp_10 * ((i != j) + 1)
        n_11 <- n_11 + n_temp_11 * ((i != j) + 1)
      }
    }
  }
  n_00 <- n_00 - n
  n_1_dot <- n_11 + n_10
  n_0_dot <- n_01 + n_00
  n_dot_1 <- n_11 + n_01
  n_dot_0 <- n_10 + n_00
  first_product <- (n_11/sqrt(n_1_dot)/sqrt(n_0_dot)) * (n_00/sqrt(n_dot_1)/sqrt(n_dot_0))
  second_product <- (n_10/sqrt(n_1_dot)/sqrt(n_0_dot)) * (n_01/sqrt(n_dot_1)/sqrt(n_dot_0))
  phi_coef <- first_product - second_product
  phi_coef
}


plot_neighborhoods <- function(network, partition) {
  plot(
    as.network(network),
    vertex.col = partition,
    edge.col = same_clusters(partition),
    main = ""
  )
}

permute_tau <- function(tau, labels) {
  new_tau <- tau
  for (i in 1:length(labels)) {
    #temp[z_memb == i] <- labels[as.numeric(names(labels))==i]
    new_tau[, i] <- tau[, labels[i]]
  }
  new_tau
}

logit <- function(x)
  qlogis(x)
inv.logit <- function(x)
  exp(x) / (1 + exp(x))

#Spectral clustering
spec_clust <- function(network, max_number) {
  network <- as.matrix(network)
  n <- nrow(network)
  n_vec = ceiling(sqrt(n))
  b <- eigen(network, symmetric = TRUE)$vectors[, 1:n_vec]
  c <- kmeans(
    b,
    centers = max_number,
    nstart = 100,
    iter.max = 20,
    algorithm = "Hartigan-Wong"
  )
  c.ind <- c$cluster
  z_memb <- factor(c.ind, levels = 1:max_number)
  z_memb
}

#Finish this function, need to get rid of empty clusters.
check_clusters <-
  function(z_memb, network, max_number, min_size = 3) {
    n <- nrow(network)
    n_vec = ceiling(sqrt(n))
    b <- eigen(network, symmetric = TRUE)$vectors[, 1:n_vec]
    z_memb <- factor(z_memb, levels = 1:max_number)
    repeat {
      z_memb_temp <- as.vector(z_memb)
      z_memb <- factor(z_memb_temp, levels = 1:max_number)
      bad_clusters <- which(table(z_memb) < min_size)
      if (length(bad_clusters) == 0)
        break

      bad_cluster <- which(table(z_memb) < 2)[1]
      donor <- which.max(table(z_memb))
      indeces <- c(bad_clusters, donor)
      temp <-
        kmeans(
          b[z_memb %in% indeces,],
          centers = 2,
          nstart = 100,
          iter.max = 20,
          algorithm = "Hartigan-Wong"
        )
      for (i in 1:length(indeces)) {
        z_memb_temp[z_memb %in% indeces][temp$cluster == i] <- indeces[i]
      }
      z_memb <- factor(z_memb_temp, levels = 1:max_number)
    }
    z_memb
  }

#Simulate hergm
simulate_hergm <- function(formula, coef_w, coef_b, z_memb, parameterization = "standard",
                           directed = TRUE, matrix_sbm = NULL)
{
  #Initialization
  directed <- is.directed(ergm.getnetwork(formula))
  n_nodes <- length(z_memb)
  block_sizes <- table(z_memb)
  memb_inds <- as.numeric(names(block_sizes))

  #Parse formula
  rhs <- paste(deparse(formula[[3]]), collapse = "")
  lhs <- paste(deparse(formula[[2]]), collapse = "")

  #mode = 0 simulating from scratch
  #mode = 1 if formula contains a network
  mode = 0
  if (exists(lhs)) {
    temp <- get0(lhs)
    if (is.network(temp)){
      mode <- 1
    }
  }
  if (parameterization == "offset") { 
    model <- ergm_model(formula, ergm.getnetwork(formula))
    etamap <- model$etamap
    param_names <- model$coef.names
    if (is.curved(model$formula)) {
      num_curved <- length(etamap$curved)
      for (ii in 1:num_curved) {
        which_curved <- etamap$curved[[ii]]$to
        param_names <- param_names[-which_curved[-c(1, 2)]]
      }
    }
    edge_loc <- which(param_names == "edges")
  }
  #we assume coef_b always has log n parameterization // adjusting this for other parameterizations 
  if (parameterization == "size") { 
    coef_b <- coef_b * log(n_nodes)
  } else if (parameterization == "offset") { 
    coef_b <- coef_b - log(n_nodes)
  }
  if (is.null(matrix_sbm)) {
    network <- matrix(rbinom(n_nodes^2, 1, inv.logit(coef_b)), n_nodes, n_nodes)
  } else {
    network <- apply(matrix_sbm, c(1,2), function(x) rbinom(1,1,inv.logit(x)))
  }
  if (directed == FALSE) {
    network[lower.tri(network)] <- t(network)[lower.tri(network)]
  }

  if (mode == 0) {
    #if (parameterization == "offset") { 
    #  model <- ergm_model(formula, ergm.getnetwork(formula))
    #  etamap <- model$etamap
    #  param_names <- model$coef.names
    #  if (is.curved(model$formula)) {
    #    num_curved <- length(etamap$curved)
    #    for (ii in 1:num_curved) {
    #      which_curved <- etamap$curved[[ii]]$to
    #      param_names <- param_names[-which_curved[-c(1, 2)]]
    #    }
    #  }
    #  edge_loc <- which(param_names == "edges")
    #}
  #Fill in within neighborhood parts
  for (i in memb_inds) {
    n_nodes_in_block = sum(z_memb == i)
    if (directed) { 
      burnin_ <- max(choose(n_nodes_in_block, 2) * 2 * 4, 1e+4)
      interval_ <- max(choose(n_nodes_in_block, 2) * 2, 1000) 
    } else { 
      burnin_ <- max(choose(n_nodes_in_block, 2) * 4, 1e+4)
      interval_ <- max(choose(n_nodes_in_block, 2), 1000) 
    }
    if (directed == TRUE) {
      lhs <- paste("network(", n_nodes_in_block, ")", sep = "")
    } else {
      lhs <- paste("network(", n_nodes_in_block, ", directed = FALSE)", sep = "")
    }
    form <- as.formula(paste(lhs, "~", rhs))
    if (parameterization == "size") {
      temp_network = simulate(form, coef = coef_w*log(n_nodes_in_block),
                              control = control.simulate(MCMC.burnin = burnin_, MCMC.interval = interval_))
    } else if (parameterization == "offset") {
      coef_w_ <- coef_w 
      coef_w_[edge_loc] <- coef_w_[edge_loc] - log(n_nodes_in_block)
      temp_network = simulate(form, coef = coef_w_,
                              control = control.simulate(MCMC.burnin = burnin_, MCMC.interval = interval_))
    } else {
      temp_network = simulate(form, coef = coef_w_,
                              control = control.simulate(MCMC.burnin = burnin_, MCMC.interval = interval_))
    }
    network[z_memb == i, z_memb == i] <- as.matrix(temp_network)
  }
  } else {
    for (i in memb_inds) {
      n_nodes_in_block = sum(z_memb == i)
      v_id <- which(z_memb == i)
      subgraph_temp <- get.inducedSubgraph(temp, v = v_id)
      lhs_sub <- "subgraph_temp"
      form <- as.formula(paste(lhs_sub, "~", rhs))
      if (parameterization == "size") {
        temp_network = simulate(form, coef = coef_w*log(n_nodes_in_block))
      } else if (parameterization == "ofset") {
        coef_w_ <- coef_w 
        coef_w_[edge_loc] <- coef_w_[edge_loc] - log(n_nodes_in_block)
        temp_network <- simulate(form, coef = coef_w_)
      } else {
        temp_network = simulate(form, coef = coef_w)
      }
      network[z_memb == i, z_memb == i] <- as.matrix(temp_network)
    }
  }
  as.network(network, directed = directed)
}


par_get_subgraph <- function(group, network, z_memb, memb_inds, formula, 
                             edge_cov_name, edge_cov_list) {
  net_list <- rep(list(NULL), length(group))
  nodes_in_group <- which(z_memb %in% group)
  group_memb <- z_memb[nodes_in_group] 
  group_net <- get.inducedSubgraph(network, v = nodes_in_group)
  for (i in 1:length(group)) {
    if (sum(group_memb == group[i]) > 1) { 
      nodes_in_clust <- which(group_memb == group[i]) 
      cur_sub <- get.inducedSubgraph(group_net, v = nodes_in_clust)
      if (sum(grepl("edgecov", formula) | grepl("dyadcov", formula)) > 0) {
        for (ii in 1:length(edge_cov_name)) {
          set.network.attribute(cur_sub, edge_cov_name[[ii]], as.matrix(edge_cov_list[[i]][[ii]]))
        }
      }
      net_list[[i]] <- cur_sub
    } else { 
      net_list[[i]] <- network.initialize(1)
    }
  } 
  return(net_list)
}

estimate_params <-
  function(formula,
           network,
           parameterization,
           z_memb,
           parallel = FALSE,
           number_cores = 3,
           max_iter = 5,
           number_groups = 1,
           verbose = 1,
           initial_estimate = NULL,
           MCMCparams,
           seeds,
           sample_size_multiplier_blocks,
           NR_step_len,
           NR_step_len_multiplier,
           NR_max_iter)
  {
    n_nodes <- length(z_memb)
    block_sizes <- table(z_memb)
    max_number = length(block_sizes)
    
    params <- mlergm(form = formula, 
                     node_memb = z_memb, 
                     theta_init = initial_estimate, 
                     verbose = verbose,
                     seed = seeds[1],
                     options = set_options(
                       burnin = MCMCparams$burnin, 
                       interval = MCMCparams$interval, 
                       sample_size = MCMCparams$samplesize,
                       NR_max_iter = NR_max_iter, 
                       MCMLE_max_iter = max_iter, 
                       do_parallel = parallel, 
                       number_cores = number_cores),
                     parameterization = parameterization)

      
     return(params)
  }


#Estimate between-neighborhood parameter
estimate_bw_param <- function(network, z_memb) {
  n_nodes <- length(z_memb)
  block_sizes <- table(z_memb)
  memb_inds <- as.numeric(names(block_sizes))
  max_number = length(block_sizes)

  indicator = matrix(TRUE, n_nodes, n_nodes)
  for (i in memb_inds)
  {
    indicator[z_memb == i, z_memb == i] <- FALSE
  }
  n <- sum(indicator)
  k <- sum(network[indicator])
  bw_param = logit(k / n)/log(n_nodes)
  bw_st_error = ((1 / k) + (1 / (n - k)))/log(n_nodes)
  list(bw_param = bw_param,
       bw_st_error = bw_st_error)
}

same_clusters <- function(z_memb) {
  n_nodes <- length(z_memb)
  block_sizes <- table(z_memb)
  memb_inds <- as.numeric(names(block_sizes))
  indicator = matrix('grey', n_nodes, n_nodes)
  for (i in memb_inds)
  {
    indicator[z_memb == i, z_memb == i] <- 'black'
  }
  indicator
}

#Estimate deinsity inside clusters
estimate_density <- function(network, z_memb, bw_param) {
  n_nodes <- length(z_memb)
  block_sizes <- table(z_memb)
  memb_inds <- as.numeric(names(block_sizes))
  max_number = length(block_sizes)

  Pi = matrix(inv.logit(bw_param), max_number, max_number)
  for (i in 1:max_number) {
    Pi[i, i] = sum(network[z_memb == memb_inds[i], z_memb == memb_inds[i]]) /
      (block_sizes[i] * (block_sizes[i] - 1))
  }
  Pi
}

#Need to add split cluster function and restart
EM_wrapper_fast <-
  function(network,
           formula,
           max_number,
           n_em_step_max = 100,
           min_size = 2,
           initialization_method,
           verbose = 1)
  {
    #Step 0: Initialization

    n_nodes <- nrow(as.matrix(network))
    
    ## Calculate network statistics needed for variational approximation of p(Z|X)
    if (verbose > 0) cat("\nStep 1: Initialize z")
    stat00 <-
      stat01 <-
      stat10 <- stat11 <- matrix(as.integer(0), n_nodes, n_nodes)

    stats <- calculateStats(network, stat00, stat01, stat10, stat11)
    stat00 <- matrix(as.double(stats[[1]]), n_nodes, n_nodes)
    stat01 <- matrix(as.double(stats[[2]]), n_nodes, n_nodes)
    stat10 <- matrix(as.double(stats[[3]]), n_nodes, n_nodes)
    stat11 <- matrix(as.double(stats[[4]]), n_nodes, n_nodes)

    #Step 1a: Get initial estimate of Z memberships using spectral clustering
    if (initialization_method == 1) # Walk trap
      {
      g <- asIgraph(as.network(network))
      b <- cluster_walktrap(g, steps = 4)
      temp_ind <- factor(b$membership, levels = 1:max_number)
      z_memb_init <- temp_ind
      z_memb <-
        factor(check_clusters(temp_ind, network, max_number, min_size),
               levels = 1:max_number)
      }
    else # Spectral clustering
      {
      z_memb <- spec_clust(network, max_number)
      z_memb_init <- z_memb
      } 

    #Step 1b: Get estimate of ERGM parameters within clusters
    bw_param <- estimate_bw_param(network, z_memb)$bw_param
    Pi = estimate_density(network, z_memb, bw_param)

    #Step 2a: Find A(Z=z) ~ P(Z=z|X=x)
    if (verbose > 0)
      {
      cat(paste("\n\nStep 2: Find variational approximation A(Z=z) ~ P(Z=z|X=x)", sep = ""))
      }
    block_sizes <- table(z_memb)
    memb_inds <- as.numeric(names(block_sizes))
    max_number = length(block_sizes)

    tau <- matrix(0, n_nodes, max_number)
    tau_prev <- matrix(0, n_nodes, max_number)
    for (i in 1:n_nodes) {
      tau[i,] <- 1
      tau[i, z_memb[i]] <- 1000
      tau[i,] <- tau[i,] / sum(tau[i,])
    }
    alpha = colMeans(tau)

    delta <- row(Pi) - col(Pi)
    counter_e_step <- 0
    repeat {
      counter_e_step <- counter_e_step + 1
      tau_prev <- tau
      tau <-
        runFixedPointEstimationEStepMM(
          as.integer(n_nodes),
          as.integer(max_number),
          as.double(alpha),
          Pi,
          stat00,
          stat01,
          stat10,
          stat11,
          tau,
          network
        )

      Pi <- easy_M_Step(
        as.integer(n_nodes),
        as.integer(max_number),
        as.double(alpha),
        Pi,
        as.matrix(network),
        tau
      )

      #Permute labels
      perm <- sample(1:max_number)
      tau <- permute_tau(tau, perm)
      Pi <- Pi[perm, perm]


      alpha = colMeans(tau)
      if(counter_e_step %% 5) gc()
      if (counter_e_step >=  n_em_step_max) {
        break
      }
    }

    z_memb <-
      factor(apply(tau, 1, which.max), levels = 1:max_number)
    z_memb_final <-
      factor(check_clusters(z_memb, network, max_number, min_size),
             levels = 1:max_number)

    list(
      Pi = Pi,
      z_memb_init = z_memb_init,
      z_memb_em = z_memb,
      z_memb_final = z_memb_final
    )
  }

hergm.large <- function(network,
           formula,
           parameterization,
           max_number,
           number_cores,
           number_groups = 1,
           indicator = NULL,
           same_between_blocks = TRUE,
           estimate_parameters = TRUE,
           verbose,
           n_em_step_max = 100,
           max_iter = max_iter, 
           initialization_method,
           initial_estimate = NULL,
           MCMCparams,
           sample_size_multiplier_blocks, 
           seeds = NULL,
           NR_step_len,
           NR_step_len_multiplier,
           NR_max_iter = NR_max_iter) {

    if (number_cores == 1) {
      parallel <- FALSE
    } else {
      parallel <- TRUE
    }

    sbm_pi <- NULL
    if (is.null(indicator)) { 
      all_indicators_fixed <- FALSE
    } else { 
      all_indicators_fixed <- TRUE
    }

    if (all_indicators_fixed == FALSE) {
      network_ <- as.matrix(network)
      set.seed(seeds[1])
      answer <- EM_wrapper_fast(
        network_,
        formula,
        max_number,
        n_em_step_max = n_em_step_max,
        min_size = 1,
        initialization_method,
        verbose = verbose
      )
      indicator = answer$z_memb_final
      sbm_pi <- answer$Pi
    } else {
      if (verbose > 0) cat("\nSkipping Steps 1 and 2: z specified")
    }

    if (estimate_parameters == TRUE) {
      if (verbose > 0) {
        cat("\n\nStep 3: Estimate parameters conditional on z")
      } 
      
      est_list <-  estimate_params(
          formula,
          network,
          parameterization,
          indicator,
          initial_estimate = initial_estimate,
          parallel = parallel,
          number_cores = number_cores,
          max_iter = max_iter,
          number_groups = number_groups,
          verbose = verbose,
          MCMCparams = MCMCparams,
          seeds = seeds,
          sample_size_multiplier_blocks = sample_size_multiplier_blocks,
          NR_step_len = NR_step_len,
          NR_step_len_multiplier = NR_step_len_multiplier,
          NR_max_iter = NR_max_iter
        )
        params <- est_list; rm(est_list) 
        estimation_status <- params$estimation_status
        mcmc_path <- params$mcmc_chain 
        parameters = t(as.matrix(params$theta))
        st.error = params$se
        labels = 1:max_number

    } else { 
    # estimate_parameters = FALSE: 
     cat("\n")
     estimation_status <- "not_estimated"
     parameters <- NULL
     st.error <- NULL
     between_parameter <- NULL
     st.error.between <- NULL
     mcmc_path <- NULL
     }

    list(
      sbm.params = sbm_pi,
      partition = indicator,
      parameters = parameters,
      st.error = st.error,
      mcmc_path = mcmc_path,
      between_parameter = params$between_theta,
      st.error.between = params$between_theta,
      labels = labels,
      estimation_status = estimation_status
    )
  }

get_coef_names <- function(model_obj, is_canonical) { 
  if(is_canonical) {
    model_obj$coef.names
  } else { 
    unlist(lapply(model_obj$terms, function(term) nvl_mimic(names(term$params), term$coef.names)))
  }
}

nvl_mimic <-  function(...) {
  for (x in list(...)) {
    if (!is.null(x)) { 
      break
    }
  }
  x
}
