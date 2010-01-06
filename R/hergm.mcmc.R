hergm.mcmc <- function(nw, model, MHproposal, MCMCparams, verbose, name, alpha_shape, alpha_rate, alpha, eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, parallel, simulate, seeds, mh_scale, output) 
{
  if (simulate == FALSE) 
    {
    t <- system.time(scalefactor <- hergm.set.mcmc(nw, model, MHproposal, MCMCparams, verbose, name, alpha_shape, alpha_rate, alpha,  eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, parallel, simulate, seeds, output, mh_scale)) # The last argument is the initial value of the scale factor
    cat("\n")
    print(t)
    }

  # Prepare
  Clist <- ergm.Cprepare(nw, model)
  if (Clist$dir == FALSE) maxedges <- Clist$n * (Clist$n - 1) / 2 # Undirected
  else maxedges <- Clist$n * (Clist$n - 1) # Directed
  hergm_list <- hergm.preprocess(nw, model, Clist, MHproposal, MCMCparams, maxedges, alpha_shape, alpha_rate, alpha,  eta_mean_mean, eta_mean_sd, eta_precision_shape, eta_precision_rate, eta_mean, eta_sd, eta, simulate, parallel, output, name, verbose)
  sample <- list()
  sample$newnwheads = maxedges + 1
  sample$mcmc = length(hergm_list$mcmc)
  sample$sample_heads = length(hergm_list$sample_heads)
  sample$sample_tails = length(hergm_list$sample_tails)
  if (simulate == TRUE) 
    {
    k <- number_clusters(hergm_list$alpha, hergm_list$Clist$n)
    cat("\nDirichlet process:")
    cat("\n- scaling parameter: ", hergm_list$alpha, sep = "")
    cat("\n- mean of number of non-empty blocks: ", formatC(k$mean, digits = 2, width = 4, format = "f", mode = "real"), sep = "") 
    cat("\n- variance of number of non-empty blocks: ", formatC(k$variance, digits = 2, width = 4, format = "f", mode = "real"), sep = "") 
    cat("\n")   
    }
  else
    {
    hergm_list$scalefactor <- scalefactor # Set scale factor
    if (verbose >= 0) cat("\nFinal scale factor:", formatC(hergm_list$scalefactor, digits = 4, width = 6, format = "f", mode = "real"), "\n")
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
  else 
    {
    sample <- hergm.wrapper(seeds[1], hergm_list) # Non-parallel computing
    if (simulate == FALSE) sample <- sample$mcmc
    }

  # Postprocess
  samplesize <- hergm_list$MCMCparams$samplesize
  if (hergm_list$simulate == TRUE) # Simulation
    {
    if (hergm_list$output == TRUE) # Full output
      {
      file_heads <- paste(sep = "", name, "_heads.out")
      file_tails <- paste(sep = "", name, "_tails.out")
      file.create(file_heads)
      file.create(file_tails)  
      h <- 0
      for (i in 1:samplesize)
        { 
        h <- h + 1
        number_edges <- sample$sample_heads[h]  # First element of sample$newnwheads and sample$newnwtails is number of edges 
        first <- h # First element corresponds to number of edges
        last <- h + number_edges
        write(sample$sample_heads[first:last], file_heads, ncolumns = number_edges + 1, append = TRUE)
        write(sample$sample_tails[first:last], file_tails, ncolumns = number_edges + 1, append = TRUE)
        h <- last
        }
      output <- sample$mcmc
      }
    else # Limited output
      {
      sample <- matrix(data = 0, nrow = samplesize, ncol = hergm_list$d)
      for (i in 1:samplesize)
        {
        for (j in 1:hergm_list$d)
          {
          index <- ((i - 1) * samplesize * hergm_list$terms) + (hergm_list$terms - hergm_list$d) + j
          sample[i,j] <- hergm_list$mcmc[index]
          }
        }
      output <- sample
      }    
    }
  else # Bayesian inference
    {
    if (hergm_list$output == TRUE)
      {
      if (is.null(name)) name <- "hergm"
      filename <- paste(sep = "", name, "_mcmc.out")
      write(sample, filename, ncolumns = hergm_list$terms)
      }
    output <- sample
    }

  if (hergm_list$output == TRUE) cat("\n")
  else cat("\n\n")

  output
}

