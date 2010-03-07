hergm.postprocess <- function(n = NULL,
                              k = NULL, 
                              d1 = 0, 
                              d2 = 0, 
                              parallel = 1,
                              burnin = NULL, 
                              samplesize = NULL, 
                              mcmc = NULL, 
                              relabel = TRUE,
                              output = TRUE,
                              name = NULL, 
                              ...)
# input: number of nodes n, number of categories k, number of non-hierarchical terms d1, number of hierarchical terms d2, 
# number of burn-in iterations, number of draws from the posterior (including number of burn-in iterations), 
# MCMC sample
# output: postprocessed MCMC sample
{
  
  # Check arguments
  d <- d1 + d2
  terms <- length_mcmc(d1, d2, k, n)  
  if (burnin > 12000) cat("\nWarning: number of burn-in iterations specified as ", burnin, ", but function hergm cannot return MCMC samples of size > 12000.\n", sep = "")
  if (samplesize > 12000) cat("\nWarning: MCMC sample size specified as ", samplesize, ", but function hergm cannot return MCMC samples of size > 12000.\n", sep = "")

  if (length(mcmc) != (parallel * samplesize * terms))
    {
    cat("Arguments:")    
    cat("\n- length(mcmc) = ", length(mcmc))
    cat("\n- parallel = ", parallel)
    cat("\n- samplesize = ", samplesize, " (note: if function hergm is called with samplesize > 12000, hergm returns 12000 draws)")
    cat("\n- n = ", n)
    cat("\n- k = ", k)
    cat("\n- d1 = ", d1)
    cat("\n- d2 = ", d2)
    cat("\n- terms = ", terms)
    cat("\n")
    error_message <- paste("hergm.postprocess: arguments inconsistent: length(mcmc) is not as expected, indicating that either mcmc or other arguments are incorrect.")
    stop(error_message, call. = FALSE)
    }
  if (k > 10) 
    {
    if (relabel == TRUE) cat("\nWarning: relabeling too time-consuming: skipping relabeling.\n")
    relabel <- FALSE
    }

  # Preprocess MCMC sample: delete burn-in iterations and transform vector into matrix, where rows correspond to MCMC draws
  mcmc_sample <- NULL
  count <- 0
  for (i in 1:parallel) 
    {
    first <- count + (burnin * terms) + 1
    last <- count + (samplesize * terms)
    mcmc_sample <- append(mcmc_sample, mcmc[first:last])
    count <- count + (samplesize * terms)
    }
  mcmc_samplesize <- parallel * (samplesize - burnin)
  if (mcmc_samplesize <= 0) 
    {
    cat("\n")
    error_message <- paste("number of burn-in iterations exceeds number of post-burn-in iterations.")
    stop(error_message, call. = FALSE)
    }
  mcmc_sample <- matrix(mcmc_sample, nrow = mcmc_samplesize, ncol = terms, byrow = TRUE)

  # Initialize arrays
  s <- list()
  s$theta <- matrix(data = 0, nrow = mcmc_samplesize, ncol = d1) 
  s$htheta_mean <- matrix(data = 0, nrow = mcmc_samplesize, ncol = d2)
  s$htheta_precision <- matrix(data = 0, nrow = mcmc_samplesize, ncol = d2)
  s$htheta <- matrix(data = 0, nrow = mcmc_samplesize, ncol = d2 * (k + 1))
  s$indicator <- matrix(data = 0, nrow = mcmc_samplesize, ncol = n)
  s$size <- matrix(data = 0, nrow = mcmc_samplesize, ncol = k)
  s$p <- matrix(data = 0, nrow = mcmc_samplesize, ncol = k)
  s$alpha <- matrix(data = 0, nrow = mcmc_samplesize, ncol = 1)
  s$pp <- matrix(data = 0, nrow = mcmc_samplesize, ncol = d)

  # Process MCMC sample
  for (row in 1:mcmc_samplesize)
    {
    column <- 0
    if (d1 > 0)
      {
      for (i in 1:d1) 
        {
        column <- column + 1
        s$theta[row,i] <- mcmc_sample[row,column]
        }
      }
    for (i in 1:d2) 
      {
      column <- column + 1
      s$htheta_mean[row,i] <- mcmc_sample[row,column]
      }
    for (i in 1:d2) 
      {
      column <- column + 1
      s$htheta_precision[row,i] <- mcmc_sample[row,column]
      }
    for (i in 1:(d2 * (k + 1))) 
      {
      column <- column + 1
      s$htheta[row,i] <- mcmc_sample[row,column]
      }
    for (i in 1:n) 
      {
      column <- column + 1
      s$indicator[row,i] <- mcmc_sample[row,column]
      }
    for (i in 1:k) 
      {
      column <- column + 1
      s$size[row,i] <- mcmc_sample[row,column]
      }
    for (i in 1:k) 
      {
      column <- column + 1
      s$p[row,i] <- mcmc_sample[row,column]
      }
    column <- column + 1
    s$alpha[row,1] <- mcmc_sample[row,column]
    for (i in 1:d) 
      {
      column <- column + 1
      s$pp[row,i] <- mcmc_sample[row,column]
      }
    }

  if (output == TRUE)
    {

    # Write MCMC sample to files
    if (is.null(name)) name <- "hergm"
    if (d1 > 0) write(t(s$theta), paste(sep = "",name,"_parameter.out"), ncolumns = d1)
    write(t(s$htheta_mean), paste(sep = "",name,"_mean_block_parameter.out"), ncolumns = d2)
    write(t(s$htheta_precision), paste(sep = "",name,"_precision_block_parameter.out"), ncolumns = d2)
    write(t(s$htheta), paste(sep = "",name,"_block_parameter.out"), ncolumns = d2 * (k + 1))
    write(t(s$indicator), paste(sep = "",name,"_indicator.out"), ncolumns = n)
    write(t(s$size), paste(sep = "",name,"_size.out"), ncolumns = k)
    write(t(s$p), paste(sep = "",name,"_p.out"), ncolumns = k)
    write(t(s$alpha), paste(sep = "",name,"_alpha.out"), ncolumns = 1)
    write(t(s$pp), paste(sep = "", name, "_statistics.out"), ncolumns = d)

    # Scaling parameter
    filename <- paste(sep = "", name, "_alpha.pdf")
    pdf(filename)
    plot(density(s$alpha, from = 0), main = "", xlab = "", ylab = "")
    abline(v = c(quantile(s$alpha, 0.05), quantile(s$alpha, 0.95)))
    abline(v = 0, lty = 2)
    dev.off()

    # Means of block parameters
    for (i in 1:d2) 
      {
      filename <- paste(sep = "", name, "_mean_block_parameter_", i, ".pdf")
      pdf(filename)
      plot(density(s$htheta_mean[,i]), main = "", xlab = "", ylab = "")
      abline(v = c(quantile(s$htheta_mean[,i], 0.05), quantile(s$htheta_mean[,i], 0.95)))
      abline(v = 0, lty = 2)
      dev.off()
      }

    # Precisions of block parameters
    for (i in 1:d2) 
      {
      filename <- paste(sep = "", name, "_precision_block_parameter_", i, ".pdf")
      pdf(filename)
      plot(density(s$htheta_precision[,i], from = 0), main = "", xlab = "", ylab = "")
      abline(v = c(quantile(s$htheta_precision[,i], 0.05), quantile(s$htheta_precision[,i], 0.95)))
      abline(v = 0, lty = 2)
      dev.off()
      }

    # Relabel sample
    if (relabel == TRUE)
      {
      minimizer <- hergm.min_loss(k, paste(sep = "", name, "_indicator.out"), 0, 100) # Specify number of iterations of post-processing algorithm
      s$q <- minimizer$p
      # hergm.pie(name, n, minimizer$p)
      index <- 0
      for (h_term in 1:d2)
        {
        index <- index + 1
        theta <- s$htheta[,index]
        for (i in 2:k) 
          {
          index <- index + 1
          theta <- cbind(theta, s$htheta[,index])
          }
        write(t(theta), paste(sep = "", name, "_block_parameter_", h_term, ".out"), ncolumns = k)
        min_theta <- hergm.permute_mcmc(theta, k, minimizer$min_permutations) 
        write(t(min_theta), paste(sep = "", name, "_block_parameter_min_", h_term, ".out"), ncolumns = k)
        for (i in 1:k)
          {
          filename <- paste(sep = "", name, "_block_parameter_", h_term, "_", i, ".pdf")
          pdf(filename)
          plot(density(min_theta[,i]), main = paste(sep = "", i), xlab = "", ylab = "", cex.main = 2)
          abline(v = c(quantile(min_theta[,i], 0.05), quantile(min_theta[,i], 0.95)))
          abline(v = 0, lty = 2)
          dev.off()
          }
        index <- index + 1
        filename <- paste(sep = "", name, "_block_parameter_", h_term, "_between.pdf")
        pdf(filename)
        plot(density(s$htheta[,index]), main = "", xlab = "", ylab = "")
        abline(v = c(quantile(s$htheta[,index], 0.05), quantile(s$htheta[,index], 0.95)))
        abline(v = 0, lty = 2)
        dev.off()
        }
      }

    # Parameters
    if (d1 > 0)
      {
      for (i in 1:d1)
        {
        filename <- paste(sep = "", name, "_parameter_", i, ".pdf")
        pdf(filename)
        if (d1 == 1) 
          {
          plot(density(s$theta), main = paste(sep = "", 1), xlab = "", ylab = "")
          abline(v = c(quantile(s$theta, 0.05), quantile(s$theta, 0.95)))
          abline(v = 0, lty = 2)
          }
        else 
          {
          plot(density(s$theta[,i]), main = paste(sep = "", i), xlab = "", ylab = "")
          abline(v = c(quantile(s$theta[,i], 0.05), quantile(s$theta[,i], 0.95)))
          abline(v = 0, lty = 2)
          }
        dev.off()
        }
      }

    # Number of relevant categories
    filename <- paste(sep = "",name,"_number_1.pdf")
    pdf(filename)
    number <- vector(length = nrow(s$size))
    for (i in 1:nrow(s$size)) number[i] <- sum(s$size[i,] >= 1)
    hist(number, 50, xlab = "", ylab = "", main = "")
    dev.off()
    filename <- paste(sep = "", name, "_number_3.pdf")
    pdf(filename)
    number <- vector(length = nrow(s$size))
    for (i in 1:nrow(s$size)) number[i] <- sum(s$size[i,] >= 3)
    hist(number, 50, xlab = "", ylab = "", main = "")
    dev.off()

    # Scatter plots
    pdf(paste(sep = "", name, "_statistics.pdf"))
    if (d == 1)
      {
      hist(s$pp, 50, main = "", xlab = "", ylab = "")
      }
    else
      {
      pairs(s$pp)
      }
    dev.off()

    }

  mean_deviance <- mean(s$htheta[,k+1]) # MCMC sample average of deviance is stored in column corresponding to first hierarchical term, block k+1
  # cat("\nMean deviance: ",mean_deviance,"\n")

  cat("\n")

  s

}

