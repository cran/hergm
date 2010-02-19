hergm.postprocess <- function(n = NULL,
                              k = NULL, 
                              d1 = 0, 
                              d2 = 0, 
                              burnin = NULL, 
                              samplesize = NULL, 
                              mcmc = NULL, 
                              output = TRUE,
                              name = NULL, 
                              ...)
# input: number of nodes n, number of categories k, number of non-hierarchical terms d1, number of hierarchical terms d2, 
# number of burn-in iterations, number of draws from the posterior (including number of burn-in iterations), 
# MCMC sample
# output: postprocessed MCMC sample
{
  d <- d1 + d2
  terms <- hergm.dimension(d1, d2, k, n)
  mcmc <- matrix(mcmc, nrow = samplesize, ncol = terms, byrow = TRUE)

  # Initialize arrays
  s <- list()
  s$theta <- matrix(data = 0, nrow = samplesize, ncol = d1) 
  s$htheta <- matrix(data = 0, nrow = samplesize, ncol = d2 * (k + 1))
  s$indicator <- matrix(data = 0, nrow = samplesize, ncol = n)
  s$size <- matrix(data = 0, nrow = samplesize, ncol = k)
  s$p <- matrix(data = 0, nrow = samplesize, ncol = k)
  s$alpha <- matrix(data = 0, nrow = samplesize, ncol = 1)
  s$pp <- matrix(data = 0, nrow = samplesize, ncol = d)

  # Process MCMC sample
  for (row in 1:samplesize)
    {
    column <- 0
    if (d1 > 0)
      {
      for (i in 1:d1) 
        {
        column <- column + 1
        s$theta[row,i] <- mcmc[row,column]
        }
      }
    for (i in 1:(d2 * (k + 1))) 
      {
      column <- column + 1
      s$htheta[row,i] <- mcmc[row,column]
      }
    for (i in 1:n) 
      {
      column <- column + 1
      s$indicator[row,i] <- mcmc[row,column]
      }
    for (i in 1:k) 
      {
      column <- column + 1
      s$size[row,i] <- mcmc[row,column]
      }
    for (i in 1:k) 
      {
      column <- column + 1
      s$p[row,i] <- mcmc[row,column]
      }
    column <- column + 1
    s$alpha[row,1] <- mcmc[row,column]
    for (i in 1:d) 
      {
      column <- column + 1
      s$pp[row,i] <- mcmc[row,column]
      }
    }

  if (output == TRUE)
    {
    # Write MCMC sample to files
    if (is.null(name)) name <- "hergm"
    if (d1 > 0) write(t(s$theta), paste(sep = "",name,"_parameter.out"), ncolumns = d1)
    write(t(s$htheta), paste(sep = "",name,"_block_parameter.out"), ncolumns = d2 * (k + 1))
    write(t(s$indicator), paste(sep = "",name,"_indicator.out"), ncolumns = n)
    write(t(s$size), paste(sep = "",name,"_size.out"), ncolumns = k)
    write(t(s$p), paste(sep = "",name,"_p.out"), ncolumns = k)
    write(t(s$alpha), paste(sep = "",name,"_alpha.out"), ncolumns = 1)
    write(t(s$pp), paste(sep = "", name, "_statistics.out"), ncolumns = d)
    }

  if (burnin > 0)
    {
    # Discard burnin iterations for label-independent entities
    for (i in 1:burnin) 
      { 
      s$alpha <- s$alpha[-1]
      if (d1 > 0)
        {
        if (d1 == 1) s$theta <- s$theta[-1]
        else s$theta <- s$theta[-1,]
        }
      s$htheta <- s$htheta[-1,]
      s$size <- s$size[-1,]
      if (d == 1) s$pp <- s$pp[-1] 
      else s$pp <- s$pp[-1,] 
      }
    }
 
  if (output == TRUE)
    {
    # Scaling parameter
    filename <- paste(sep = "", name, "_alpha.pdf")
    pdf(filename)
    plot(density(s$alpha, from = 0), main = "", xlab = "", ylab = "")
    dev.off()

    # Relabel sample
    if (k <= 20) 
      {
      minimizer <- hergm.min_loss(k, paste(sep = "", name, "_indicator.out"), burnin, 100) # Specify number of iterations of post-processing algorithm
      s$q <- minimizer$p
      hergm.pie(name, n, minimizer$p)
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
    else cat("Sample not relabeld: number of permutations too large.")

    # Parameters
    if (d1 > 0)
      {
      for (i in 1:d1)
        {
        filename <- paste(sep = "", name, "_parameter_", i, ".pdf")
        pdf(filename)
        if (d1 == 1) plot(density(s$theta), main = paste(sep = "", 1), xlab = "", ylab = "")
        else plot(density(s$theta[,i]), main = paste(sep = "", i), xlab = "", ylab = "")
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

    cat("\n")

  s
}

