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

hergm.loss <- function(n_nodes, z, p)
# Loss function (Stephens, 2000, discussion): minus log joint classification probability
# input: number of nodes, categories of nodes, classification probabilities
# output: loss 
{
  loss <- 0
  for (i in 1:n_nodes) # Node
    {
    k <- z[i] # Category of node
    loss <- loss - log(p[i,k]) # Loss associated with classification of node to category
    }
  loss 
}

hergm.permute <- function(n_nodes, z, nu)
# Permute labels of category indicators z in accordance with permutation nu
# input: number of nodes, categories of nodes, permutation
# output: permuted categories of nodes
{
  nu_z <- vector(mode = "numeric", length = n_nodes)
  for (i in 1:n_nodes)
    {
    k <- z[i] # Present value
    nu_z[i] <- nu[k] # Permuted value
    }
  nu_z
}

hergm.permute_indicator <- function(n_nodes, n_categories, n_sample, indicator, permutations)
{
  nu_indicator <- matrix(data = 0, nrow = n_sample, ncol = n_nodes) # Classification probabilities
  z <- vector(mode = "numeric", length = n_nodes)
  nu <- vector(mode = "numeric", length = n_categories) 
  nu_z <- vector(mode = "numeric", length = n_nodes)
  for (i in 1:n_sample)
    {
    z <- indicator[i,]
    nu <- permutations[i,] # Set nu to minimizing permutation
    nu_indicator[i,] <- hergm.permute(n_nodes,z,nu) # Permute category indicators z in accordance with permutation nu
    }
  nu_indicator
}

hergm.sa <- function(n_nodes, n_categories, n_permutations, permutations, z, p, last_minimizer) 
# Simulated Annealing
# input: number of nodes, number categories, number of permutations, MCMC sample point of categories of nodes, classification probabilities, last minimizer
# output: sample minimum, minimizer
{
  nu <- permutations[last_minimizer,] # Set nu to i-th permutation
  nu_z <- hergm.permute(n_nodes,z,nu) # Permute z in accordance with nu
  last_loss <- hergm.loss(n_nodes,nu_z,p) # Compute loss under nu_z and p
  last <- last_minimizer # Permutation
  minimum <- last_loss # Minimum of loss function
  minimizer <- last_minimizer # Minimizer of loss function
  t <- 2 # Temperature
  factor <- 0.9 # Factor by which temperature is decreased
  n_iterations <- n_categories # Number of iterations
  while (t > 0.1) # Number of batches with constant temperature
    {
    t <- factor * t # Decrease temperature
    for (iteration in 1:n_iterations) # Number of iterations per batch with constant temperature 
      {
      # Candidate: walk on "ring"
      u <- runif(1, min = -10, max = 10) # Uniform(-10,10)
      i <- last + round(u)
      if (i < 1) i <- i + n_permutations 
      else if (i > n_permutations) i <- i - n_permutations 
      nu <- permutations[i,] # Set nu to i-th permutation
      nu_z <- hergm.permute(n_nodes,z,nu) # Permute z in accordance with nu
      loss <- hergm.loss(n_nodes,nu_z,p) # Compute loss under nu_z and p
      # To accept or not
      if (is.infinite(last_loss)) last_loss <- 1e5
      if (is.infinite(loss)) loss <- 1e5
      d <- loss - last_loss 
      if (d < 0) accept <- 1
      else 
        {
        pr <- exp(- d / t)
        if (u < pr) accept <- 1
        else accept <- 0
        }
      # If accepted, store
      if (accept == 1)
        {
        if (loss < minimum) # Store sample minimum and minimizer
          {
          minimum <- loss
          minimizer <- i
          }
        last_loss <- loss # Store loss
        last <- i # Store permutation
        }
      }
    }
  SA <- list()
  SA$minimum <- minimum
  SA$minimizer <- minimizer
  SA
}

hergm.step_1 <- function(n_nodes, n_categories, n_sample, nu_indicator)
# Step 1 of minimization algorithm
# input: number of nodes, number of categories, MCMC sample size, MCMC sample of categories of nodes
# output: classification probabilities which minimize loss function
{
  p <- matrix(data = 0, nrow = n_nodes, ncol = n_categories) # Classification probabilities
  for (i in 1:n_nodes)
    {
    for (k in 1:n_categories)
      {
      p[i,k] <- sum(nu_indicator[,i] == k) / n_sample # Sample proportion of times node i is member of category k
      }
    }
  p
}

hergm.step_2 <- function(n_sample, n_nodes, n_categories, n_permutations, permutations, indicator, p, last_minimizers)
# Step 2 of minimization algorithm
# input: MCMC sample size, number of nodes, number categories, number of permutations, permutations, MCMC sample of categories of nodes, classification probabilities
# output: permutations which minimize loss function
{
  min_permutations <- matrix(data = 0, nrow = n_sample, ncol = n_categories)
  z <- vector(mode = "numeric", length = n_nodes)
  nu <- vector(mode = "numeric", length = n_categories) 
  present_loss <- 0
  for (h in 1:n_sample) 
    {
    z <- indicator[h,]
    if (n_categories <= 5) # If number of categories small: complete enumeration
      {
      # cat("Relabeling algorithm with complete enumeration")
      min_loss <- Inf
      for (i in 1:n_permutations) 
        { 
        nu <- permutations[i,] # Set nu to i-th permutation
        nu_z <- hergm.permute(n_nodes,z,nu) # Permute z in accordance with nu
        loss <- hergm.loss(n_nodes,nu_z,p) # Compute loss under nu_z and p
        if (loss < min_loss)
          {
          min_loss <- loss
          minimizer <- i
          }
        }
      }
    else # If number of categories large: simulated annealing
      {
      # cat("\nRelabeling algorithm with simulated annealing")
      last_minimizer <- last_minimizers[h]
      m <- hergm.sa(n_nodes,n_categories,n_permutations,permutations,z,p,last_minimizer) 
      min_loss <- m$minimum
      minimizer <- m$minimizer
      last_minimizers[h] <- minimizer
      }
    present_loss <- present_loss + min_loss # Value of loss function at the end of step 2
    min_permutations[h,] <- permutations[minimizer,] # Store minimizing permutation for h-th MCMC sample point
    }
  step_2 <- list()
  step_2$present_loss <- present_loss
  step_2$min_permutations <- as.matrix(min_permutations)
  step_2
}

hergm.min_loss <- function(n_categories, filename, n_burnin, stop_criterion)
# Minimize posterior expected loss
# input: number of categories, convergence criterion, stop criterion
# output: minimizing permutations, classification probabilities, relabeled MCMC output of categories of nodes
{
  cat("\n\nRelabeling algorithm")
  cat("\n--------------------")
  n_permutations <- factorial(n_categories)
  permutations <- hergm.permutation.wrapper(n_categories) # Generate possible permutations of the category labels
  input <- as.matrix(read.table(filename)) # MCMC sample of category indicators
  if (n_burnin >= 1)
    {
    for (i in 1:n_burnin) # Category indicators: discard burn-in iterations
      { 
      input <- input[-1,] 
      }
    }
  indicator <- input
  if (min(indicator) == 0) # Category indicators: if category labels are 0..n_categories-1, then translate by 1 to obtain category labels 1..n_categories
    {
    indicator <- indicator + 1 
    # cat("\nInputted category indicators translated by 1")
    }
  n_sample <- nrow(indicator) # MCMC sample size
  n_nodes <- ncol(indicator) # Number of nodes
  p <- matrix(data = 0, nrow = n_nodes, ncol = n_categories) # Classification probabilities
  min_permutations <- matrix(data = 0, nrow = n_sample, ncol = n_categories)
  for (i in 1:n_sample) 
    {
    k <- round(runif(n = 1, min = 1, max = n_permutations)) # Sample permutation
    min_permutations[i,] <- permutations[k,] 
    }
  min_min_loss <- Inf
  min_min_permutations <- min_permutations
  last_minimizers <- vector(mode = "numeric", length = n_sample) # Minimizers
  last_minimizers[] <- 1
  last_loss <- Inf 
  change <- Inf
  i <- 0
  while (i < stop_criterion) # Stop criterion 
    {
    i <- i + 1
    cat("\n\nIteration ",i)
    nu_indicator <- hergm.permute_indicator(n_nodes,n_categories,n_sample,indicator,min_permutations)
    p <- hergm.step_1(n_nodes,n_categories,n_sample,nu_indicator) # Minimize loss with respect to classification probabilites
    step_2 <- hergm.step_2(n_sample,n_nodes,n_categories,n_permutations,permutations,indicator,p,last_minimizers) # Minimize loss with respect to permutations of category labels
    min_permutations <- step_2$min_permutations
    present_loss <- step_2$present_loss
    if (present_loss < min_min_loss) 
      {
      min_min_loss <- present_loss
      min_min_permutations <- min_permutations
      }
    cat("\nLoss:",present_loss)
    if (is.finite(present_loss - last_loss)) change <- present_loss - last_loss
    else change <- -Inf
    cat("\nChange",change)
    last_loss <- present_loss
    if (abs(change) <= 0.01) stop_criterion <- i
    }
  cat("\n\nMinimum loss:",min_min_loss)
  min_permutations <- min_min_permutations
  indicator <- hergm.permute_indicator(n_nodes,n_categories,n_sample,indicator,min_permutations)
  cat("\n\nNumber of iterations:",i)
  cat("\n\nMCMC sample relabeled\n")
  name <- sub(pattern = "_indicator.out", replacement = "", filename)
  x <- c("Local minimum of loss function:", min_min_loss)
  write(x, paste(sep = "", name, "_min_loss.out"), ncolumns = 2)
  write(t(p), paste(sep = "", name, "_classification.out"), ncolumns = n_categories)
  write(t(indicator), paste(sep = "", name, "_indicator_min.out"), ncolumns = n_nodes)
  minimizer <- list()
  minimizer$min_permutations <- min_permutations
  minimizer$p <- p
  minimizer$indicator <- indicator
  minimizer
}	

hergm.permute_mcmc <- function(mcmc, n_categories, min_permutations)
# Relabel MCMC sample
# input: name of file containing MCMC sample,  minimizing permutations
# output: MCMC sample relabeled
{
  min_mcmc <- matrix(data = 0, nrow = nrow(mcmc), ncol = n_categories)
  for (h in 1:nrow(mcmc))
    {
    for (i in 1:n_categories)
      {
      k <- min_permutations[h,i] # Present permutation: identity permutation 
      min_mcmc[h,k] <- mcmc[h,i] # Minimizing permutation
      }
    }
  min_mcmc
}

hergm.pie <- function(filename, n_nodes, p)
# Make one pie chart for each node, representing its classification probabilities
# input: name of file to which output ought to be written, number of nodes, classification probabilities
# output: pie charts
{
  for (i in 1:n_nodes)
    {
    file_i <- paste(filename, "_", i, ".pdf", sep = "")
    pdf(file_i)  
    pie(p[i,], radius = 1, clockwise = TRUE,  col = c("white", "orange", "brown", "red", "black"), main = "", labels = c(1:100), cex = 3.5)
    dev.off()
    }
}



