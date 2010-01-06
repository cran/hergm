length_mcmc <- function(d1, d2, k, n)
# input: number of ergm terms d1, number of hergm terms d2, number of blocks k, number of nodes n
# output: number of elements stored on iteration of MCMC algorithm
{
  d <- d1 + d2
  terms <- d1 + (2 * d2) + (d2 * (k + 1)) + n + k + k + 1 + d 
  #terms <- d1 # Number of ergm terms
  #       + (2 * d2) # Number of mean and precision parameters of Gaussian baseline distribution of Dirichlet / stick-breaking prior
  #       + ((d2 + 1) * k) # Number of hergm terms
  #       + n # Number of category indicators  
  #       + k # Number of category sizes 
  #       + k # Number of category probabilities
  #       + 1 # Clustering parameter
  #       + d # Number of posterior predictions
  #print("terms")
  #print(terms)
  terms
}

number_clusters <- function(alpha, n)
# input: scaling parameter of Dirichlet process <alpha>, number of units <n>
# output: first two moments of number of non-empty clusters
# see: Teh (2010), Encyclopedia of Machine Learning
{
  mean <- alpha * (digamma(alpha + n) - digamma(alpha))
  variance <- mean + (alpha * alpha) * (trigamma(alpha + n) - trigamma(alpha))
  moments <- list()
  moments$mean <- mean
  moments$variance <- variance
  moments
}

bernoulli_map_mean_to_natural <- function(n, mu, directed)
# input: number of nodes <n>, mean-value parameter of Bernoulli random graph model <mu>, indicator of whether network is directed
# output: natural parameter of Bernoulli random graph model
{
  if (directed == FALSE) df <- n * (n - 1) / 2
  else df <- n * (n - 1)
  theta <- - log((df - mu) / mu)
  theta
}

