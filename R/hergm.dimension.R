hergm.dimension <- function(d1, d2, k, n)
# input: number of ergm terms d1, number of hergm terms d2, number of blocks k, number of nodes n
# output: number of elements stored on iteration of MCMC algorithm
{
  d <- d1 + d2
  terms <- d1 + (d2 * (k + 1)) + n + k + k + 1 + d 
  #terms <- d1 # Number of ergm terms
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
