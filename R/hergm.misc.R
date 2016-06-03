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

length_mcmc <- function(d1, d2, k, n, predictions)
# input: number of ergm terms d1, number of hergm terms d2, number of blocks k, number of nodes n
# output: number of elements stored on iteration of MCMC algorithm
{
  terms <- d1 + (2 * d2) + (d2 * (k + 1)) 
  if (d2 > 0) terms <- terms + n + k + k + 1 
  if (predictions == TRUE) terms <- terms + d1 + d2
  #terms <- d1 # Number of ergm terms
  #       + (2 * d2) # Number of mean and precision parameters of Gaussian baseline distribution of Dirichlet / stick-breaking prior
  #       + ((d2 + 1) * k) # Number of hergm terms
  #       + n # Number of category indicators  
  #       + k # Number of category sizes 
  #       + k # Number of category probabilities
  #       + 1 # Clustering parameter
  #       + d1 + d2 # Number of posterior predictions
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

summary_sample_network <- function(sample, n, output, i)
# input: sample network in the form of edge list, number of nodes, output, number of sample networks
# output: summary of sample network
{
  sample.edgelist <- cbind(sample$heads, sample$tails) # Edge list of simulated network
  sample.network <- as.network(sample.edgelist, directed = FALSE, matrix.type = "edgelist") # The simulated network as network object; note: simulated_network$gal$n <- maximum vertex number
  if (sample.network$gal$n < n) add.vertices(sample.network, nv = n-sample.network$gal$n) # If simulated$gal$n < n, add isolates 
  components <- component.dist(sample.network)
  output$component.number[i] <- length(components$csize) # Number of components
  output$max.component.size[i] <- max(components$csize) # Size of largest component
  output$distances <- geodist(sample.network)
  output$distances <- output$distances$gdist
  output$frequencies_distances <- table(output$distances) # First column: frequency of self loops; columns 2:number: frequencies finite and (last column) infinite distances
  output$number <- length(output$frequencies_distances) - 1 - (sum(output$distances == Inf) > 0) # Number of distances minus 0-distance minus Inf-distance
  output$distance.label[i,1:output$number] <- rownames(output$frequencies_distances)[2:(output$number+1)]
  output$distance[i,1:output$number] <- output$frequencies_distances[2:(output$number+1)] # Frequencies of finite distances
  if (is.directed(sample.network)) # Directed networks
    {
    output$edges[i] <- summary(sample.network ~ edges)
    output$degree[i,] <- summary(sample.network ~ odegree(1:n-1)) # Degree distribution             
    output$stars[i] <- summary(sample.network ~ ostar(2))
    output$triangles[i] <- summary(sample.network ~ ttriple)
    }
  else # Undirected networks
    {
    output$edges[i] <- summary(sample.network ~ edges)
    output$degree[i,] <- summary(sample.network ~ degree(1:n-1)) # Degree distribution             
    output$stars[i] <- summary(sample.network ~ kstar(2))
    output$triangles[i] <- summary(sample.network ~ triangles)
    }
  output
}

