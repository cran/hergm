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

silent <- function(f) 
# Shup up loud function f that insists on printing to the screen
{ 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(f)) 
}

summary_sample_network <- function(edgelists, sample_size, directed, n)
{
  # Simulated networks
  output <- list()
  output$component.number <- vector(length = sample_size)
  output$max.component.size <- vector(length = sample_size)
  output$distance.label <- matrix(0, nrow = sample_size, ncol = n)
  output$distance <- matrix(0, nrow = sample_size, ncol = n)
  output$edges <- vector(length = sample_size)
  output$degree <- matrix(0, nrow = sample_size, ncol = n)
  output$stars <- vector(length = sample_size)
  output$triangles <- vector(length = sample_size)
  for (i in 1:sample_size)
    {
    edgelist <- edgelists[[i]] # Edge list of simulated network
    simulated.network <- as.network(edgelist, directed = directed, matrix.type = "edgelist") # The simulated network as simulated.network; note: simulated_network$gal$n <- maximum vertex number
    if (simulated.network$gal$n < n) add.vertices(simulated.network, nv = n-simulated.network$gal$n) # If simulated$gal$n < n, add isolates 
    components <- silent(component.dist(simulated.network)) # Problem: in sna/src/components.c, output is printed that cannot be surppressed.
    output$component.number[i] <- length(components$csize) # Number of components
    output$max.component.size[i] <- max(components$csize) # Size of largest component
    output$distances <- geodist(simulated.network)
    output$distances <- output$distances$gdist
    output$frequencies_distances <- table(output$distances) # First column: frequency of self loops; columns 2:number: frequencies finite and (last column) infinite distances
    output$number <- length(output$frequencies_distances) - 1 - (sum(output$distances == Inf) > 0) # Number of distances minus 0-distance minus Inf-distance
    output$distance.label[i,1:output$number] <- rownames(output$frequencies_distances)[2:(output$number+1)]
    output$distance[i,1:output$number] <- output$frequencies_distances[2:(output$number+1)] # Frequencies of finite distances
    output$edges[i] <- summary(simulated.network ~ edges)
    if (is.directed(simulated.network)) # Directed networks
      {
      output$degree[i,] <- summary(simulated.network ~ odegree(1:n-1)) # Degree distribution             
      output$stars[i] <- summary(simulated.network ~ ostar(2))
      output$triangles[i] <- summary(simulated.network ~ ttriple)
      }
    else # Undirected networks
      {
      output$degree[i,] <- summary(simulated.network ~ degree(1:n-1)) # Degree distribution             
      output$stars[i] <- summary(simulated.network ~ kstar(2))
      output$triangles[i] <- summary(simulated.network ~ triangles)
      }
    }
  output
}

plot_weighted_graph <- function(network, p) 
{
  n <- network$gal$n
  for (i in 1:n)
    {
    for (j in 1:n)
      {
      if (p[i,j] <= .5) p[i,j] <- 0
      }
    }
  edge.col<- matrix(paste0("grey", round(exp(10*p) / exp(10), digits = 2) * 100), nrow=nrow(p), ncol=ncol(p)) # Create a color matrix based on various hues of greys
  # Plot
  coordinates <- gplot(network, gmode="graph", mode="fruchtermanreingold", vertex.cex=1, vertex.col=0, vertex.border=0, displaylabels=TRUE, label=c(1:n), label.cex=0.8)
  gplot(as.network(p), coord=coordinates, edge.col=edge.col, gmode="graph")
}

