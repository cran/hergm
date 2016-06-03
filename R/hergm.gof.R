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

gof <- function(sample = NULL,
                      ...)
 {
 UseMethod("gof")
 }

gof.hergm <- function(sample = NULL, 
                      ...)
# input: postprocessed sample, number of nodes, number of blocks, number covariates, observed network
# output: goodness of fit data
{
  public <- FALSE
  if (public == TRUE)
  {
  # Extract
  network <- sample$network
  n <- sample$n
  model <- sample$model
  max_number <- sample$max_number
  indicator <- sample$indicator
  ergm_theta <- sample$ergm_theta
  hergm_theta <- sample$hergm_theta
  sample_size <- nrow(indicator)
  
  # Observed network
  d <- network # Observed network
  n <- d$gal$n # Number of nodes
  # Observed values of statistics
  observed.degree <- summary(d ~ degree(1:n-1)) # Degree distribution
  observed.components <- component.dist(d)
  observed.component.number <- length(observed.components$csize) # Number of components
  observed.max.component.size <- max(observed.components$csize) # Size of largest component
  observed.distances <- geodist(d)
  observed.distances <- observed.distances$gdist
  observed.frequencies_distances <- table(observed.distances) # First column: frequency of self loops; columns 2:number: frequencies finite and (last column) infinite distances
  number <- length(observed.frequencies_distances) - 1 - (sum(observed.distances == Inf) > 0) # Number of distances minus 0-distance minus Inf-distance
  observed.distance.label <- rownames(observed.frequencies_distances)[2:(number+1)]
  observed.distance <- observed.frequencies_distances[2:(number+1)] # Frequencies of finite distances
  observed.star <- summary(d ~ kstar(2))
  observed.triangle <- summary(d ~ triangle)
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

  # Simulated networks
  edgelists <- simulate.hergm(  
      model = model,
      network = network,
      max_number = max_number, 
      indicator = indicator, 
      sample_size = sample_size,
      verbose = 1
      )
  for (i in 1:sample_size) 
    {
    sample <- edgelists$edgelist[[i]]
    output <- summary_sample_network(sample=sample, n=n, output, i) 
    }

  # Goodness-of-fit plots
  par(mfrow = c(2, 2))
  hist(output$component.number, 50, prob = T, xlim = c(-max(abs(output$component.number)), max(abs(output$component.number))), main = "", xlab = "number of components", ylab = "", cex.lab=1.5) 
  abline(v = c(quantile(output$component.number, 0.025), quantile(output$component.number, 0.975)), col = "blue")
  abline(v = observed.component.number, col="red")
  hist(output$max.component.size, 50, prob = T, xlim = c(-max(abs(output$max.component.size)), max(abs(output$max.component.size))), main = "", xlab = "max(component size)", ylab = "", cex.lab=1.5)
  abline(v = c(quantile(output$max.component.size, 0.025), quantile(output$max.component.size, 0.975)), col = "blue")
  abline(v = observed.max.component.size, col="red")
  boxplot(output$distance[,1:10], main = "", xlab = "geodesic distances", ylab = "", cex.lab=1.5)
  points(x = c(observed.distance[1:10]), col = "red", type = "b")
  boxplot(output$degree[,1:max(observed.degree)], ylim = c(0, max(output$degree)), main = "", xlab = "degrees", ylab = "", cex.lab=1.5)
  points(x = c(observed.degree[1:max(observed.degree)]), col = "red", type = "b")
  # hist(output$star, 50, prob = T, xlim = c(-max(abs(output$star)), max(abs(output$star))), main = "", xlab = "number of 2-stars", ylab = "")
  # abline(v = c(quantile(output$star, 0.025), quantile(output$star, 0.975)), col = "blue")
  # abline(v = observed.star, col="red")
  # hist(output$triangle, 50, prob = T, xlim = c(-max(abs(output$triangle)), max(abs(output$triangle))), main = "", xlab = "number of triangles", ylab = "", cex.lab=1.5)
  # abline(v = c(quantile(output$triangle, 0.025), quantile(output$triangle, 0.975)), col = "red")
  # abline(v = observed.triangle, col="red")

  output
  }
}

