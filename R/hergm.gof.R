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

gof <- function(object,
                ...)
 {
 UseMethod("gof")
 }

gof.hergm <- function(object, 
                      ...)
{
  # Extract
  network <- object$network
  d <- network 
  directed <- is.directed(d)
  n <- d$gal$n
  indicator <- object$indicator
  sample_size <- nrow(indicator)
  verbose <- object$verbose
 
  # Observed network
  observed.components <- silent(component.dist(d))
  observed.component.number <- length(observed.components$csize) # Number of components
  observed.max.component.size <- max(observed.components$csize) # Size of largest component
  observed.distances <- geodist(d)
  observed.distances <- observed.distances$gdist
  observed.frequencies_distances <- table(observed.distances) # First column: frequency of self loops; columns 2:number: frequencies finite and (last column) infinite distances
  number <- length(observed.frequencies_distances) - 1 - (sum(observed.distances == Inf) > 0) # Number of distances minus 0-distance minus Inf-distance
  observed.distance.label <- rownames(observed.frequencies_distances)[2:(number+1)]
  observed.distance <- observed.frequencies_distances[2:(number+1)] # Frequencies of finite distances
  observed.edges <- summary(d ~ edges)
  if (directed == TRUE)
    {
    observed.degree <- summary(d ~ odegree(1:n-1)) # Degree distribution
    observed.star <- summary(d ~ ostar(2))
    observed.triangle <- summary(d ~ ttriple)
    }
  else
    {
    observed.degree <- summary(d ~ degree(1:n-1)) # Degree distribution
    observed.star <- summary(d ~ kstar(2))
    observed.triangle <- summary(d ~ triangle)
    }

  # Simulated networks
  object.hergm <- simulate.hergm(object, verbose=verbose)
  output <- summary_sample_network(edgelists=object.hergm$edgelists, sample_size=sample_size, directed=directed, n=n)

  # Goodness-of-fit plots
  par(mfrow = c(2, 3))
  #hist(output$component.number, 50, prob = T, xlim = c(0, max(abs(output$component.number))), main = "", xlab = "number of components", ylab = "", cex.lab=1.25) 
  #abline(v = c(quantile(output$component.number, 0.025), quantile(output$component.number, 0.975)), col = "blue")
  #abline(v = observed.component.number, col="red")
  hist(output$max.component.size, 50, prob = T, xlim = c(0, max(abs(output$max.component.size))), main = "", xlab = "size of largest component", ylab = "", cex.lab=1.25)
  #abline(v = c(quantile(output$max.component.size, 0.025), quantile(output$max.component.size, 0.975)), col = "blue")
  abline(v = observed.max.component.size, col="red")
  boxplot(output$distance[,1:(n-2)], main = "", xlab = "geodesic distances", ylab = "", cex.lab=1.25)
  points(x = c(observed.distance[1:(n-2)]), col = "red", type = "b")
  boxplot(output$degree[,1:(n-1)], ylim = c(0,n-1), main = "", xlab = "degrees", ylab = "", cex.lab=1.25)
  points(x = c(observed.degree), col = "red", type = "b")
  hist(output$edges, 50, prob = T, main = "", xlab = "number of edges", ylab = "", cex.lab=1.25)
  #abline(v = c(quantile(output$star, 0.025), quantile(output$star, 0.975)), col = "blue")
  abline(v = observed.edges, col="red")
  if (directed == TRUE) xlab <- "number of 2-out-stars"
  else xlab <- "number of 2-stars"
  hist(output$star, 50, prob = T, main = "", xlab = xlab, ylab = "", cex.lab=1.25)
  #abline(v = c(quantile(output$star, 0.025), quantile(output$star, 0.975)), col = "blue")
  abline(v = observed.star, col="red")
  if (directed == TRUE) xlab <- "number of transitive triples"
  else xlab <- "number of triangles"
  hist(output$triangle, 50, prob = T, main = "", xlab = xlab, ylab = "", cex.lab=1.25)
  #abline(v = c(quantile(output$triangle, 0.025), quantile(output$triangle, 0.975)), col = "blue")
  abline(v = observed.triangle, col="red")

  output
}

