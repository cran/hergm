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

plot.hergm <- function(x,
                       threshold = c(.7, .8, .9),
                       ...)
{
  par(mfrow=c(1,1))
  formula <- x$formula
  max_number <- x$max_number
  network <- hergm.getnetwork(formula, max_number)
  n <- network$gal$n
  if (x$method == "ml") # ML method
    {
    plot_neighborhoods(network, x$results$partition)
    }
  else # Bayes methods  
    {
    if (max_number > 1)
      { 
      if (is.directed(network)) gmode <- "digraph" 
      else gmode <- "graph"
      if ((x$relabel %in% c(1, 2)) && (!(is.null(x$p_i_k)))) # relabel = 1
        {
        p <- gplot(network, gmode=gmode, mode="fruchtermanreingold", vertex.cex=1, vertex.col=0, vertex.border=0, displaylabels=TRUE, label=c(1:n), label.cex=0.8)
        for(i in 1:nrow(p))
          {
          ergmm.drawpie(center=p[i,], radius=0.25, probs=x$p_i_k[i,])
          }
        }
      else # Choose relabel = 3 instead, which is always possible, though the computing time is quadratic in number of nodes
        { 
        # Extract
        indicator <- x$indicator
        sample_size <- x$sample_size
        # Estimate same-block-membership posterior probabilities
        if (sample_size > 0)
          {
          p <- matrix(0, n, n)
          for (i in 1:n)
            {
            for (j in 1:n)
              {
              for (h in 1:sample_size)
                {
                if (indicator[h,i] == indicator[h,j]) p[i,j] <- p[i,j] + 1
                }
              if (i == j) p[i,j] <- 0
              else p[i,j] <- p[i,j] / sample_size
              }
           }
          # Plot
          #plot_weighted_graph(network, p) 
          threshold <- unique(threshold)
          threshold <- sort(threshold)
          count <- length(threshold) + 1
          par(mfrow=c(ceiling(count/2), 2))
          coordinates <- gplot(network, gmode=gmode, mode="fruchtermanreingold", vertex.cex=1, vertex.col="black", vertex.border=0, displaylabels=TRUE, label=c(1:n), label.cex=0.8, edge.col="black", xlab="observed network", cex.lab=1.25)
          for (i in 1:length(threshold))
            {
            same.neighborhood.graph <- (p >= threshold[i])
            same.neighborhood.graph <- as.network(same.neighborhood.graph, directed=is.directed(network))
            gplot(same.neighborhood.graph, coord=coordinates, gmode="graph", vertex.col = "black", vertex.cex=1, vertex.border=0, displaylabels=TRUE, label.cex=0.8, edge.col="gray", xlab=paste("same-neighborhood-network (", threshold[i], ")", sep=""), cex.lab=1.25) 
            }
          }
        }
      }
  }
}

