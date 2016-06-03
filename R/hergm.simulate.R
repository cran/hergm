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

simulate <- function(model = NULL,
                     network = NULL,
                     eta = NULL,
                     max_number = NULL,
                     indicator = NULL,
                     sample_size = 1000,
                     sample = NULL,
                     verbose = 1, 
                     ...)
 {
 UseMethod("simulate")
 }

simulate.hergm <- function(model = NULL,
                           network = NULL,
                           eta = NULL,
                           max_number = NULL, 
                           indicator = NULL,
                           sample_size = 1000, 
                           sample = NULL,
                           verbose = 1, 
                           ...)
# input: postprocessed sample, number of nodes, number of blocks, number covariates, observed network
# output: goodness of fit data
{
  # Extract
  if (!is.null(sample))
    { 
    network <- sample$network
    n <- sample$network$gal$n
    model <- sample$model
    max_number <- sample$max_number
    indicator <- sample$indicator
    ergm_theta <- sample$ergm_theta
    hergm_theta <- sample$hergm_theta
    sample_size <- sample$sample_size
    }
  else 
    {
    n <- network$gal$n
    }

  # Initialize
  output <- list()
  output$component.number <- vector(length = sample_size) 
  output$max.component.size <- vector(length = sample_size)
  output$distance.label <- matrix(0, nrow = sample_size, ncol = n)
  output$distance <- matrix(0, nrow = sample_size, ncol = n)
  output$edges <- vector(length = sample_size)
  output$degree <- matrix(0, nrow = sample_size, ncol = n)
  output$stars <- vector(length = sample_size)
  output$triangles <- vector(length = sample_size)

  # Sample

print(indicator)

  edgelists <- list()
  edgelists$edgelist <- list()
  for (i in 1:sample_size)
    { 
    eta <- c(ergm_theta[i,], hergm_theta[i,]) 
    verbose <- 1
    if (verbose == 1) 
      {
      cat("\nSample", i)
      cat("\n------------------------------------------------------------------")
      cat("\nInput:")
      if (length(ergm_theta) > 0) cat("\n- parameters:", ergm_theta[i,])
      if (length(hergm_theta) > 0) 
        {
        cat("\n- block parameters:", hergm_theta[i,])
        cat("\n- block memberships:", indicator[i,])
        }
      cat("\n")
      }
    edgelist <- hergm(model$formula, 
      max_number = max_number, 
      eta = eta, 
      indicator = indicator[i,], 
      simulate = TRUE, 
      samplesize = 1,
      verbose = verbose
      )
    edgelists$edgelist[[i]] <- edgelist
    }

  edgelists
}

