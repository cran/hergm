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

simulate <- function(object, 
                     max_number = NULL,
                     indicator = NULL,
                     eta = NULL,
                     sample_size = 1,
                     verbose = 0, 
                     ...)
 {
 UseMethod("simulate")
 }

simulate.hergm <- function(object,
                           max_number = NULL, 
                           indicator = NULL,
                           eta = NULL,
                           sample_size = 1,
                           verbose = 0, 
                           ...)
{
  # Extract
  condition <- (is(object) == "hergm")
  condition <- min(condition)
  if (condition == TRUE)
    { 
    formula <- object$formula
    max_number <- object$max_number
    indicator <- object$indicator
    ergm_theta <- object$ergm_theta
    hergm_theta <- object$hergm_theta
    eta <- cbind(ergm_theta, hergm_theta)
    if (is.null(sample_size)) sample_size <- object$sample_size
    verbose <- object$verbose
    }
  else formula <- object
  if (sample_size == 1)
    {
    if (is.vector(indicator)) indicator <- t(as.matrix(indicator))
    if (is.vector(indicator)) eta <- t(as.matrix(eta))
    }
  else 
    {
    indicator_matrix <- NULL
    eta_matrix <- NULL
    if (is.vector(indicator)) 
      {
      for (i in 1:sample_size) indicator_matrix <- rbind(indicator_matrix, indicator)
      indicator <- indicator_matrix
      }
    if (is.vector(eta)) 
      {
      for (i in 1:sample_size) eta_matrix <- rbind(eta_matrix, eta)
      eta <- eta_matrix
      }
    }
  n <- ncol(indicator)

  # Sample
  edgelists <- list()
  edgelists$edgelist <- list()
  for (i in 1:sample_size)
    { 
    if (verbose == 1) 
      {
      cat("\nSample", i)
      cat("\n------------------------------------------------------------------")
      cat("\nInput:")
      cat("\n- parameters:", eta[i,]) 
      if (length(indicator) > 0) 
        {
        cat("\n- block memberships:", indicator[i,])
        }
      cat("\n")
      }
    object.hergm <- hergm(formula, 
      max_number = max_number, 
      eta = eta[i,], 
      indicator = indicator[i,], 
      simulate = TRUE, 
      sample_size = 1,
      verbose = verbose
      )
    edgelists$edgelist[[i]] <- object.hergm$edgelist
    }
  
  edgelists
}

