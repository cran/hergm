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

summary <- function(sample = NULL,
                    ...)
 {
 UseMethod("summary")
 }

summary.hergm <- function(sample = NULL,
                    ...)
# input: output of function hergm
# output: summary of output of function hergm
{
  if (sample$simulate == TRUE) # MCMC sample of networks from model
    {
    # Initialize
    sample_size <- sample$sample_size
    n <- sample$network$gal$n
    output <- list()
    output$component.number <- vector(length = sample_size)
    output$max.component.size <- vector(length = sample_size)
    output$distance.label <- matrix(0, nrow = sample_size, ncol = n)
    output$distance <- matrix(0, nrow = sample_size, ncol = n)
    output$edges <- vector(length = sample_size)
    output$degree <- matrix(0, nrow = sample_size, ncol = n)
    output$stars <- vector(length = sample_size)
    output$triangles <- vector(length = sample_size)
    for (i in 1:sample$sample_size)
      {
      summary_sample_network(sample, n, output, i)
      }
    }
  else # MCMC sample of clustering and parameters from posterior
    {
    cat("\nSummary of posterior of clustering of nodes:")
    sample$p_i_k
    output <- hergm.plot(sample)
    cat("\nSummary of marginal posteriors of parameters can be obtained by using running \"hergm$parameter\", where \"parameter\" is replaced by the parameter of interest: see \"?hergm\".")
    }
  output
}

