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

mcmc.diagnostics.hergm <- function(object, ...)
{
  if (object$method == "ml") 
    {
    cat("\nConvergence checks have been implemented for method = bayes but not method = ml.")
    cat("\nPlease use gof() to assess the goodness-of-fit of the model.\n\n")
    }
  else
    {
    sample_size <- object$sample_size
    hyper_prior <- object$hyper_prior
    samle_size <- object$sample_size
    network <- object$network
    n <- network$gal$n
    verbose <- object$verbose
    if (verbose >= 0) cat("\nConvergence check using R function mcgibbsit()...")

    output <- list()
    warn <- FALSE
    if (!(is.null(object$hergm_theta)))
      {
      nonfixed <- NULL
      if (object$relabel %in% c(1, 2)) hergm_theta <- object$relabeled.hergm_theta
      else hergm_theta <- object$hergm_theta
      for (i in 1:ncol(hergm_theta))
        {
        if (sd(hergm_theta[,i]) > 0)
          {
          nonfixed <- cbind(nonfixed, i)
          }
        }
      }
    if (sample_size < 600) 
      {
      warning("mcgibbsit() cannot be used with a sample size smaller than 600, therefore no convergence checks have been run.")
      }
    else
      {
      if (!(is.null(object$ergm_theta))) 
        { 
        output$mcmc.ergm_theta <- mcgibbsit(object$ergm_theta)
        if (max(output$mcmc.ergm_theta$resmatrix[,3]) > sample_size) warn <- TRUE
        if (max(output$mcmc.ergm_theta$resmatrix[,5]) > 10) warn <- TRUE
        }
      if (!(is.null(hergm_theta)))
      #if (((object$relabel %in% c(1, 2)) && (!(is.null(object$relabeled.hergm_theta)))) || ((object$max_number == 1) && ((!(is.null(object$relabeled.hergm_theta)))))) 
        {
        output$mcmc.hergm_theta <- mcgibbsit(hergm_theta[,nonfixed])
        if ((max(output$mcmc.hergm_theta$resmatrix[,3])) > sample_size) warn <- TRUE
        if ((max(output$mcmc.hergm_theta$resmatrix[,5])) > 10) warn <- TRUE
        }
      if (hyper_prior == 1)
        {
        output$mcmc.alpha <- mcgibbsit(object$alpha)
        if ((max(output$mcmc.alpha$resmatrix[,3])) > sample_size) warn <- TRUE
        if ((max(output$mcmc.alpha$resmatrix[,5])) > 10) warn <- TRUE
        output$mcmc.eta_mean <- mcgibbsit(object$eta_mean)
        if ((max(output$mcmc.eta_mean$resmatrix[,3])) > sample_size) warn <- TRUE
        if ((max(output$mcmc.eta_mean$resmatrix[,5])) > 10) warn <- TRUE
        output$mcmc.eta_precision <- mcgibbsit(object$eta_precision)
        if ((max(output$mcmc.eta_precision$resmatrix[,3])) > sample_size) warn <- TRUE
        if ((max(output$mcmc.eta_precision$resmatrix[,5])) > 10) warn <- TRUE
        }
      }
    count <- 0 # Count the number of trace plots
    if (!(is.null(object$ergm_theta))) count <- count + 1
    if (!(is.null(hergm_theta))) count <- count + 1
    if (hyper_prior == 1) count <- count + 3
    oldpar <- par(no.readonly = TRUE)    
    on.exit(par(oldpar))           
    par(mfrow=c(ceiling(count/2), 2))
    if (!(is.null(object$ergm_theta))) 
      {
      output$ergm_theta <- matplot(object$ergm_theta, type="l", xlab="ergm term parameters", ylab="", main="", cex.lab=1.25)
      }
    if (!(is.null(hergm_theta)))
      {
      output$hergm_theta <- matplot(hergm_theta[,nonfixed], type="l", xlab="hergm term parameters", ylab="", main="", cex.lab=1.25) 
      }
    if (hyper_prior == 1)
      {
      output$alpha <- plot(object$alpha, type="l", xlab="concentration parameter alpha", ylab="", main="", cex.lab=1.25)
      output$eta_mean <- matplot(object$eta_mean, type="l", xlab="means of hergm term parameters", ylab="", main="", cex.lab=1.25)
      output$eta_precision <- matplot(object$eta_precision, type="l", xlab="precisions of hergm term parameters", ylab="", main="", cex.lab=1.25)
      }
    if (verbose >= 0)
      {
      if (warn == TRUE) 
        {
        warning("There are signs of non-convergence: to view details, enter\n         \'print(object$mcmc.diagnostics)\'\n         where object is the object returned by function hergm().\n\n")
        }
      else cat("OK\n\n")
      if ((object$max_number > 1) && (object$number_fixed < object$network$gal$n) && (!(object$relabel %in% c(1, 2)))) warning("The label-switching problem has not been solved so that the interpretation of block-dependent hergm terms parameters is problematic.")
      }
  
    output
    }
}
