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

print <- function(object,
                  ...)
 {
 UseMethod("print")
 }

print.hergm <- function (object, 
                         ...)
{
 if (object$simulate == FALSE) 
   {
   max_number <- object$max_number
   if (object$relabel %in% c(1, 2)) hergm_theta <- object$relabeled.hergm_theta
   else hergm_theta <- object$hergm_theta
   cat("\n====================\n")
   cat("Summary of model fit\n")
   cat("====================\n")
   cat("\nFormula: ")
   formula <- deparse(object$formula)
   cat(formula, "\n")
   cat("\nSize of MCMC sample from posterior: ", object$sample_size)
   cat("\n")
   cat("\nPosterior quantiles                         2.5%        50%      97.5%")
   cat("\n----------------------------------------------------------------------")
   if (object$hyper_prior == 1)
     {
     if (!(is.null(object$alpha))) 
       {
       cat("\nConcentration parameter alpha:           ")
       cat(formatC(quantile(object$alpha, .025), format="f", width=7, digits=3), "   ")
       cat(formatC(quantile(object$alpha, .500), format="f", width=7, digits=3), "   ")
       cat(formatC(quantile(object$alpha, .975), format="f", width=7, digits=3))
       }
     if (!(is.null(object$eta_mean))) 
       {
       for (i in 1:ncol(object$eta_mean))
         {
         cat("\nMean of parameters of hergm term ", i, ":      ", sep="")
         cat(formatC(quantile(object$eta_mean[,i], .025), format="f", width=7, digits=3), "   ")
         cat(formatC(quantile(object$eta_mean[,i], .500), format="f", width=7, digits=3), "   ")
         cat(formatC(quantile(object$eta_mean[,i], .975), format="f", width=7, digits=3))
         }
       }
     if (!(is.null(object$eta_precision)))
       {
       for (i in 1:ncol(object$eta_precision))
         {
         cat("\nPrecision of parameters of hergm term ", i, ": ", sep="")
         cat(formatC(quantile(object$eta_precision[,i], .025), format="f", width=7, digits=3), "   ")
         cat(formatC(quantile(object$eta_precision[,i], .500), format="f", width=7, digits=3), "   ")
         cat(formatC(quantile(object$eta_precision[,i], .975), format="f", width=7, digits=3))
         }
       cat("\n----------------------------------------------------------------------")
       }
     }
   if (!(is.null(hergm_theta)))
     {
     for (i in 1:ncol(hergm_theta))
       {
       term <- ceiling(i / (max_number + 1))
       k <- i - (term - 1) * (max_number + 1)
       if (k <= max_number) # Within-block parameter 
         {
         cat("\nhergm term ", term, ": parameter of block ", k, ":      ", sep="")
         cat(formatC(quantile(hergm_theta[,i], .025), format="f", width=7, digits=3), "   ")
         cat(formatC(quantile(hergm_theta[,i], .500), format="f", width=7, digits=3), "   ")
         cat(formatC(quantile(hergm_theta[,i], .975), format="f", width=7, digits=3))
         }
       else # Between-block parameter
         {
         if (max_number > 2) 
           {
           cat("\nhergm term ", term, ": between-block parameter:   ", sep="")
           cat(formatC(quantile(hergm_theta[,i], .025), format="f", width=7, digits=3), "   ")
           cat(formatC(quantile(hergm_theta[,i], .500), format="f", width=7, digits=3), "   ")
           cat(formatC(quantile(hergm_theta[,i], .975), format="f", width=7, digits=3))
           }
         } 
       }
     cat("\n----------------------------------------------------------------------")
     }
   if (!(is.null(object$ergm_theta)))
     {
     for (i in 1:ncol(object$ergm_theta)) 
       {
       cat("\nergm term ", i, " parameter:                   ", sep="")
       cat(formatC(quantile(object$ergm_theta[,i], .025), format="f", width=7, digits=3), "   ")
       cat(formatC(quantile(object$ergm_theta[,i], .500), format="f", width=7, digits=3), "   ")
       cat(formatC(quantile(object$ergm_theta[,i], .975), format="f", width=7, digits=3))
       }
     cat("\n----------------------------------------------------------------------")
     } 
   cat("\n\n")
  }
}

