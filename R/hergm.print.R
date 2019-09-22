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

print.hergm <- function (x, 
                         ...)
{
 if (x$simulate == FALSE)
   {
   max_number <- x$max_number
   if (x$relabel %in% c(1, 2)) hergm_theta <- x$relabeled.hergm_theta
   else hergm_theta <- x$hergm_theta
   cat("\n====================\n")
   cat("Summary of model fit\n")
   cat("====================\n")
   cat("\nFormula: ")
   formula <- deparse(x$formula)
   cat(formula, "\n")
   if (x$method == "bayes") # Bayes method
     {
     cat("\nSize of MCMC sample from posterior: ", x$sample_size)
     cat("\n")
     cat("\nPosterior quantiles                         2.5%        50%      97.5%")
     cat("\n----------------------------------------------------------------------")
     if (x$hyper_prior == 1)
       {
       if (!(is.null(x$alpha))) 
         {
         cat("\nConcentration parameter alpha:           ")
         cat(formatC(quantile(x$alpha, .025), format="f", width=7, digits=3), "   ")
         cat(formatC(quantile(x$alpha, .500), format="f", width=7, digits=3), "   ")
         cat(formatC(quantile(x$alpha, .975), format="f", width=7, digits=3))
         }
       if (!(is.null(x$eta_mean))) 
         {
         for (i in 1:ncol(x$eta_mean))
           {
           cat("\nMean of parameters of hergm term ", i, ":      ", sep="")
           cat(formatC(quantile(x$eta_mean[,i], .025), format="f", width=7, digits=3), "   ")
           cat(formatC(quantile(x$eta_mean[,i], .500), format="f", width=7, digits=3), "   ")
           cat(formatC(quantile(x$eta_mean[,i], .975), format="f", width=7, digits=3))
           }
         }
       if (!(is.null(x$eta_precision)))
         {
         for (i in 1:ncol(x$eta_precision))
           {
           cat("\nPrecision of parameters of hergm term ", i, ": ", sep="")
           cat(formatC(quantile(x$eta_precision[,i], .025), format="f", width=7, digits=3), "   ")
           cat(formatC(quantile(x$eta_precision[,i], .500), format="f", width=7, digits=3), "   ")
           cat(formatC(quantile(x$eta_precision[,i], .975), format="f", width=7, digits=3))
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
           if (max_number >= 2) 
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
     if (!(is.null(x$ergm_theta)))
       {
       for (i in 1:ncol(x$ergm_theta)) 
         {
         cat("\nergm term ", i, " parameter:                   ", sep="")
         cat(formatC(quantile(x$ergm_theta[,i], .025), format="f", width=7, digits=3), "   ")
         cat(formatC(quantile(x$ergm_theta[,i], .500), format="f", width=7, digits=3), "   ")
         cat(formatC(quantile(x$ergm_theta[,i], .975), format="f", width=7, digits=3))
         }
       cat("\n----------------------------------------------------------------------")
       } 
     cat("\n\n")
    }
  else # ML method
    {
    summary(x$mlergm_out)
    #if (x$estimate_parameters == TRUE)
    #  {
    # if (x$parameterization == "size") cat("\nSize-dependent parameterization of the form theta * log(size of block).\n")
    #  else if (x$parameterization == "offset") cat("\nSize-dependent parameterization of the form theta - log(size of block).\n")
    #  else cat("\nSize-independent parameterization: all within-block parameters are constant across blocks.\n")
    #  cat("\nEstimates of theta (S.E)                           ")
    #  cat("\n------------------------------------------")
    #  for (i in 1:length(x$results$parameters))
    #    {
    #    cat("\nwithin-block parameter ", i, ":   ", sep="")
    #    cat(formatC(x$results$parameters[[i]], format="f", width=6, digits=3), " ")
    #    if (x$results$st.error[[i]] < .001) cat("<0.001")
    #    else cat(formatC(x$results$st.error[[i]], format="f", width=6, digits=3), " ")
    #    }
    #  if (max_number >= 2) 
    #    {
    #    cat("\nbetween-block parameter:    ", sep="")
    #    cat(formatC(x$results$between_parameter, format="f", width=6, digits=3), " ")
    #    if (x$results$st.error.between < .001) cat("<0.001")
    #    else cat(formatC(x$results$st.error.between, format="f", width=6, digits=3), " ")
    #    }
    #  cat("\n------------------------------------------")
    #  cat("\n\n")
    #  if ((max_number >= 2) && ((x$parameterization == "size") || (x$parameterization == "offset"))) 
    #    {
    #    cat("Offset log(size of block) of block 1, 2, ...:")
    #    for (i in 1:max_number) cat(" ", log(sum(x$results$partition==i)))
    #    cat("\n\n")
    #    }
    #  }
    #} else { cat("\nThe parameters of the model have not been estimated.\n\n")
    }
  }
}

