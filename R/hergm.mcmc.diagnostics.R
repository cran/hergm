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

mcmc.diagnostics <- function(sample = NULL,
                  ...)
{
  hyper_prior <- sample$hyper_prior
  output <- list()
  if (hyper_prior == 1)
    {
    output$alpha <- mcgibbsit(sample$alpha)
    output$eta_mean <- mcgibbsit(sample$eta_mean)
    output$eta_precision <- mcgibbsit(sample$eta_precision)
    }
  output$hergm_theta <- mcgibbsit(sample$hergm_theta[,1:ncol(sample$hergm_theta)])
  if (hyper_prior == 1)
    {
    par(mfrow=c(4,1))
    plot(sample$alpha, type="l", xlab=expression(alpha), ylab="", main="", cex.lab=1.5)
    matplot(sample$eta_mean, type="l", xlab=expression(mu), ylab="", main="", cex.lab=1.5)
    matplot(sample$eta_precision, type="l", xlab=expression(Sigma^{-1}), ylab="", main="", cex.lab=1.5)
    }
  else matplot(sample$hergm_theta, type="l", xlab=expression(theta), ylab="", main="", cex.lab=1.5)
  output
}

