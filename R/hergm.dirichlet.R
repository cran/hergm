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

hergm.dirichlet <- function(n, number, alpha, eta_mean, eta_sd) 
# Note: works in one-dimensional case
{
  indicator <- vector(length=n)
  eta <- vector(length=number)
  output <- .C("Dirichlet",
                     as.integer(n),
                     as.integer(number),
                     as.numeric(alpha),
                     as.numeric(eta_mean),
                     as.numeric(eta_sd),
                     indicator = as.integer(indicator),
                     eta = as.numeric(eta),
                     PACKAGE="hergm")
  output
}

