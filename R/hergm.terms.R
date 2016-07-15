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

# Upon encountering a model term such as [name](args), ergm and related
# routines will call a function of the form InitErgmTerm.[name].  The specific
# function call will be of the form
#   InitErgmTerm.[name](network, model, arguments, ...)
# where network is the network object, model is a model object that should be
# updated and then returned, arguments is the list (if any) of arguments
# passed to the model term by the user, and ... includes any arguments
# passed to the InitErgmTerm function from within the program.
#
# Arguments of
# the latter type include such items as drop (a logical flag telling whether
# degenerate terms should be dropped) and expanded (a logical flag used
# by curved exponential family terms).  Because such arguments are usually
# passed to ALL InitErgmTerm functions, regardless of whether they are used,
# it is important that each InitErgmTerm function declaration include the
# dot-dot-dot (...) argument.  Finally, such arguments are not guaranteed
# to be passed when the InitErgmTerm function is called, so any InitErgmTerm function
# requiring such an argument should supply a default value.
#
# An example:  If drop=TRUE is passed from inside the program,
# then the statement
#     ergm(network ~ triangle + kstar (2:4) + nodematch("sex"))
# results in the following function calls:
#     InitErgmTerm.triangle (network, model, list(), drop=TRUE)
#     InitErgmTerm.kstar (network, model, list(2:4), drop=TRUE)
#     InitErgmTerm.nodematch (network, model, list("sex"), drop=TRUE)
#
# Each InitErgmTerm.[name] function should check its argument list for errors, 
# then set termnumber to 1+length(model$terms).
# Next, it should add the names of the statistics that
# will be computed to the vector model$coef.names.  These names must be
# concatenated onto model$coef.names in the same order they will be produced
# by the changestat function.
# Finally, it should create 
# model$terms[[termnumber]] , a list with the following elements, some
# required and some optional:
#                                                                                                    
# Required arguments of model$terms[[termnumber]]
# -----------------------------------------------
#    name: This is the (text) name of the term.  It is expected that there
#          is a C function called d_[name].
#  soname: This is the (text) name of the package containing the C function
#          called d_[name].
#  inputs: This is a (numeric) vector with at least 3 elements, as described
#          below:
#    Element 1 -- For functions that require a vector of covariates, either
#                 nodal or dyadic, this optional value is the number of
#                 input parameters BEFORE the beginning of the covariate
#                 vector.  For instance, if there are no input parameters
#                 passed before the covariate vector, this value should be
#                 set to zero.  The changestat function in C will be passed a
#                 pointer to the start of this vector of covariates, though
#                 the changestat function may choose to ignore this pointer,
#                 in which case the value of element 1 is arbitrary.
#    Element 2 -- The number of change statistics returned by the function.
#    Element 3 -- The total number of input parameters and covariates
#                 to be passed to the function.  If there are no nodal or 
#                 dyadic covariates, the value of element 1 is arbitrary.
#   Element 4+ -- The input parameters to be passed to the function.
#                 For example, if element 3 equals 3, then elements
#                 4, 5, 6 are the parameters to be passed.  No 4th element
#                 is necessary if element 3==0.  If there are nodal or
#                 dyadic covariates, they should be appended after any other
#                 input parameters (and element 1 may then be set to the
#                 number of other input parameters excluding the covariates).
#
# Optional arguments of model$terms[[termnumber]]
# -----------------------------------------------
#    dependence: Logical variable telling whether addition of this term to
#                the model makes the model into a dyadic dependence model.
#                If none of the terms sets dependence==TRUE, then the model
#                is assumed to be a dyadic independence model, which means
#                that the pseudolikelihood estimate coincides with the
#                maximum likelihood estimate.  Default value:  TRUE
#        params: For curved exponential family models, this argument must be
#                a list:  Each item in the list should be named with the
#                corresponding parameter name (one or more of these will
#                probably coincide with the coef.names used when
#                initialfit=TRUE; the initial values of such parameter values
#                will be set by MPLE, so their values in params are ignored.)
#                Any parameter not having its initial value set by MPLE
#                should be given its initial value in this params list.
#           eta: A function that gives the map from theta (the canonical
#                parameters associated with the statistics for this term)
#                to eta (the corresponding curved parameters).  The length
#                of eta is the same as the length of the params list above.
#                This function takes two args:  theta and length(eta).
#      gradient: A function that gives the gradient of the eta map above.
#                If theta has length p and eta has length q, then gradient
#                should return a p by q matrix.
#                This function takes two args:  theta and length(eta).
#  emptynetworkstats: Vector of values (if nonzero) for the statistics evaluated
#                on the empty network.  If all are zero for this term, this
#                argument may be omitted.  Example:  If the degree0 term is
#                among the statistics, this argument is necessary because
#                degree0 = number of nodes for the empty network.


######################################################### 
InitErgmTerm.edges_i <- function(nw, arglist, ...) # Michael 
{
  n <- nw$gal$n # Number of nodes
  #print("InitErgmTerm.edges_i")
  #print(nw$gal$n)
  #print(nw$terms)
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=FALSE, 
    varnames = c("number", "indicator", "theta"),
    vartypes = c("numeric", "numeric", "numeric"),
    defaultvalues = list(NULL, NULL, NULL),
    required = c(FALSE, FALSE, FALSE)) 
  termnumber <- 1 + length(nw$terms)
  #print("InitErgmTerm.edges_i")
  number <- nw$max_number # (Maximum) number of categories
  #print("number")
  #print(number)
  if (is.null(a$indicator)) 
    {
    indicator <- vector(mode = "numeric", length = n) # Category indicators  
    for (i in 1:length(indicator)) indicator[i] <- 1
    }
  else indicator <- a$indicator
  #print(indicator)
  if (is.null(a$theta)) 
    {
    theta <- vector(mode = "numeric", length = number + 1) # Within- and between-category parameters
    for (i in 1:length(theta)) theta[i] <- 0
    }
  else theta <- a$theta
  #print(theta)
  list(name = "edges_i", coef.names="edges_i", 
                                inputs = c(0, 1, 1+length(indicator)+length(theta), c(number, indicator, theta)),
                                dependence = FALSE)
  #print(nw$terms[[termnumber]])
}

######################################################### 
InitErgmTerm.arcs_i <- function(nw, arglist, ...) # Michael 
{
  n <- nw$gal$n # Number of nodes
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=FALSE,
    varnames = c("number", "indicator", "theta"),
    vartypes = c("numeric", "numeric", "numeric"),
    defaultvalues = list(NULL, NULL, NULL),
    required = c(FALSE, FALSE, FALSE))
  termnumber <- 1 + length(nw$terms)
  #print("InitErgmTerm.arcs_i")
  indicator <- vector(mode = "numeric", length = n) # Category indicators  
  for (i in 1:length(indicator)) indicator[i] <- 1
  #print(indicator)
  number <- nw$max_number # (Maximum) number of categories
  #print(number)
  theta <- vector(mode = "numeric", length = number + 1) # Within- and between-category parameters
  for (i in 1:length(theta)) theta[i] <- 0 
  #print(theta)
  list(name = "arcs_i", coef.names="arcs_i", 
                                inputs = c(0, 1, 1+length(indicator)+length(theta), c(number, indicator, theta)),
                                dependence = FALSE)
  #print(nw$terms[[termnumber]])
}

######################################################### 
InitErgmTerm.arcs_j <- function(nw, arglist, ...) # Michael 
{
  n <- nw$gal$n # Number of nodes
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=FALSE,
    varnames = c("number", "indicator", "theta"),
    vartypes = c("numeric", "numeric", "numeric"),
    defaultvalues = list(NULL, NULL, NULL),
    required = c(FALSE, FALSE, FALSE))
  termnumber <- 1 + length(nw$terms)
  #print("InitErgmTerm.arcs_j")
  indicator <- vector(mode = "numeric", length = n) # Category indicators  
  for (i in 1:length(indicator)) indicator[i] <- 1
  #print(indicator)
  number <- nw$max_number # (Maximum) number of categories
  #print(number)
  theta <- vector(mode = "numeric", length = number + 1) # Within- and between-category parameters
  for (i in 1:length(theta)) theta[i] <- 0 
  #print(theta)
  list(name = "arcs_j", coef.names="arcs_j", 
                                inputs = c(0, 1, 1+length(indicator)+length(theta), c(number, indicator, theta)),
                                dependence = FALSE)
  #print(nw$terms[[termnumber]])
}

######################################################### 
InitErgmTerm.edges_ij <- function(nw, arglist, ...) # Michael 
{
  n <- nw$gal$n # Number of nodes
  a <- check.ErgmTerm(nw, arglist, bipartite=FALSE,
    varnames = c("number", "indicator", "theta"),
    vartypes = c("numeric", "numeric", "numeric"),
    defaultvalues = list(NULL, NULL, NULL),
    required = c(FALSE, FALSE, FALSE))
  termnumber <- 1 + length(nw$terms)
  #print("InitErgmTerm.edges_ij")
  indicator <- vector(mode = "numeric", length = n) # Category indicators  
  for (i in 1:length(indicator)) indicator[i] <- 1
  #print(indicator)
  number <- nw$max_number # (Maximum) number of categories
  #print(number)
  theta <- vector(mode = "numeric", length = number + 1) # Within- and between-category parameters
  for (i in 1:length(theta)) theta[i] <- 0
  #print(theta)
  list(name = "edges_ij", coef.names="edges_ij",
                                inputs = c(0, 1, 1+length(indicator)+length(theta), c(number, indicator, theta)),
                                dependence = FALSE)
  #print(nw$terms[[termnumber]])
}

######################################################### 
InitErgmTerm.mutual_i <- function(nw, arglist, ...) # Michael 
{
  n <- nw$gal$n # Number of nodes
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=FALSE,
    varnames = c("number", "indicator", "theta"),
    vartypes = c("numeric", "numeric", "numeric"),
    defaultvalues = list(NULL, NULL, NULL),
    required = c(FALSE, FALSE, FALSE))
  termnumber <- 1 + length(nw$terms)
  #print("InitErgmTerm.mutual_i")
  indicator <- vector(mode = "numeric", length = n) # Category indicators  
  for (i in 1:length(indicator)) indicator[i] <- 1
  #print(indicator)
  number <- nw$max_number # (Maximum) number of categories
  #print(number)
  theta <- vector(mode = "numeric", length = number + 1) # Within- and between-category parameters
  for (i in 1:length(theta)) theta[i] <- 0 
  #print(theta)
  list(name = "mutual_i", coef.names="mutual_i", 
                                inputs = c(0, 1, 1+length(indicator)+length(theta), c(number, indicator, theta)),
                                dependence = FALSE) 
  #print(nw$terms[[termnumber]])
}

######################################################### 
InitErgmTerm.mutual_ij <- function(nw, arglist, ...) # Michael 
{
  n <- nw$gal$n # Number of nodes
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=FALSE,
    varnames = c("number", "indicator", "theta"),
    vartypes = c("numeric", "numeric", "numeric"),
    defaultvalues = list(NULL, NULL, NULL),
    required = c(FALSE, FALSE, FALSE))
  termnumber <- 1 + length(nw$terms)
  #print("InitErgmTerm.mutual_ij")
  indicator <- vector(mode = "numeric", length = n) # Category indicators  
  for (i in 1:length(indicator)) indicator[i] <- 1
  #print(indicator)
  number <- nw$max_number # (Maximum) number of categories
  #print(number)
  theta <- vector(mode = "numeric", length = number + 1) # Within- and between-category parameters
  for (i in 1:length(theta)) theta[i] <- 0 
  #print(theta)
  list(name = "mutual_ij", coef.names="mutual_ij", 
                                inputs = c(0, 1, 1+length(indicator)+length(theta), c(number, indicator, theta)),
                                dependence = FALSE) 
  #print(nw$terms[[termnumber]])
}

######################################################### 
InitErgmTerm.triangle_ijk <- function(nw, arglist, ...) # Michael 
{
  n <- nw$gal$n # Number of nodes
  #print("InitErgmTerm.triangle_ijk")
  a <- check.ErgmTerm(nw, arglist, bipartite=FALSE,
    varnames = c("number", "indicator", "theta"),
    vartypes = c("numeric", "numeric", "numeric"),
    defaultvalues = list(NULL, NULL, NULL),
    required = c(FALSE, FALSE, FALSE))
  termnumber <- 1 + length(nw$terms)
  #print("InitErgmTerm.triangle_ijk")
  indicator <- vector(mode = "numeric", length = n) # Category indicators  
  for (i in 1:length(indicator)) indicator[i] <- 1
  #print(indicator)
  number <- nw$max_number # (Maximum) number of categories
  #print(number)
  theta <- vector(mode = "numeric", length = number + 1) # Within- and between-category parameters
  for (i in 1:length(theta)) theta[i] <- 0 
  #print(theta)
  list(name = "triangle_ijk", coef.names="triangle_ijk", 
                                inputs = c(0, 1, 1+length(indicator)+length(theta), c(number, indicator, theta)),
                                dependence = TRUE)
  #print(nw$terms[[termnumber]])
}

######################################################### 
InitErgmTerm.ttriple_ijk <- function(nw, arglist, ...) # Michael 
{
  n <- nw$gal$n # Number of nodes
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=FALSE,
    varnames = c("number", "indicator", "theta"),
    vartypes = c("numeric", "numeric", "numeric"),
    defaultvalues = list(NULL, NULL, NULL),
    required = c(FALSE, FALSE, FALSE))
  termnumber <- 1 + length(nw$terms)
  #print("InitErgmTerm.ttriple_ijk")
  indicator <- vector(mode = "numeric", length = n) # Category indicators  
  for (i in 1:length(indicator)) indicator[i] <- 1
  #print(indicator)
  number <- nw$max_number # (Maximum) number of categories
  #print(number)
  theta <- vector(mode = "numeric", length = number + 1) # Within- and between-category parameters
  for (i in 1:length(theta)) theta[i] <- 0 
  #print(theta)
  list(name = "ttriple_ijk", coef.names="ttriple_ijk", 
                                inputs = c(0, 1, 1+length(indicator)+length(theta), c(number, indicator, theta)),
                                dependence = TRUE)
  #print(nw$terms[[termnumber]])
}

######################################################### 
InitErgmTerm.ctriple_ijk <- function(nw, arglist, ...) # Michael 
{
  n <- nw$gal$n # Number of nodes
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=FALSE,
    varnames = c("number", "indicator", "theta"),
    vartypes = c("numeric", "numeric", "numeric"),
    defaultvalues = list(NULL, NULL, NULL),
    required = c(FALSE, FALSE, FALSE))
  termnumber <- 1 + length(nw$terms)
  #print("InitErgmTerm.ctriple_ijk")
  indicator <- vector(mode = "numeric", length = n) # Category indicators  
  for (i in 1:length(indicator)) indicator[i] <- 1
  #print(indicator)
  number <- nw$max_number # (Maximum) number of categories
  #print(number)
  theta <- vector(mode = "numeric", length = number + 1) # Within- and between-category parameters
  for (i in 1:length(theta)) theta[i] <- 0 
  #print(theta)
  list(name = "ctriple_ijk", coef.names="ctriple_ijk", 
                                inputs = c(0, 1, 1+length(indicator)+length(theta), c(number, indicator, theta)),
                                dependence = TRUE)
  #print(nw$terms[[termnumber]])
}

######################################################### 
InitErgmTerm.twostar_ijk <- function(nw, arglist, ...) # Michael 
{
  n <- nw$gal$n # Number of nodes
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=FALSE,
    varnames = c("number", "indicator", "theta"),
    vartypes = c("numeric", "numeric", "numeric"),
    defaultvalues = list(NULL, NULL, NULL),
    required = c(FALSE, FALSE, FALSE))
  termnumber <- 1 + length(nw$terms)
  #print("InitErgmTerm.twostar_ijk")
  indicator <- vector(mode = "numeric", length = n) # Category indicators  
  for (i in 1:length(indicator)) indicator[i] <- 1
  #print(indicator)
  number <- nw$max_number # (Maximum) number of categories
  #print(number)
  theta <- vector(mode = "numeric", length = number + 1) # Within- and between-category parameters
  for (i in 1:length(theta)) theta[i] <- 0 
  #print(theta)
  list(name = "twostar_ijk", coef.names="twostar_ijk", 
                                inputs = c(0, 1, 1+length(indicator)+length(theta), c(number, indicator, theta)),
                                dependence = TRUE)
  #print(nw$terms[[termnumber]])
}

