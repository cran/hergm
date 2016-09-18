\name{hergm-terms}

\alias{hergm.terms}

\alias{terms-hergm}

\alias{terms.hergm}

\alias{edges_i}

\alias{arcs_i}

\alias{arcs_j}

\alias{edges_ij}

\alias{mutual_i}

\alias{mutual_ij}

\alias{twostar_ijk}

\alias{transitiveties_ijk}

\alias{triangle_ijk}

\alias{ttriple_ijk}

\alias{ctriple_ijk}

\docType{methods}

\title{Model terms}

\description{

Hierarchical exponential-family random graph models with local dependence can be specified by calling the function \code{hergm(formula)},
where formula is a formula of the form \code{network ~ terms}.
By specifying suitable terms, 
it is possible to specify a wide range of models: see \code{hergm}.
\code{hergm.terms} can be found here.
In addition, \code{ergm.terms} can be used to include covariates.

}

\arguments{

\item{edges_i (undirected network)}{adding the term \code{edges_i} to the model adds node-dependent edge terms to the model.}

\item{arcs_i (directed network)}{adding the term \code{arcs_i} to the model adds node-dependent outdegree terms to the model.}

\item{arcs_j (directed network)}{adding the term \code{arcs_j} to the model adds node-dependent indegree terms to the model.}

\item{edges_ij (undirected, directed network)}{adding the term \code{edges_ij} to the model adds block-dependent edge terms to the model.}

\item{mutual_i (directed network)}{adding the term \code{mutual_i} to the model adds additive, block-dependent mutual edge terms to the model.}

\item{mutual_ij (directed network)}{adding the term \code{mutual_ij} to the model adds block-dependent mutual edge terms to the model.}

\item{twostar_ijk (undirected network)}{adding the term \code{twostar_ijk} to the model adds block-dependent two-star terms to the model;
}

\item{transitiveties_ijk (directed network)}{adding the term \code{transitiveties_ijk} to the model adds block-dependent transitive ties terms to the model.}

\item{triangle_ijk (undirected, directed network)}{adding the term \code{triangle_ijk} to the model adds block-dependent triangle terms to the model.}

\item{ttriple_ijk (directed network)}{adding the term \code{ttriple_ijk} to the model adds block-dependent transitive triple terms to the model.}
 
\item{ctriple_ijk (directed network)}{adding the term \code{ctriple_ijk} to the model adds block-dependent cyclic triple terms to the model.}

}

\references{

Handcock, M. S. (2003). Assessing degeneracy in statistical models of social networks. Technical report, Center for Statistics and the Social Sciences, University of Washington, Seattle, http://www.csss.washington.edu/Papers.

Holland, P. W. and S. Leinhardt (1981). An exponential family of probability distributions for directed graphs. Journal of the American Statistical Association, Theory \& Methods, 76, 33--65.

Nowicki, K. and T. A. B. Snijders (2001). Estimation and prediction for stochastic blockstructures. Journal of the American Statistical Association, Theory \& Methods, 96, 1077--1087.

Snijders, T. A. B. and K. Nowicki (1997). Estimation and prediction for stochastic blockmodels for graphs with latent block structure. Journal of Classification 14, 75--100.

Schweinberger, M. (2011). Instability, sensitivity, and degeneracy of discrete exponential families. Journal of the American Statistical Association, Theory & Methods, 106, 1361--1370.

Schweinberger, M. and M. S. Handcock (2015). Local dependence in random graph models: characterization, properties, and statistical Inference. Journal of the Royal Statistical Society, Series B (Statistical Methodology), 7, 647-676.

Schweinberger, M., Petrescu-Prahova, M. and D. Q. Vu (2014). Disaster response on September 11, 2001 through the lens of statistical network analysis. Social Networks, 37, 42--55.

Vu, D. Q., Hunter, D. R. and M. Schweinberger (2013). Model-based clustering of large networks. Annals of Applied Statistics, 7, 1010--1039.

}

\seealso{hergm, ergm.terms}

\examples{
\dontrun{
data(example)
# p_1 model: undirected network
hergm(d ~ edges_i)

data(sampson) 
# p_1 model: directed network
hergm(samplike ~ arcs_i + arcs_j + mutual)

data(example)
# Stochastic block model: undirected network 
hergm(d ~ edges_ij)

data(sampson) 
# Stochastic block model: directed network 
hergm(samplike ~ edges_ij + mutual)

data(example)
# Exponential-family random graph model with local dependence: undirected network
hergm(d ~ edges_ij + triangle_ijk)

data(sampson) 
# Exponential-family random graph model with local dependence: directed network
hergm(samplike ~ edges + mutual + ttriple_ijk)
}
}

