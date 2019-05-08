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

\item{edges_i (undirected network)}{adding the term \code{edges_i} to the model adds node-dependent edge terms to the model;
please note: the term \code{edges_i} can be used with \code{method = "bayes"} but cannot be used with the default \code{method = "ml"}.}

\item{arcs_i (directed network)}{adding the term \code{arcs_i} to the model adds node-dependent outdegree terms to the model;
please note: the term \code{arcs_i} can be used with \code{method = "bayes"} but cannot be used with the default \code{method = "ml"}.}

\item{arcs_j (directed network)}{adding the term \code{arcs_j} to the model adds node-dependent indegree terms to the model;
please note: the term \code{arcs_j} can be used with \code{method = "bayes"} but cannot be used with the default \code{method = "ml"}.}

\item{edges_ij (undirected, directed network)}{adding the term \code{edges_ij} to the model adds block-dependent edge terms to the model.}

\item{mutual_i (directed network)}{adding the term \code{mutual_i} to the model adds additive, block-dependent mutual edge terms to the model.}

\item{mutual_ij (directed network)}{adding the term \code{mutual_ij} to the model adds block-dependent mutual edge terms to the model.}

\item{twostar_ijk (undirected network)}{adding the term \code{twostar_ijk} to the model adds block-dependent two-star terms to the model;
}

\item{transitiveties_ijk (directed network)}{adding the term \code{transitiveties_ijk} to the model adds block-dependent transitive ties terms to the model.}

\item{triangle_ijk (undirected, directed network)}{adding the term \code{triangle_ijk} to the model adds block-dependent triangle terms to the model.}

\item{ttriple_ijk (directed network)}{adding the term \code{ttriple_ijk} to the model adds block-dependent transitive triple terms to the model;
please note: the term \code{ttriple_ijk} can be used with \code{method = "bayes"} but cannot be used with the default \code{method = "ml"}.}
 
\item{ctriple_ijk (directed network)}{adding the term \code{ctriple_ijk} to the model adds block-dependent cyclic triple terms to the model;
please note: the term \code{ctriple_ijk} can be used with \code{method = "bayes"} but cannot be used with the default \code{method = "ml"}.}

}

\references{

Schweinberger, M. and P. Luna (2018). HERGM: Hierarchical exponential-family random graph models. Journal of Statistical Software, 85, 1--39.

}

\seealso{hergm, ergm.terms}

\examples{
\dontrun{data(sampson) 
hergm(samplike ~ edges + transitiveties_ijk)
}
}
