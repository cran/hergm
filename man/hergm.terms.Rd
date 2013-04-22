
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

\alias{triangle_ijk}

\alias{ttriple_ijk}

\alias{ctriple_ijk}

\docType{methods}

\title{Hierarchical Exponential-Family Random Graph Models: Terms}

\description{Hierarchical Exponential-Family Random Graph Models can be specified by calling the function \code{hergm(formula)},
  where formula is a formula of the form \code{network ~ terms}.

  By using suitable terms, it is possible to specify 

- the p_1 model for directed networks of Holland and Leinhardt (1981) and its extension to undirected random graph models with Dirichlet process priors (see arcs_i, arcs_j, mutual_i, and edges_i).
While the p_1 model for undirected and directed networks with parametric priors contains O(n) parameters (n = number of nodes) and therefore is not parsimonious,
the non-parametric Dirichlet process prior encourages a small number of unique parameters and therefore represents an elegant alternative to parametric priors.

- the stochastic block model of Snijders and Nowicki (1997) and Nowicki and Snijders (2001) with natural parameterization (restricted between-block parameters) and Dirichlet process priors (see edges_ij).

- hierarchical exponential-family models with stick-breaking priors (see mutual_ij, twostar_ijk, triangle_ijk, ttriple_ijk, ctriple_ijk).
Hierarchical exponential-family models replace the strong dependence of simple exponential-family models by weak dependence with an eye to solving the near-degeneracy problem of simple exponential-family model.

\code{hergm.terms} can be found here.
Additional terms, e.g. covariate-dependent terms, can be found in \code{ergm.terms}.
}

\arguments{

\item{edges_i(k) (undirected network)}{adding the term \code{edges_i} to the model adds node-dependent edge terms to the model;
the optional argument \code{k} is the maximum number of blocks (default: the number of nodes).
}

\item{arcs_i(k) (directed network)}{adding the term \code{arcs_i} to the model adds node-dependent outdegree terms to the model;
the optional argument \code{k} is the maximum number of blocks (default: the number of nodes).
}

\item{arcs_j(k) (directed network)}{adding the term \code{arcs_j} to the model adds node-dependent indegree terms to the model;
the optional argument \code{k} is the maximum number of blocks (default: the number of nodes).
}

\item{edges_ij(k) (undirected, directed network)}{adding the term \code{edges_ij} to the model adds block-dependent edge terms to the model.
the optional argument \code{k} is the maximum number of blocks (default: the number of nodes).
}

\item{mutual_i(k) (directed network)}{adding the term \code{mutual_i} to the model adds additive, block-dependent mutual edge terms to the model.
the optional argument \code{k} is the maximum number of blocks (default: the number of nodes).
}

\item{mutual_ij(k) (directed network)}{adding the term \code{mutual_ij} to the model adds block-dependent mutual edge terms to the model.
the optional argument \code{k} is the maximum number of blocks (default: the number of nodes).
}

\item{twostar_ijk(k) (undirected network)}{adding the term \code{twostar_ijk} to the model adds block-dependent two-star terms to the model.
the optional argument \code{k} is the maximum number of blocks (default: the number of nodes).
}

\item{triangle_ijk(k) (undirected, directed network)}{adding the term \code{triangle_ijk} to the model adds block-dependent triangle terms to the model.
the optional argument \code{k} is the maximum number of blocks (default: the number of nodes).
}

\item{ttriple_ijk(k) (directed network)}{adding the term \code{ttriple_ijk} to the model adds block-dependent transitive triple terms to the model.
the optional argument \code{k} is the maximum number of blocks (default: the number of nodes).
}
 
\item{ctriple_ijk(k) (directed network)}{adding the term \code{ctriple_ijk} to the model adds block-dependent cyclic triple terms to the model.
the optional argument \code{k} is the maximum number of blocks (default: the number of nodes).
}

}

\references{

Holland, P. W. and S. Leinhardt (1981). An exponential family of probability distributions for directed graphs. Journal of the American Statistical
  Association 76 (373), 33--65.

Nowicki, K. and T. A. B. Snijders (2001). Estimation and prediction for
  stochastic blockstructures. Journal of the American Statistical Association 96 (455), 1077--1087.

Snijders, T. A. B. and K. Nowicki (1997). Estimation and prediction for
  stochastic blockmodels for graphs with latent block structure. Journal of
  Classification 14, 75--100.

}

\seealso{network, hergm, ergm.terms, hergm.postprocess}

\examples{\dontrun{# Load undirected network with 15 nodes (see ?example)
data(example)
# p_1 model for undirected network with Dirichlet process prior
hergm(d ~ edges_i)

# Load directed network with 18 nodes (see ?sampson)
data(sampson) 
# p_1 model for directed network with Dirichlet process prior
hergm(samplike ~ arcs_i + arcs_j + mutual)

# Load undirected network with 15 nodes (see ?example)
data(example)
# Stochastic block model for undirected network 
# with natural parameterization and Dirichlet process prior
hergm(d ~ edges_ij)

# Load directed network with 18 nodes (see ?sampson)
data(sampson) 
# Stochastic block model for directed network 
# with natural parameterization and Dirichlet process prior
hergm(samplike ~ edges_ij + mutual)

# Load undirected network with 15 nodes (see ?example)
data(example)
# Hierarchical exponential-family model with stick-breaking prior
hergm(d ~ edges + mutual + ttriple_ijk)

# Load directed network with 18 nodes (see ?sampson)
data(sampson) 
# Hierarchical exponential-family model with stick-breaking prior
hergm(samplike ~ edges + mutual + ttriple_ijk)
}
}

