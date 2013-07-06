
\name{hergm-dirichlet}

\alias{hergm.dirichlet}

\title{Hierarchical Exponential-Family Random Graph Models: Sample Truncated Dirichlet Process Prior}

\usage{

hergm.dirichlet(n, 
                number, 
                alpha, 
                eta_mean, 
                eta_sd)

}

\description{

Sample truncated Dirichlet process prior.

}

\arguments{

\item{n}{number of nodes.}

\item{number}{number of blocks.}

\item{alpha}{scaling parameter.}

\item{eta_mean}{mean of Gaussian base distribution of Dirichlet process prior.}

\item{eta_sd}{standard deviation of Gaussian base distribution of Dirichlet process prior.}

}

\value{Indicators of block memberships of nodes \code{indicator} and block parameters \code{eta}.}

\examples{\dontrun{# Load undirected network with 15 nodes
hergm.dirichlet(n=100, number=100, alpha=1, eta_mean=-1, eta_sd=1)
}
}

