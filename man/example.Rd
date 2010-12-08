
\name{example}

\alias{example}

\alias{d}

\title{Example data set}

\usage{

data(example)

}

\description{

Example data set: synthetic, undirected network with 15 nodes.

}

\value{Undirected network.}

\seealso{network, hergm, ergm.terms, hergm.terms, hergm.postprocess}

\examples{\dontrun{# Load undirected network with 15 nodes
data(example)

# p_1 model for undirected network with Dirichlet process prior
hergm(d ~ edges_i)

# Hierarchical exponential-family model with stick-breaking prior
hergm(d ~ edges + triangle_ijk)
}
}

