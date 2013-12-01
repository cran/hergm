
\name{bali}

\alias{bali}

\title{Terrorist network behind Bali bombing in 2002}

\usage{

data(bali)

}

\description{

The network corresponds to the contacts between the 17 terrorists who carried out the bombing in Bali, Indonesia in 2002.
The network is taken from Koschade (2006).

}

\value{Undirected network.}

\seealso{network, hergm, ergm.terms, hergm.terms, hergm.postprocess}

\references{

Koschade, S. (2006). A social network analysis of Jemaah Islamiyah: The applications to counter-terrorism and intelligence.
Studies in Conflict and Terrorism, 29, 559--575.

}

\examples{\dontrun{# Load undirected network with 17 nodes
data(bali)

# Hierarchical exponential-family model with stick-breaking prior
hergm(bali ~ edges_ij(5) + triangle_ijk(5))

}
}

