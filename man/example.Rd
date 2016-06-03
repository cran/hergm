
\name{example}

\alias{example}

\alias{d}

\title{Example network}

\usage{

data(example)

}

\description{

Example data set: synthetic, undirected network with 15 nodes.

}

\value{Undirected network.}

\seealso{network, hergm, ergm.terms, hergm.terms}

\examples{
\dontrun{
data(example)

hergm(d ~ edges_i)

hergm(d ~ edges_ij + triangle_ijk)
}
}

