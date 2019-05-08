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
\dontrun{data(example)
object <- hergm(d ~ edges_i, max_number=3)
summary(object)
}
}

