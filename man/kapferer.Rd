\name{kapferer}

\alias{kapferer}

\title{Kapferer collaboration network}

\usage{

data(kapferer)

}

\description{

The network corresponds to collaborations between 39 workers in a tailor shop in Africa:
an undirected edge between workers \code{i} and \code{j} indicates that the workers collaborated.
The network is taken from Kapferer (1972).

}

\value{Undirected network.}

\seealso{network, hergm, ergm.terms, hergm.terms}

\references{

Kapferer, B. (1972). Strategy and Transaction in an African Factory. Manchester University
Press, Manchester, U.K.

}

\examples{
\dontrun{data(kapferer)
object <- hergm(kapferer ~ edges_ij + transitiveties_ijk, max_number=3)
summary(object)
}
}

