\name{bali}

\alias{bali}

\title{Bali terrorist network}

\usage{

data(bali)

}

\description{

The network corresponds to the contacts between the 17 terrorists who carried out the bombing in Bali, Indonesia in 2002.
The network is taken from Koschade (2006).

}

\value{Undirected network.}

\seealso{network, hergm, ergm.terms, hergm.terms}

\references{

Koschade, S. (2006). A social network analysis of Jemaah Islamiyah: The applications to counter-terrorism and intelligence.
Studies in Conflict and Terrorism, 29, 559--575.

}

\examples{
\dontrun{data(bali)
object <- hergm(bali ~ edges_ij + triangle_ijk, max_number=3, method="bayes")
summary(object)
}
}

