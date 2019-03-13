\name{bunt}

\alias{bunt}

\title{Van de Bunt friendship network}

\usage{

data(bunt)

}

\description{

Van de Bunt (1999) and Van de Bunt et al. (1999) collected data on friendships between 32 freshmen at a European university at 7 time points.
Here,
the last time point is used.
A directed edge from student \code{i} to \code{j} indicates that student \code{i} considers student \code{j} to be a ``friend" or ``best friend".

}

\value{Directed network.}

\seealso{network, hergm, ergm.terms, hergm.terms}

\references{

Van de Bunt, G. G. (1999). Friends by choice. An Actor-Oriented Statistical Network Model for Friendship Networks through Time. Thesis Publishers,
Amsterdam.

Van de Bunt, G. G., Van Duijn, M. A. J., and T. A. B. Snijders (1999). Friendship Networks Through Time: An Actor-Oriented Statistical Network Model. Computational and Mathematical Organization Theory, 5, 167--192.

}

\examples{
\dontrun{data(bunt)
hergm(bunt ~ edges_ij + ttriple_ijk)
}
}

