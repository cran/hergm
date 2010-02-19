
\name{hergm-postprocess}

\alias{hergm.postprocess}

\alias{postprocess-hergm}

\alias{postprocess.hergm}

\title{Hierarchical Exponential-Family Random Graph Models: Postprocessing of Bayesian MCMC samples}

\usage{

hergm.postprocess(n = NULL,
                  k = NULL, 
                  d1 = 0, 
                  d2 = 0, 
                  burnin = NULL, 
                  samplesize = NULL, 
                  mcmc = NULL, 
                  output = TRUE,
                  name = NULL, 
                  \dots)

}

\description{

If called with the option \code{simulate = FALSE},
the function \code{hergm} returns a raw MCMC sample of parameters from the posterior in the form of a vector.
The function \code{hergm.process} postprocesses the raw MCMC sample:
if called with the \code{output = FALSE},
\code{hergm.postprocess} extracts information of interest,
otherwise it solves, in addition, the so-called label-switching problem.
The label-switching problem is rooted in the invariance of the likelihood function to permutations of the labels of blocks, and implies that the raw MCMC sample cannot be used to infer to block-dependent entities.
The label-switching problem can be solved in a Bayesian decision-theoretic framework: by defining a loss function and minimizing the posterior expected loss.
Calling \code{hergm.process} minimizes the posterior expected loss using a simple and convenient loss function.
The required computations can be time-consuming when the number of blocks \code{k} is large.

}

\arguments{

\item{n}{number of nodes.}

\item{k}{number of blocks.}

\item{d1}{number of ergm terms.}

\item{d2}{number of hergm terms.}

\item{burnin}{number of burn-in iterations.}

\item{samplesize}{total sample size, including number of burn-in iterations and, if parallel > 1, summed across computing nodes.}

\item{mcmc}{MCMC sample in the form of vector.}

\item{output}{if TRUE, full output, including relabeled MCMC sample; otherwise limited output.}

\item{name}{name of project; if output = TRUE, name of project is used to name output files.}

\item{\dots}{additional arguments, to be passed to lower-level functions in the future.}

}

\value{Postprocessed MCMC sample.}

\seealso{network, hergm.example, hergm, ergm.terms, hergm.terms}

\examples{\dontrun{# Load undirected network with 15 nodes (see ?example)
data(example)

# Number of blocks: truncation value of Dirichlet process prior
k <- 10

# Generate MCMC sample: total samplesize is 1 * 1,000 = 1,000
mcmc <- hergm(d ~ edges_i, alpha_shape = 1, alpha_rate = 1, 
eta_mean = -1, eta_sd = 2, seed = 2010, parallel = 1, samplesize = 1000)
# Postprocess MCMC sample of size 1,000
processed_mcmc <- hergm.postprocess(n = 15, k = k, d1 = 0, d2 = 1, 
burnin = 0, samplesize = 1000, mcmc = mcmc)

# Generate MCMC sample: total samplesize is 20 * 100 = 2,000
mcmc <- hergm(d ~ edges + triangle_ijk, alpha_shape = 1, alpha_rate = 1, 
eta_mean = c(-1, 0), eta_sd = c(2, 1), seed = 2010, parallel = 20, samplesize = 100)
# Postprocess MCMC sample of size 2,000
processed_mcmc <- hergm.postprocess(n = 15, k = k, d1 = 0, d2 = 1, 
burnin = 0, samplesize = 2000, mcmc = mcmc)
}
}
