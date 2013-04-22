
\name{hergm}

\alias{hergm}

\title{Hierarchical Exponential-Family Random Graph Models: Simulation and Bayesian inference}

\description{

The package \code{hergm} implements exponential-family random graph models with Dirichlet process / stick-breaking priors,
including

- the p_1 model for directed networks of Holland and Leinhardt (1981) and its extension to undirected random graph models with Dirichlet process priors (see arcs_i, arcs_j, mutual_i, and edges_i).
While the p_1 model for undirected and directed networks with parametric priors contains O(n) parameters (n = number of nodes) and therefore is not parsimonious,
the non-parametric Dirichlet process prior encourages a small number of unique parameters and therefore represents an elegant alternative to parametric priors.

- the stochastic block model of Snijders and Nowicki (1997) and Nowicki and Snijders (2001) with natural parameterization (restricted between-block parameters) and Dirichlet process priors (see edges_ij).

- hierarchical exponential-family models with stick-breaking priors (see mutual_ij, twostar_ijk, triangle_ijk, ttriple_ijk, ctriple_ijk).
Hierarchical exponential-family models replace the strong dependence of simple exponential-family models by weak dependence with an eye to solving the near-degeneracy problem of simple exponential-family model.

The package \code{hergm} implements simulation and Bayesian inference for the mentioned models.

}

\usage{

hergm(formula,
      alpha = NULL,
      alpha_shape = NULL,
      alpha_rate = NULL,
      eta = NULL,
      eta_mean = NULL,
      eta_sd = NULL,
      eta_mean_mean = NULL,
      eta_mean_sd = NULL,
      eta_precision_shape = NULL,
      eta_precision_rate = NULL,
      indicator = NULL,
      parallel = 1,
      simulate = FALSE,
      seeds = NULL,
      samplesize = 1e+5,
      burnin = 1e+4,
      interval = 1e+2,
      mh_scale = NULL,
      output = TRUE,
      verbose = -1,
      name = NULL, 
      \dots)

}

\arguments{

\item{formula}{formula of the form \code{network ~ terms}.
          Networks can be created by calling the function \code{network}.
          Possible terms can be found in \code{ergm.terms} and \code{hergm.terms}.
}

\item{alpha}{scaling parameter of truncated Dirichlet process / stick-breaking prior of natural parameters of exponential-family model.}

\item{alpha_shape, alpha_rate}{shape and rate (inverse scale) parameter of Gamma prior of scaling parameter.}

\item{eta}{natural parameters of exponential-family model.}

\item{eta_mean, eta_sd}{means and standard deviations of Gaussian baseline distribution of Dirichlet process / stick-breaking prior of natural parameters.}

\item{eta_mean_mean, eta_mean_sd}{means and standard deviations of Gaussian prior of mean of Gaussian baseline distribution of Dirichlet process / stick-breaking prior.}

\item{eta_precision_shape, eta_precision_rate}{shape and rate (inverse scale) parameter of Gamma prior of precision parameter of Gaussian baseline distribution of Dirichlet process / stick-breaking prior.}

\item{indicator}{if \code{simulate == TRUE} and indicator of block membership is specified, condition on specified indicator (of block membership).}

\item{parallel}{number of processors; if more than one, computing is parallel.}

\item{simulate}{if \code{TRUE}, simulation of networks, otherwise Bayesian inference.}

\item{seeds}{seed of pseudo-random number generator; if computing is parallel, number of seeds must equal number of processors.}

\item{samplesize}{if \code{simulate == TRUE}, number of networks to be sampled, otherwise number of draws from posterior; if computing is parallel, number of draws per processors.}

\item{burnin}{if \code{simulate == TRUE}, number of burn-in iterations.}

\item{interval}{if \code{simulate == TRUE}, number of proposals between sampled networks.}

\item{mh_scale}{if \code{simulate == FALSE}, scale factor of candicate-generating distribution of Metropolis-Hastings algorithm.}

\item{output}{if \code{TRUE}, full output, otherwise limited output.}

\item{name}{name of project; if \code{output == TRUE}, name of project is used to name output files.}

\item{verbose}{console output: -1: no output; 0: short output; +1: long output.}

\item{\dots}{additional arguments, to be passed to lower-level functions in the future.}

}

\value{If called with the option \code{simulate = TRUE}, 
the function \code{hergm} returns a sample of networks, 
otherwise a raw MCMC sample from the posterior.
To postprocess the sample,
call the function \code{hergm.postprocess}.}

\seealso{network, ergm.terms, hergm.terms, hergm.postprocess}

\references{

Holland, P. W. and S. Leinhardt (1981). An exponential family of probability distributions for directed graphs. Journal of the American Statistical
  Association 76 (373), 33--65.

Nowicki, K. and T. A. B. Snijders (2001). Estimation and prediction for
  stochastic blockstructures. Journal of the American Statistical Association 96 (455), 1077--1087.

Snijders, T. A. B. and K. Nowicki (1997). Estimation and prediction for
  stochastic blockmodels for graphs with latent block structure. Journal of
  Classification 14, 75--100.

}

\examples{\dontrun{# Load undirected network with 15 nodes (see ?example)
data(example)

# p_1 model for undirected network with Dirichlet process prior
hergm(d ~ edges_i)

# Stochastic block model for undirected network 
# with natural parameterization and Dirichlet process prior
hergm(d ~ edges_ij)

# Hierarchical exponential-family model with stick-breaking prior
hergm(d ~ edges + triangle_ijk)
}
}

