
\name{hergm}

\alias{hergm}

\title{Estimate and simulate hierarchical exponential-family random graph models with local dependence}

\description{

The function \code{hergm} estimates and simulates three classes of hierarchical exponential-family random graph models:

1. The p_1 model of Holland and Leinhardt (1981) in exponential-family form and extensions by Vu, Hunter, and Schweinberger (2013) and Schweinberger, Petrescu-Prahova, and Vu (2014) to both directed and undirected random graphs with additional model terms, with and without covariates, and with parametric and nonparametric priors (see \code{arcs_i}, \code{arcs_j}, \code{edges_i}, \code{edges_ij}, \code{mutual_i}, \code{mutual_ij}).

2. The stochastic block model of Snijders and Nowicki (1997) and Nowicki and Snijders (2001) in exponential-family form and extensions by Vu, Hunter, and Schweinberger (2013) and Schweinberger, Petrescu-Prahova, and Vu (2014) with additional model terms, with and without covariates, and with parametric and nonparametric priors (see \code{arcs_i}, \code{arcs_j}, \code{edges_i}, \code{edges_ij}, \code{mutual_i}, \code{mutual_ij}).

3. The exponential-family random graph models with local dependence of Schweinberger and Handcock (2015), with and without covariates, and with parametric and nonparametric priors (see \code{arcs_i}, \code{arcs_j}, \code{edges_i}, \code{edges_ij}, \code{mutual_i}, \code{mutual_ij}, \code{twostar_ijk}, \code{triangle_ijk}, \code{ttriple_ijk}, \code{ctriple_ijk}).
The exponential-family random graph models with local dependence replace the long-range dependence of conventional exponential-family random graph models by short-range dependence. 
Therefore, the exponential-family random graph models with local dependence replace the strong dependence of conventional exponential-family random graph models by weak dependence, 
reducing the problem of model degeneracy (Handcock, 2003; Schweinberger, 2011) and improving goodness-of-fit (Schweinberger and Handcock, 2015).

The function \code{hergm} postprocesses the output returned by the function \code{hergm}.
If \code{hergm} is called with \code{relabel > 0},
it solves the so-called label-switching problem.
The label-switching problem is rooted in the invariance of the likelihood function to permutations of the labels of blocks, and implies that raw MCMC samples from the posterior cannot be used to infer to block-dependent entities.
The label-switching problem can be solved in a Bayesian decision-theoretic framework: by choosing a loss function and minimizing the posterior expected loss.
Two loss functions are implemented in \code{hergm}, the loss function of Schweinberger and Handcock (2015) (\code{relabel == 1}) and the loss function of Peng and Carvalho (2015) (\code{relabel == 2}).
The first loss function seems to be superior in terms of the reported clustering probabilities, but is more expensive in terms of computing time.
A rule of thumb is to use the first loss function when \code{max_number < 15} and use the second loss function otherwise.

}

\usage{

hergm(formula,
                  max_number = NULL,
                  hierarchical = TRUE,
                  parametric = FALSE,
                  initialize = FALSE,
                  perturb = FALSE,
                  scaling = NULL,
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
                  mean_between = NULL,
                  all_indicators_fixed = FALSE,
                  indicators_fixed = FALSE,
                  indicator = NULL,
                  parallel = 1,
                  simulate = FALSE,
                  seeds = NULL,
                  samplesize = 1e+5,
                  interval = 1024,
                  burnin = 16*interval,
                  mh.scale = 0.25,
                  variational = FALSE,
                  temperature = c(1,100),
                  predictions = FALSE,
                  posterior.burnin = 0,
                  posterior.thinning = 1,
                  relabel = 0,
                  number.runs = 1,
                  verbose = 1,
                  \dots)

}

\arguments{

\item{formula}{formula of the form \code{network ~ terms}.
          Networks can be created by calling the function \code{network}.
          Possible terms can be found in \code{ergm.terms} and \code{hergm.terms}.
}

\item{max_number}{maximum number of blocks.}

\item{hierarchical}{hierarchical prior; if \code{hierarchical == TRUE}, prior is hierarchical (i.e., the means and variances of block parameters are governed by a hyper-prior), otherwise non-hierarchical (i.e., the means and variances of block parameters are fixed).}

\item{parametric}{parametric prior; if \code{parametric == FALSE}, prior is truncated Dirichlet process prior, otherwise parametric Dirichlet prior.}

\item{initialize}{if \code{initialize == TRUE}, initialize block memberships by spectral clustering.}

\item{perturb}{if \code{initialize == TRUE} and \code{perturb == TRUE}, initialize block memberships by spectral clustering and perturb.}

\item{scaling}{if \code{scaling == TRUE}, use size-dependent parameterizations which ensure that the scaling of between- and within-neighborhood terms is consistent with sparse edge terms.}

\item{alpha}{concentration parameter of truncated Dirichlet process prior of natural parameters of exponential-family model.}

\item{alpha_shape, alpha_rate}{shape and rate parameter of Gamma prior of scaling parameter.}

\item{eta}{natural parameters of exponential-family random graph model.}

\item{eta_mean, eta_sd}{means and standard deviations of Gaussian baseline distribution of Dirichlet process prior of natural parameters.}

\item{eta_mean_mean, eta_mean_sd}{means and standard deviations of Gaussian prior of mean of Gaussian baseline distribution of Dirichlet process prior.}

\item{eta_precision_shape, eta_precision_rate}{shape and rate (inverse scale) parameter of Gamma prior of precision parameter of Gaussian baseline distribution of Dirichlet process prior.}

\item{mean_between}{if \code{simulate == TRUE} and \code{eta == NULL}, then \code{mean_between} specifies the mean-value parameter of edges between blocks.}

\item{all_indicators_fixed}{indicates whether all indicators of the block memberships are fixed at the specified indicators; if some block memberships are unsspecified, spectral clustering is used to initialize all block memberships.}

\item{indicators_fixed}{indicates whether the indicators of the block memberships are fixed at the specified indicators.}

\item{indicator}{if the indicators of block memberships are specified as numbers between 1 and \code{max_number}, the specified indicators are either used as starting values (\code{indicator_fixed == FALSE}) or the indicators are fixed at the specified indicators (\code{indicator_fixed == TRUE}), which is useful when indicators of block memberhips are observed (which is the case in multilevel networks).}

\item{parallel}{number of computing nodes; if \code{parallel > 1}, \code{hergm} is run on \code{parallel} computing nodes.}

\item{simulate}{if \code{simulate == TRUE}, simulation of networks, otherwise Bayesian inference.}

\item{seeds}{seed of pseudo-random number generator; if \code{parallel > 1}, number of seeds must equal number of computing nodes.}

\item{samplesize}{if \code{simulate == TRUE}, number of network draws, otherwise number of posterior draws; if \code{parallel > 1}, number of draws on each computing node.}

\item{interval}{if \code{simulate == TRUE}, number of proposals between sampled networks.}

\item{burnin}{if \code{simulate == TRUE}, number of burn-in iterations.}

\item{mh.scale}{if \code{simulate == FALSE}, scale factor of candicate-generating distribution of Metropolis-Hastings algorithm.}

\item{variational}{if \code{simulate == FALSE} and \code{variational == TRUE}, variational methods are used to construct the proposal distributions of block memberships; limited to selected models.}

\item{temperature}{if \code{simulate == FALSE} and \code{variational == TRUE}, minimum and maximum temperature; the temperature is used to melt down the proposal distributions of indicators, which are based on the full conditional distributions of indicators but can have low entropy, resulting in slow mixing of the Markov chain; the temperature is a function of the entropy of the full conditional distributions and is designed to increase the entropy of the proposal distributions, and the minimum and maximum temperature are user-defined lower and upper bounds on the temperature.}

\item{predictions}{if \code{predictions == TRUE} and \code{simulate == FALSE}, returns posterior predictions of statistics in the model.}

\item{posterior.burnin}{number of burn-in iterations; if computing is parallel, number of burn-in iterations per processor.}

\item{posterior.thinning}{if \code{thinning > 1}, every \code{thinning}-th sample point is used while all others discarded.}

\item{relabel}{if \code{relabel > 0}, relabel MCMC sample by minimizing the posterior expected loss of Schweinberger and Handcock (2015) (\code{relabel == 1}) or Peng and Carvalho (2015) (\code{relabel == 2}).}

\item{number.runs}{if \code{relabel == 1}, number of runs of relabeling algorithm.}

\item{verbose}{if \code{verbose == -1}, no console output; if \code{verbose == 0}, short console output; if \code{verbose == +1}, long console output.}

\item{\dots}{additional arguments, to be passed to lower-level functions in the future.}

}

\value{

If called with the option \code{simulate == TRUE}, 
the function \code{hergm} returns a sample of networks, 
otherwise a MCMC sample from the posterior.

\item{ergm_theta}{parameters of \code{ergm-terms}.}

\item{alpha}{concentration parameter of truncated Dirichlet process prior of parameters of \code{hergm-terms}.}

\item{eta_mean}{mean parameters of Gaussian base distribution of parameters of \code{hergm-terms}.}

\item{eta_precision}{precision parameters of Gaussian base distribution of parameters of \code{hergm-terms}.}

\item{hergm_theta}{parameters of \code{hergm-terms}.}

\item{loss}{if \code{relabel == TRUE}, local minimum of loss function.}

\item{p_k}{probabilities of membership to blocks.}

\item{indicator}{indicators of memberships of nodes.}

\item{p_i_k}{probabilities of membership of nodes to blocks.}

\item{prediction}{posterior predictions of statistics.}

}

\seealso{network, ergm.terms, hergm.terms, hergm.gof, hergm.plot, summary}

\references{

Handcock, M. S. (2003). Assessing degeneracy in statistical models of social networks. Technical report, Center for Statistics and the Social Sciences, University of Washington, Seattle, http://www.csss.washington.edu/Papers.

Holland, P. W. and S. Leinhardt (1981). An exponential family of probability distributions for directed graphs. Journal of the American Statistical Association, Theory \& Methods, 76, 33--65.

Nowicki, K. and T. A. B. Snijders (2001). Estimation and prediction for stochastic blockstructures. Journal of the American Statistical Association, Theory \& Methods, 96, 1077--1087.

Peng, L. and L. Carvalho (2015). Bayesian degree-corrected stochastic block models for community detection. Technical report, Boston University, arXiv:1309.4796v1.

Snijders, T. A. B. and K. Nowicki (1997). Estimation and prediction for stochastic blockmodels for graphs with latent block structure. Journal of Classification 14, 75--100.

Schweinberger, M. (2011). Instability, sensitivity, and degeneracy of discrete exponential families. Journal of the American Statistical Association, Theory & Methods, 106, 1361--1370.

Schweinberger, M. and M. S. Handcock (2015). Local dependence in random graph models: characterization, properties, and statistical Inference. Journal of the Royal Statistical Society, Series B (Statistical Methodology), 7, 1-30, in press.

Schweinberger, M., Petrescu-Prahova, M. and D. Q. Vu (2014). Disaster response on September 11, 2001 through the lens of statistical network analysis. Social Networks, 37, 42--55.

Vu, D. Q., Hunter, D. R. and M. Schweinberger (2013). Model-based clustering of large networks. Annals of Applied Statistics, 7, 1010--1039.

}

\examples{
\dontrun{
data(example)

hergm(d ~ edges_i)

hergm(d ~ edges_ij)

hergm(d ~ edges_ij + triangle_ijk)

data(sampson)

hergm(samplike ~ arcs_i + arcs_j)

hergm(samplike ~ edges_ij + mutual_ij)

hergm(samplike ~ edges_ij + mutual_ij + ttriple_ijk)
}
}

