\name{hergm}

\alias{hergm}

\title{Hierarchical exponential-family random graph models with local dependence}

\description{

The function \code{hergm} estimates and simulates three classes of hierarchical exponential-family random graph models:

1. The p_1 model of Holland and Leinhardt (1981) in exponential-family form and extensions by Vu, Hunter, and Schweinberger (2013) and Schweinberger, Petrescu-Prahova, and Vu (2014) to both directed and undirected random graphs with additional model terms, with and without covariates, and with parametric and nonparametric priors (see \code{arcs_i}, \code{arcs_j}, \code{edges_i}, \code{edges_ij}, \code{mutual_i}, \code{mutual_ij}).

2. The stochastic block model of Snijders and Nowicki (1997) and Nowicki and Snijders (2001) in exponential-family form and extensions by Vu, Hunter, and Schweinberger (2013) and Schweinberger, Petrescu-Prahova, and Vu (2014) with additional model terms, with and without covariates, and with parametric and nonparametric priors (see \code{arcs_i}, \code{arcs_j}, \code{edges_i}, \code{edges_ij}, \code{mutual_i}, \code{mutual_ij}).

3. The exponential-family random graph models with local dependence of Schweinberger and Handcock (2015), with and without covariates, and with parametric and nonparametric priors (see \code{arcs_i}, \code{arcs_j}, \code{edges_i}, \code{edges_ij}, \code{mutual_i}, \code{mutual_ij}, \code{twostar_ijk}, \code{triangle_ijk}, \code{ttriple_ijk}, \code{ctriple_ijk}).
The exponential-family random graph models with local dependence replace the long-range dependence of conventional exponential-family random graph models by short-range dependence. 
Therefore, the exponential-family random graph models with local dependence replace the strong dependence of conventional exponential-family random graph models by weak dependence, 
reducing the problem of model degeneracy (Handcock, 2003; Schweinberger, 2011) and improving goodness-of-fit (Schweinberger and Handcock, 2015).

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
      indicator = NULL,
      parallel = 1,
      simulate = FALSE,
      seeds = NULL,
      sample_size = 1e+5,
      interval = 1024,
      burnin = 16*interval,
      mh.scale = 0.25,
      variational = FALSE,
      temperature = c(1,100),
      predictions = FALSE,
      posterior.burnin = 2000,
      posterior.thinning = 1,
      relabel = 1,
      number_runs = 1,
      verbose = 0,
      \dots)

}

\arguments{

\item{formula}{formula of the form \code{network ~ terms}.
          \code{network} is an object of class \code{network} and can be created by calling the function \code{network}.
          Possible terms can be found in \code{ergm.terms} and \code{hergm.terms}.}

\item{max_number}{maximum number of blocks.}

\item{hierarchical}{hierarchical prior; if \code{hierarchical == TRUE}, prior is hierarchical (i.e., the means and variances of block parameters are governed by a hyper-prior), otherwise non-hierarchical (i.e., the means and variances of block parameters are fixed).}

\item{parametric}{parametric prior; if \code{parametric == FALSE}, prior is truncated Dirichlet process prior, otherwise parametric Dirichlet prior.}

\item{initialize}{if \code{initialize == TRUE}, initialize block memberships of nodes by spectral clustering.}

\item{perturb}{if \code{initialize == TRUE} and \code{perturb == TRUE}, initialize block memberships of nodes by spectral clustering and perturb.}

\item{scaling}{if \code{scaling == TRUE}, use size-dependent parameterizations which ensure that the scaling of between- and within-neighborhood terms is consistent with sparse edge terms.}

\item{alpha}{concentration parameter of truncated Dirichlet process prior of natural parameters of exponential-family model.}

\item{alpha_shape, alpha_rate}{shape and rate parameter of Gamma prior of concentration parameter.}

\item{eta}{the parameters of \code{ergm.terms} and \code{hergm.terms}; the parameters of \code{hergm.terms} must consist of \code{max_number} within-neighborhood parameters and one between-neighborhood parameter.}

\item{eta_mean, eta_sd}{means and standard deviations of Gaussian baseline distribution of Dirichlet process prior of natural parameters.}

\item{eta_mean_mean, eta_mean_sd}{means and standard deviations of Gaussian prior of mean of Gaussian baseline distribution of Dirichlet process prior.}

\item{eta_precision_shape, eta_precision_rate}{shape and rate (inverse scale) parameter of Gamma prior of precision parameter of Gaussian baseline distribution of Dirichlet process prior.}

\item{mean_between}{if \code{simulate == TRUE} and \code{eta == NULL}, then \code{mean_between} specifies the mean-value parameter of edges between blocks.}

\item{indicator}{if the indicators of block memberships of nodes are specified as integers between \code{1} and \code{max_number}, the specified indicators are fixed, which is useful when indicators of block memberhips are observed (e.g., in multilevel networks).}

\item{parallel}{number of computing nodes; if \code{parallel > 1}, \code{hergm} is run on \code{parallel} computing nodes.}

\item{simulate}{if \code{simulate == TRUE}, simulation of networks, otherwise Bayesian inference.}

\item{seeds}{seed of pseudo-random number generator; if \code{parallel > 1}, number of seeds must equal number of computing nodes.}

\item{sample_size}{if \code{simulate == TRUE}, number of network draws, otherwise number of posterior draws; if \code{parallel > 1}, number of draws on each computing node.}

\item{interval}{if \code{simulate == TRUE}, number of proposals between sampled networks.}

\item{burnin}{if \code{simulate == TRUE}, number of burn-in iterations.}

\item{mh.scale}{if \code{simulate == FALSE}, scale factor of candicate-generating distribution of Metropolis-Hastings algorithm.}

\item{variational}{if \code{simulate == FALSE} and \code{variational == TRUE}, variational methods are used to construct the proposal distributions of block memberships of nodes; limited to selected models.}

\item{temperature}{if \code{simulate == FALSE} and \code{variational == TRUE}, minimum and maximum temperature; the temperature is used to melt down the proposal distributions of indicators, which are based on the full conditional distributions of indicators but can have low entropy, resulting in slow mixing of the Markov chain; the temperature is a function of the entropy of the full conditional distributions and is designed to increase the entropy of the proposal distributions, and the minimum and maximum temperature are user-defined lower and upper bounds on the temperature.}

\item{predictions}{if \code{predictions == TRUE} and \code{simulate == FALSE}, returns posterior predictions of statistics in the model.}

\item{posterior.burnin}{number of posterior burn-in iterations; if computing is parallel, \code{posterior.burnin} is applied to the sample generated by each processor;
please note that \code{hergm} returns min(\code{sample_size}, 10000) sample points and the burn-in is applied to the sample of size min(\code{sample_size}, 10000), therefore \code{posterior.burnin} should be smaller than min(\code{sample_size}, 10000).}

\item{posterior.thinning}{if \code{posterior.thinning > 1}, every \code{posterior.thinning}-th sample point is used while all others discarded; if computing is parallel, \code{posterior.thinning} is applied to the sample generated by each processor; please note that \code{hergm} returns min(\code{sample_size}, 10000) sample points and the thinning is applied to the sample of size min(\code{sample_size}, 10000) - \code{posterior.burnin}, therefore \code{posterior.thinning} should be smaller than min(\code{sample_size}, 10000) - \code{posterior.burnin}.}

\item{relabel}{if \code{relabel > 0}, relabel MCMC sample by minimizing the posterior expected loss of Schweinberger and Handcock (2015) (\code{relabel == 1}) or Peng and Carvalho (2015) (\code{relabel == 2}).}

\item{number_runs}{if \code{relabel == 1}, number of runs of relabeling algorithm.}

\item{verbose}{if \code{verbose == -1}, no console output; if \code{verbose == 0}, short console output; if \code{verbose == +1}, long console output.
If, e.g., \code{simulate == FALSE} and \code{verbose == 1}, then \code{hergm} reports the following console output:

Progress: 50.00\% of 1000000

...

means of block parameters:  -0.2838   1.3323

precisions of block parameters:   0.9234   1.4682

block parameters:

  -0.2544  -0.2560  -0.1176  -0.0310  -0.1915  -1.9626

   0.4022   1.8887   1.9719   0.6499   1.7265   0.0000

block indicators: 1 3 1 1 1 1 3 1 1 2 2 2 2 2 1 1 1

block sizes:  10   5   2   0   0

block probabilities: 0.5396 0.2742 0.1419 0.0423 0.0020

block probabilities prior parameter: 0.4256

posterior prediction of statistics:     66    123

where ... indicates additional information about the Markov chain Monte Carlo algorithm that is omitted here.
The console output corresponds to:

- "means of block parameters" correspond to the mean parameters of the Gaussian base distribution of parameters of \code{hergm-terms}.

- "precisions of block parameters" correspond to the precision parameters of the Gaussian base distribution of parameters of \code{hergm-terms}.

- "block parameters" correspond to the parameters of \code{hergm-terms}.

- "block indicators" correspond to the indicators of block memberships of nodes.

- "block sizes" correspond to the block sizes.

- "block probabilities" correspond to the prior probabilities of block memberships of nodes.

- "block probabilities prior parameter" corresponds to the concentration parameter of truncated Dirichlet process prior of parameters of \code{hergm-terms}.

- if \code{predictions == TRUE}, "posterior prediction of statistics" correspond to posterior predictions of sufficient statistics.
}

\item{\dots}{additional arguments, to be passed to lower-level functions in the future.}

}

\value{

The function \code{hergm} returns an object of class \code{hergm} with components:
          
\item{network}{\code{network} is an object of class \code{network} and can be created by calling the function \code{network}.}

\item{formula}{formula of the form \code{network ~ terms}.
          \code{network} is an object of class \code{network} and can be created by calling the function \code{network}.
          Possible terms can be found in \code{ergm.terms} and \code{hergm.terms}.}

\item{n}{number of nodes.}

\item{hyper_prior}{indicator of whether hyper prior has been specified, i.e., whether the parameters \code{alpha}, \code{eta_mean}, and \code{eta_precision} are estimated.}

\item{alpha}{concentration parameter of truncated Dirichlet process prior of parameters of \code{hergm-terms}.}

\item{ergm_theta}{parameters of \code{ergm-terms}.}

\item{eta_mean}{mean parameters of Gaussian base distribution of parameters of \code{hergm-terms}.}

\item{eta_precision}{precision parameters of Gaussian base distribution of parameters of \code{hergm-terms}.}

\item{d1}{total number of parameters of \code{ergm} terms.}

\item{d2}{total number of parameters of \code{hergm} terms.}

\item{hergm_theta}{parameters of \code{hergm-terms}.}

\item{relabeled.hergm_theta}{relabeled parameters of \code{hergm-terms} by using \code{relabel = 1} or \code{relabel = 2}.}

\item{number_fixed}{number of fixed indicators of block memberships of nodes.}

\item{indicator}{indicators of block memberships of nodes.}

\item{relabel}{if \code{relabel > 0}, relabel MCMC sample by minimizing the posterior expected loss of Schweinberger and Handcock (2015) (\code{relabel == 1}) or Peng and Carvalho (2015) (\code{relabel == 2}).}

\item{relabeled.indicator}{relabeled indicators of block memberships of nodes by using \code{relabel = 1} or \code{relabel = 2}.}

\item{size}{the size of the blocks, i.e., the number of nodes of blocks.}

\item{parallel}{number of computing nodes; if \code{parallel > 1}, \code{hergm} is run on \code{parallel} computing nodes.}

\item{p_i_k}{posterior probabilities of block membership of nodes.}

\item{p_k}{probabilities of block memberships of nodes.}

\item{predictions}{if \code{predictions == TRUE} and \code{simulate == FALSE}, returns posterior predictions of statistics in the model.}

\item{simulate}{if \code{simulate == TRUE}, simulation of networks, otherwise Bayesian inference.}

\item{prediction}{posterior predictions of statistics.}

\item{edgelist}{edge list of simulated network.}

\item{sample_size}{if \code{simulate == TRUE}, number of network draws, otherwise number of posterior draws minus number of burn-in iterations; if \code{parallel > 1}, number of draws on each computing node.}

\item{extract}{indicator of whether function \code{hergm.postprocess} has postprocessed the object of class \code{hergm} generated by function \code{hergm} and thus whether the MCMC sample generated by function \code{hergm} has been extracted from the object of class \code{hergm}.}

\item{convergence.diagnostics}{MCMC diagnostics generated by function \code{mcmc.diagnostics}, which in turn relies on function \code{mcgibbsit} in \code{R} package \code{mcgibbsit}; see \code{?mcgibbsit}.} 

\item{verbose}{if \code{verbose == -1}, no console output; if \code{verbose == 0}, short console output; if \code{verbose == +1}, long console output.}

}

\seealso{network, ergm.terms, hergm.terms, hergm.postprocess, mcmc.diagnostics.hergm, summary.hergm, print.hergm, plot.hergm, gof.hergm, simulate.hergm}

\references{

Handcock, M. S. (2003). Assessing degeneracy in statistical models of social networks. Technical report, Center for Statistics and the Social Sciences, University of Washington, Seattle, http://www.csss.washington.edu/Papers.

Holland, P. W. and S. Leinhardt (1981). An exponential family of probability distributions for directed graphs. Journal of the American Statistical Association, Theory \& Methods, 76, 33--65.

Nowicki, K. and T. A. B. Snijders (2001). Estimation and prediction for stochastic blockstructures. Journal of the American Statistical Association, Theory \& Methods, 96, 1077--1087.

Peng, L. and L. Carvalho (2015). Bayesian degree-corrected stochastic block models for community detection. Technical report, Boston University, arXiv:1309.4796v1.

Snijders, T. A. B. and K. Nowicki (1997). Estimation and prediction for stochastic blockmodels for graphs with latent block structure. Journal of Classification 14, 75--100.

Schweinberger, M. (2011). Instability, sensitivity, and degeneracy of discrete exponential families. Journal of the American Statistical Association, Theory & Methods, 106, 1361--1370.

Schweinberger, M. and M. S. Handcock (2015). Local dependence in random graph models: characterization, properties, and statistical Inference. Journal of the Royal Statistical Society, Series B (Statistical Methodology), 7, 647-676.

Schweinberger, M., Petrescu-Prahova, M. and D. Q. Vu (2014). Disaster response on September 11, 2001 through the lens of statistical network analysis. Social Networks, 37, 42--55.

Vu, D. Q., Hunter, D. R. and M. Schweinberger (2013). Model-based clustering of large networks. Annals of Applied Statistics, 7, 1010--1039.

}

\examples{
\dontrun{data(example)
hergm(d ~ edges_i)
hergm(d ~ edges_ij)
hergm(d ~ edges_ij + triangle_ijk)

data(sampson)
hergm(samplike ~ arcs_i + arcs_j)
hergm(samplike ~ edges_ij + mutual_ij)
hergm(samplike ~ edges_ij + mutual_ij + ttriple_ijk)
}
}

