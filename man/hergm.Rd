\name{hergm}

\alias{hergm}

\title{Hierarchical exponential-family random graph models with local dependence}

\description{

The function \code{hergm} estimates and simulates three classes of hierarchical exponential-family random graph models:

1. The p_1 model of Holland and Leinhardt (1981) in exponential-family form and extensions by Vu, Hunter, and Schweinberger (2013) and Schweinberger, Petrescu-Prahova, and Vu (2014) to both directed and undirected random graphs with additional model terms, with and without covariates, and with parametric and nonparametric priors (see \code{arcs_i}, \code{arcs_j}, \code{edges_i}, \code{edges_ij}, \code{mutual_i}, \code{mutual_ij}).

2. The stochastic block model of Snijders and Nowicki (1997) and Nowicki and Snijders (2001) in exponential-family form and extensions by Vu, Hunter, and Schweinberger (2013) and Schweinberger, Petrescu-Prahova, and Vu (2014) with additional model terms, with and without covariates, and with parametric and nonparametric priors (see \code{arcs_i}, \code{arcs_j}, \code{edges_i}, \code{edges_ij}, \code{mutual_i}, \code{mutual_ij}).

3. The exponential-family random graph models with local dependence of Schweinberger and Handcock (2015), with and without covariates, and with parametric and nonparametric priors (see \code{arcs_i}, \code{arcs_j}, \code{edges_i}, \code{edges_ij}, \code{mutual_i}, \code{mutual_ij}, \code{twostar_ijk}, \code{triangle_ijk}, \code{ttriple_ijk}, \code{ctriple_ijk}).
The exponential-family random graph models with local dependence replace the long-range dependence of conventional exponential-family random graph models by short-range dependence. 
Therefore, exponential-family random graph models with local dependence replace the strong dependence of conventional exponential-family random graph models by weak dependence, 
reducing the problem of model degeneracy (Handcock, 2003; Schweinberger, 2011) and improving goodness-of-fit (Schweinberger and Handcock, 2015).
In addition, exponential-family random graph models with local dependence satisfy a weak form of self-consistency in the sense that these models are self-consistent under neighborhood sampling (Schweinberger and Handcock, 2015), 
which enables consistent estimation of neighborhood-dependent parameters (Schweinberger and Stewart, 2017; Schweinberger, 2017).

}

\usage{

hergm(formula,
      max_number = 2,
      hierarchical = TRUE,
      parametric = FALSE,
      parameterization = "offset",
      initialize = FALSE,
      initialization_method = 1,
      estimate_parameters = TRUE,
      initial_estimate = NULL,
      n_em_step_max = 100,
      max_iter = 4,
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
      method = "ml",
      seeds = NULL,
      sample_size = NULL,
      sample_size_multiplier_blocks = 20,
      NR_max_iter = 200,
      NR_step_len = NULL,
      NR_step_len_multiplier = 0.2, 
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

\item{hierarchical}{hierarchical prior; if \code{hierarchical = TRUE}, prior is hierarchical (i.e., the means and variances of block parameters are governed by a hyper-prior), otherwise non-hierarchical (i.e., the means and variances of block parameters are fixed).}

\item{parametric}{parametric prior; if \code{parametric = FALSE}, prior is truncated Dirichlet process prior, otherwise parametric Dirichlet prior.}

\item{parameterization}{There are three possible parameterizations of within-block terms when using \code{method == "ml"}.
Please note that between-block terms do not use these parameterizations,
and \code{method == "bayes"} allows the parameters of all within-block terms to vary across blocks and hence does not use them either.
\itemize{
\item \code{standard}: 
The parameters of all within-block terms are constant across blocks.
\item \code{offset}: 
The offset \code{log(n[k])} is subtracted from the parameters of the within-block edge terms and is added to the parameters of the within-block mutual edge terms along the lines of Krivitsky, Handcock, and Morris (2011),
Krivitsky and Kolaczyk (2015),
and Stewart, Schweinberger, Bojanowski, and Morris (2019),
where \code{n[k]} is the number of nodes in block \code{k}.
The parameters of all other within-block terms are constant across blocks.
\item \code{size}: 
The parameters of all within-block terms are multiplied by \code{log(n[k])} along the lines of Babkin et al. (2020),
where \code{n[k]} is the number of nodes in block \code{k}.
}}

\item{initialize}{if \code{initialize = TRUE}, initialize block memberships of nodes.}

\item{initialization_method}{if \code{initialization_method = 1}, block memberships of nodes are initialized by walk trap; if \code{initialization_method = 2}, block memberships of nodes are initalized by spectral clustering.}

\item{estimate_parameters}{if \code{method = "ml"} and \code{estimate_parameters = TRUE}, estimate parameters.}

\item{initial_estimate}{if \code{method = "ml"} and \code{estimate_parameters = TRUE}, specifies starting point.}

\item{n_em_step_max}{if \code{method = "ml"}, maximum number of iterations of Generalized Expectation Maximization algorithm estimating the block structure.}

\item{max_iter}{if \code{method = "ml"}, maximum number of iterations of Monte Carlo maximization algorithm estimating parameters given block structure.}

\item{perturb}{if \code{initialize = TRUE} and \code{perturb = TRUE}, initialize block memberships of nodes by spectral clustering and perturb.}

\item{scaling}{if \code{scaling = TRUE}, use size-dependent parameterizations which ensure that the scaling of between- and within-block terms is consistent with sparse edge terms.}

\item{alpha}{concentration parameter of truncated Dirichlet process prior of natural parameters of exponential-family model.}

\item{alpha_shape, alpha_rate}{shape and rate parameter of Gamma prior of concentration parameter.}

\item{eta}{the parameters of \code{ergm.terms} and \code{hergm.terms}; the parameters of \code{hergm.terms} must consist of \code{max_number} within-block parameters and one between-block parameter.}

\item{eta_mean, eta_sd}{means and standard deviations of Gaussian baseline distribution of Dirichlet process prior of natural parameters.}

\item{eta_mean_mean, eta_mean_sd}{means and standard deviations of Gaussian prior of mean of Gaussian baseline distribution of Dirichlet process prior.}

\item{eta_precision_shape, eta_precision_rate}{shape and rate (inverse scale) parameter of Gamma prior of precision parameter of Gaussian baseline distribution of Dirichlet process prior.}

\item{mean_between}{if \code{simulate = TRUE} and \code{eta = NULL}, then \code{mean_between} specifies the mean-value parameter of edges between blocks.}

\item{indicator}{if the indicators of block memberships of nodes are specified as integers between \code{1} and \code{max_number}, the specified indicators are fixed, which is useful when indicators of block memberhips are observed (e.g., in multilevel networks).}

\item{parallel}{number of computing nodes; if \code{parallel > 1}, \code{hergm} is run on \code{parallel} computing nodes.}

\item{simulate}{if \code{simulate = TRUE}, simulate networks from model, otherwise estimate model given observed network.}

\item{method}{if \code{method = "bayes"}, Bayesian methods along the lines of Schweinberger and Handcock (2015) and Schweinberger and Luna (2018) are used; otherwise, if \code{method = "ml"}, then approximate maximum likelihood methods along the lines of Babkin et al. (2020) are used; note that Bayesian methods are the gold standard but are too time-consuming to be applied to networks with more than 100 nodes, whereas the approximate maximum likelihood methods can be applied to networks with thousands of nodes.}

\item{seeds}{seed of pseudo-random number generator; if \code{parallel > 1}, number of seeds must equal number of computing nodes.}

\item{sample_size}{if \code{simulate = TRUE}, number of network draws, otherwise number of posterior draws; if \code{parallel > 1}, number of draws on each computing node.}

\item{sample_size_multiplier_blocks}{if \code{method = "ml"}, multiplier of the number of network draws from within-block subgraphs; the total number of network draws from within-block subgraphs is \code{sample_size_multiplier_blocks} * number of possible edges of largest within-block subgraph; if \code{sample_size_multiplier_blocks = NULL}, then total number of network draws from within-block subgraphs is \code{sample_size}.}

\item{NR_max_iter}{if \code{method = "ml"}, the maximum number of iterationns to be used in the estimation of parameters.}

\item{NR_step_len}{if \code{method = "ml"}, the step-length to be used for increments in the estimation of parameters. If set to NULL (default), then an adaptive step length procedure is used.}

\item{NR_step_len_multiplier}{if \code{method = "ml"}, multiplier for adjusting the step-length in the estimation procedure after a divergent increment.} 

\item{interval}{if \code{simulate = TRUE}, number of proposals between sampled networks.}

\item{burnin}{if \code{simulate = TRUE}, number of burn-in iterations.}

\item{mh.scale}{if \code{simulate = FALSE}, scale factor of candicate-generating distribution of Metropolis-Hastings algorithm.}

\item{variational}{if \code{simulate = FALSE} and \code{variational = TRUE}, variational methods are used to construct the proposal distributions of block memberships of nodes; limited to selected models.}

\item{temperature}{if \code{simulate = FALSE} and \code{variational = TRUE}, minimum and maximum temperature; the temperature is used to melt down the proposal distributions of indicators, which are based on the full conditional distributions of indicators but can have low entropy, resulting in slow mixing of the Markov chain; the temperature is a function of the entropy of the full conditional distributions and is designed to increase the entropy of the proposal distributions, and the minimum and maximum temperature are user-defined lower and upper bounds on the temperature.}

\item{predictions}{if \code{predictions = TRUE} and \code{simulate = FALSE}, returns posterior predictions of statistics in the model.}

\item{posterior.burnin}{number of posterior burn-in iterations; if computing is parallel, \code{posterior.burnin} is applied to the sample generated by each processor;
please note that \code{hergm} returns min(\code{sample_size}, 10000) sample points and the burn-in is applied to the sample of size min(\code{sample_size}, 10000), therefore \code{posterior.burnin} should be smaller than min(\code{sample_size}, 10000).}

\item{posterior.thinning}{if \code{posterior.thinning > 1}, every \code{posterior.thinning}-th sample point is used while all others discarded; if computing is parallel, \code{posterior.thinning} is applied to the sample generated by each processor; please note that \code{hergm} returns min(\code{sample_size}, 10000) sample points and the thinning is applied to the sample of size min(\code{sample_size}, 10000) - \code{posterior.burnin}, therefore \code{posterior.thinning} should be smaller than min(\code{sample_size}, 10000) - \code{posterior.burnin}.}

\item{relabel}{if \code{relabel > 0}, relabel MCMC sample by minimizing the posterior expected loss of Schweinberger and Handcock (2015) (\code{relabel = 1}) or Peng and Carvalho (2016) (\code{relabel = 2}).}

\item{number_runs}{if \code{relabel = 1}, number of runs of relabeling algorithm.}

\item{verbose}{if \code{verbose = -1}, no console output; if \code{verbose = 0}, short console output; if \code{verbose = +1}, long console output.
If, e.g., \code{simulate = FALSE} and \code{verbose = 1}, then \code{hergm} reports the following console output:

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

- if \code{predictions = TRUE}, "posterior prediction of statistics" correspond to posterior predictions of sufficient statistics.
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

\item{relabel}{if \code{relabel > 0}, relabel MCMC sample by minimizing the posterior expected loss of Schweinberger and Handcock (2015) (\code{relabel = 1}) or Peng and Carvalho (2016) (\code{relabel = 2}).}

\item{relabeled.indicator}{relabeled indicators of block memberships of nodes by using \code{relabel = 1} or \code{relabel = 2}.}

\item{size}{the size of the blocks, i.e., the number of nodes of blocks.}

\item{parallel}{number of computing nodes; if \code{parallel > 1}, \code{hergm} is run on \code{parallel} computing nodes.}

\item{p_i_k}{posterior probabilities of block membership of nodes.}

\item{p_k}{probabilities of block memberships of nodes.}

\item{predictions}{if \code{predictions = TRUE} and \code{simulate = FALSE}, returns posterior predictions of statistics in the model.}

\item{simulate}{if \code{simulate = TRUE}, simulation of networks, otherwise Bayesian inference.}

\item{prediction}{posterior predictions of statistics.}

\item{edgelist}{edge list of simulated network.}

\item{sample_size}{if \code{simulate = TRUE}, number of network draws, otherwise number of posterior draws minus number of burn-in iterations; if \code{parallel > 1}, number of draws on each computing node.}

\item{extract}{indicator of whether function \code{hergm.postprocess} has postprocessed the object of class \code{hergm} generated by function \code{hergm} and thus whether the MCMC sample generated by function \code{hergm} has been extracted from the object of class \code{hergm}.}

\item{verbose}{if \code{verbose = -1}, no console output; if \code{verbose = 0}, short console output; if \code{verbose = +1}, long console output.}

}

\seealso{network, ergm.terms, hergm.terms, hergm.postprocess, summary, print, plot, gof, simulate}

\references{

Babkin, S., Stewart, J., Long, X., and M. Schweinberger (2020). Large-scale estimation of random graph models with local dependence. Computational Statistics and Data Analysis, 152, 1--19.

Cao, M., Chen, Y., Fujimoto, K., and M. Schweinberger (2018). A two-stage working model strategy for network analysis under hierarchical exponential random graph models. Proceedings of the 2018 IEEE/ACM International Conference on Advances in Social Networks Analysis and Mining,
290--298. 

Handcock, M. S. (2003). Assessing degeneracy in statistical models of social networks. Technical report, Center for Statistics and the Social Sciences, University of Washington, Seattle.
http://www.csss.washington.edu/Papers.

Holland, P. W. and S. Leinhardt (1981). An exponential family of probability distributions for directed graphs. Journal of the American Statistical Association, Theory & Methods, 76, 33--65.

Krivitsky, P. N., Handcock, M. S., & Morris, M. (2011). 
Adjusting for network size and composition effects in exponential-family random graph models. 
Statistical Methodology, 8(4), 319-339.

Krivitsky, P.N, and Kolaczyk, E. D. (2015). 
On the question of effective sample size in network modeling: An asymptotic inquiry. 
Statistical science: a review journal of the Institute of Mathematical Statistics, 30(2), 184.

Nowicki, K. and T. A. B. Snijders (2001). Estimation and prediction for stochastic blockstructures. Journal of the American Statistical Association, Theory & Methods, 96, 1077--1087.

Peng, L. and L. Carvalho (2016). Bayesian degree-corrected stochastic block models for community detection. Electronic Journal of Statistics 10, 2746--2779.

Schweinberger, M. (2011). Instability, sensitivity, and degeneracy of discrete exponential families. Journal of the American Statistical Association, Theory & Methods, 106, 1361--1370.

Schweinberger, M. (2020). Consistent structure estimation of exponential-family random graph models with block structure. Bernoulli, 26, 1205--1233.

Schweinberger, M. and M. S. Handcock (2015). Local dependence in random graph models: characterization, properties, and statistical inference. Journal of the Royal Statistical Society, Series B (Statistical Methodology), 7, 647-676.

Schweinberger, M., Krivitsky, P. N., Butts, C.T. and J. Stewart (2020). Exponential-family models of random graphs: Inference in finite, super, and infinite population scenarios. Statistical Science, 35, 627-662.

Schweinberger, M. and P. Luna (2018). HERGM: Hierarchical exponential-family random graph models. Journal of Statistical Software, 85, 1--39.

Schweinberger, M., Petrescu-Prahova, M. and D. Q. Vu (2014). Disaster response on September 11, 2001 through the lens of statistical network analysis. Social Networks, 37, 42--55.

Schweinberger, M. and J. Stewart (2020). Concentration and consistency results for canonical and curved exponential-family random graphs. The Annals of Statistics, 48, 374--396.

Snijders, T. A. B. and K. Nowicki (1997). Estimation and prediction for stochastic blockmodels for graphs with latent block structure. Journal of Classification, 14, 75--100.

Stewart, J., Schweinberger, M., Bojanowski, M., and M. Morris (2019). Multilevel network data facilitate statistical inference for curved {ERGM}s with geometrically weighted terms. Social Networks, 59, 98--119.

Vu, D. Q., Hunter, D. R. and M. Schweinberger (2013). Model-based clustering of large networks. Annals of Applied Statistics, 7, 1010--1039.

}

\examples{

\donttest{
data(example) 
m <- summary(d ~ edges)
} 

}
