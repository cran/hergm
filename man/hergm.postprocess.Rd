\name{hergm.postprocess}

\alias{hergm.postprocess}

\alias{postprocess.hergm}

\title{Postprocess object of class \code{hergm}}

\description{

The function \code{hergm.postprocess} postprocesses an object of class \code{hergm}.
Please note that the function \code{hergm} calls the function \code{hergm.postprocess} with \code{relabel = 0} by default or with other values of \code{relabel} specified by the user,
therefore users do not need to call the function \code{hergm.postprocess} unless it is desired to postprocess an object of class \code{hergm} with a value of \code{relabel} that was not used by function \code{hergm}.

If \code{hergm.postprocess} is called with \code{relabel > 0},
it solves the so-called label-switching problem.
The label-switching problem is rooted in the invariance of the likelihood function to permutations of the labels of blocks, and implies that raw MCMC samples from the posterior cannot be used to infer to block-dependent entities.
The label-switching problem can be solved in a Bayesian decision-theoretic framework: by choosing a loss function and minimizing the posterior expected loss.
Two loss functions are implemented in \code{hergm.postprocess}, the loss function of Schweinberger and Handcock (2015) (\code{relabel == 1}) and the loss function of Peng and Carvalho (2016) (\code{relabel == 2}).
The first loss function seems to be superior in terms of the reported clustering probabilities, but is more expensive in terms of computing time.
A rule of thumb is to use the first loss function when \code{max_number < 15} and use the second loss function otherwise.

}

\usage{

hergm.postprocess(object,
             burnin = 2000, 
             thinning = 1,
             relabel = 1,
             number_runs = 1,
             \dots)

}

\arguments{

\item{object}{object of class \code{hergm}; objects of class \code{hergm} can be generated by function \code{hergm}.}

\item{burnin}{number of posterior burn-in iterations; if computing is parallel, \code{burnin} is applied to the sample generated by each processor;
please note that \code{hergm} returns min(\code{sample_size}, 10000) sample points and the burn-in is applied to the sample of size min(\code{sample_size}, 10000), therefore \code{burnin} should be smaller than min(\code{sample_size}, 10000).}

\item{thinning}{if \code{thinning > 1}, every \code{thinning}-th sample point is used while all others discarded; if computing is parallel, \code{thinning} is applied to the sample generated by each processor; please note that \code{hergm} returns min(\code{sample_size}, 10000) sample points and the thinning is applied to the sample of size min(\code{sample_size}, 10000) - \code{burnin}, therefore \code{thinning} should be smaller than min(\code{sample_size}, 10000) - \code{burnin}.}

\item{relabel}{if \code{relabel > 0}, relabel MCMC sample by minimizing the posterior expected loss of Schweinberger and Handcock (2015) (\code{relabel == 1}) or Peng and Carvalho (2016) (\code{relabel == 2}).}

\item{number_runs}{if \code{relabel == 1}, number of runs of relabeling algorithm.}

\item{\dots}{additional arguments, to be passed to lower-level functions in the future.}

}

\value{

\item{ergm_theta}{parameters of \code{ergm-terms}.}

\item{alpha}{concentration parameter of truncated Dirichlet process prior of parameters of \code{hergm-terms}.} 

\item{eta_mean}{mean parameters of Gaussian base distribution of parameters of \code{hergm-terms}.}

\item{eta_precision}{precision parameters of Gaussian base distribution of parameters of \code{hergm-terms}.}

\item{hergm_theta}{parameters of \code{hergm-terms}.}

\item{loss}{if \code{relabel == TRUE}, local minimum of loss function.}

\item{p_k}{probabilities of block memberships of nodes.}

\item{indicator}{indicators of block memberships of nodes.}

\item{p_i_k}{posterior probabilities of block memberships of nodes.}

\item{prediction}{posterior predictions of statistics.}

}

\seealso{hergm}

\references{

Peng, L. and L. Carvalho (2016). Bayesian degree-corrected stochastic block models for community detection. Electronic Journal of Statistics 10, 2746--2779.

Schweinberger, M. and M. S. Handcock (2015). Local dependence in random graph models: characterization, properties, and statistical Inference. Journal of the Royal Statistical Society, Series B (Statistical Methodology), 7, 647-676.

Schweinberger, M. and P. Luna (2018). HERGM: Hierarchical exponential-family random graph models. Journal of Statistical Software, 85, 1--39.

}
