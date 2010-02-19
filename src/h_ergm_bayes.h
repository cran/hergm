#include "h_ergm_basics.h"

typedef struct
/*
Prior of clustering parameter: Gamma
*/
{
  double alpha_shape; /* Shape */
  double alpha_rate; /* Rate (inverse scale) */
}
  priorstructure_ls;

typedef struct
/*
Joint prior of parameters generating data: Gaussian
*/
{
  double *mean1; /* Mean of marginal Gaussian prior of non-hierarchical ergm terms */
  double *mean2; /* Mean of marginal Gaussian prior of hierarchical ergm terms */
  double **b; 
  double **cf1; /* Cholesky factor of covariance matrix of conditional Gaussian prior of non-structural parameters */
  double **cf2; /* Cholesky factor of covariance matrix of conditional Gaussian prior of structural parameters */
  double **precision1; /* Precision (inverse covariance) matrix of conditional Gaussian prior of non-structural parameters */
  double **precision2; /* Precision (inverse covariance) matrix of conditional Gaussian prior of structural parameters */
} 
  priorstructure;

typedef struct
/*
MCMC sample: structure
*/
{
  /* Law generating structure: */
  double *accept_theta; /* Acceptance rate of Metropolis-Hastings algorithm for updating structural parameters */
  double accept_indicator; /* Acceptance rate of Metropolis-Hastings algorithm for updating indicators */
  double *alpha; /* Shape parameter of beta distribution */
  double **p; /* Category-bound probability */
  /* Structure: */
  int **size; /* Category-bound variable: number of nodes belonging to category */
  int **indicator; /* Node-bound variable: category to which node belongs */
  /* Law generating data: */
  double ***theta; /* Category-bound parameter */
}
  mcmcstructure_ls;

typedef struct
/*
MCMC sample: non-structural parameters
*/
{
  double accept; /* Acceptance rate of Metropolis-Hastings algorithm for updating non-structural parameters */
  double **theta; /* Non-structural parameters */
}
  mcmcstructure;

priorstructure_ls* Initialize_Prior_ls(double shape, double rate);
/*
input: shape, rate (inverse scale) of Gamma prior of clustering parameter
output: prior
*/

priorstructure* Initialize_Priorstructure(int d1, int d2);
/*
input: number of non-hierchical, hierarchical ergm terms
output: prior
*/

priorstructure* Initialize_Prior(int d1, int d2, double *mean1, double *mean2, double *b, double *cf1, double *cf2, double *precision1, double *precision2);
/* 
input: number of non-hierarchical, hierarchical ergm terms, R input in the form of vectors and (by vec operator) vectorized matrices
output: prior of non-structural, structural parameters
*/

