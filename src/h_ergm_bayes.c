#include "h_ergm_bayes.h"

priorstructure_ls* Initialize_Prior_ls(double shape, double rate)
/*
input: shape, rate (inverse scale) of Gamma prior of clustering parameter
output: prior
*/
{
  priorstructure_ls *prior_ls;
  prior_ls = (priorstructure_ls*) S_alloc(1,sizeof(priorstructure_ls));
  prior_ls->alpha_shape = shape; /* Shape */
  prior_ls->alpha_rate = rate; /* Rate (inverse scale) of Gamma prior */
  return prior_ls;
}

priorstructure* Initialize_Priorstructure(int d1, int d2)
/*
input: number of non-hierchical, hierarchical ergm terms
output: prior
*/
{
  priorstructure *prior;
  prior = (priorstructure*) S_alloc(1,sizeof(priorstructure));
  prior->mean1 = D(d1); /* Mean of marginal Gaussian prior of non-structural parameters */
  prior->mean2 = D(d2); /* Mean of marginal Gaussian prior of structural parameters */
  prior->cf1 = DD(d1,d1); /* Cholesky factor of covariance matrix of conditional Gaussian prior of non-structural parameters */
  prior->precision1 = DD(d1,d1); /* Precision (inverse covariance) matrix of conditional Gaussian prior of non-structural parameters */ 
  prior->b = DD(d2,d1);
  prior->cf2 = DD(d2,d2); /* Cholesky factor of covariance matrix of conditional Gaussian prior of structural parameters */
  prior->precision2 = DD(d2,d2); /* Precision (inverse covariance) matrix of conditional Gaussian prior of structural parameters */ 
  return prior;
}

priorstructure* Initialize_Prior(int d1, int d2, double *mean1, double *mean2, double *b, double *cf1, double *cf2, double *precision1, double *precision2)
/* 
input: number of non-hierarchical, hierarchical ergm terms, R input in the form of vectors and (by vec operator) vectorized matrices
output: prior of non-structural, structural parameters
*/
{
  int i, j, k;
  priorstructure *prior;
  prior = Initialize_Priorstructure(d1,d2);
  Set_D_D(d1,prior->mean1,mean1); /* Mean of marginal Gaussian prior of non-structural parameters  */
  Set_D_D(d2,prior->mean2,mean2); /* Mean of marginal Gaussian prior of structural parameters */
  k = 0;
  for (j = 0; j < d1; j++) /* d1 columns */
    {
    for (i = 0; i < d1; i++) /* d1 rows */
      {
      prior->cf1[i][j] = cf1[k]; /* Cholesky factor of covariance matrix of conditional Gaussian prior of non-structural parameters */
      prior->precision1[i][j] = precision1[k]; /* Precision (inverse covariance) matrix of conditional Gaussian prior of non-structural parameters */
      k = k + 1;
      }
    }
  k = 0;
  for (j = 0; j < d1; j++) /* d1 columns */
    {
    for (i = 0; i < d2; i++) /* d2 rows */
      {
      prior->b[i][j] = b[k];
      k = k + 1;      
      }
    }
  k = 0;
  for (j = 0; j < d2; j++) /* d2 columns */
    {
    for (i = 0; i < d2; i++) /* d2 rows */
      {
      prior->cf2[i][j] = cf2[k]; /* Cholesky factor of covariance matrix of conditional Gaussian prior of structural parameters */
      prior->precision2[i][j] = precision2[k]; /* Precision (inverse covariance) matrix of conditional Gaussian prior of structural parameters */
      k = k + 1;      
      }
    }
  return prior;
}

