#include "h_ergm_utils.h"

int Sample_Discrete(double *p)
/*
input: probability vector p
output: value of discrete random variable with pmf p
*/
{
  int i;
  double sum, u;
  /* GetRNGstate(); */
  u = unif_rand(); /* Sample uniform[0,1] */
  i = 0;
  sum = p[0];
  while (sum < u) /* Interval [0,1] partioned into subintervals; find subinterval in which u falls */
    {
    i = i + 1;
    sum = sum + p[i];
    }
  /* PutRNGstate(); */
  return i;
}

double Multinomial_PMF(int n, int d, int *x, double *p)
/*
input: sample size n, dimension d, random vector x, probability vector p
output: multinomial(n,p) pmf on log scale 
*/
{
  int i;
  double log_pmf;
  log_pmf = Stirling(n); /* Stirling's approximation of n! on log scale */
  for (i = 0; i < d; i++) 
    {
    log_pmf = log_pmf - Stirling(x[i]); 
    }
  for (i = 0; i < d; i++) 
    {
    log_pmf = log_pmf + (x[i] * ln(p[i])); 
    }
  return log_pmf;
}

void Mean_Conditional_MVN(int d_x, int d_y, double *mean_x, double *mean_y, double *x, double **b, double *conditional_mean_y)
/* 
input: dimensions of vectors x, y, mean vectors of marginal (multivariate) normal distributions of vectors x, y, vector x 
output: mean of conditional (multivariate) normal distribution of vector y given vector x
*/
{
  int i, j;
  double sum, *z;
  if (d_x == 0.0) Set_D_D(d_y,conditional_mean_y,mean_y);
  else 
    {
    z = D(d_x);
    for (i = 0; i < d_x; i++)
      {
      z[i] = x[i] - mean_x[i]; /* Center x */
      }
    for (i = 0; i < d_y; i++)
      {
      sum = 0.0;
      for (j = 0; j < d_x; j++)
        {
        sum = sum + (b[i][j] * z[j]);
        }
      conditional_mean_y[i] = mean_y[i] + sum;
      }
    }
}

void Sample_MVN(int d, double *m, double **C, double *x)
/* 
input: dimension d, mean vector m, Cholesky factor C of covariance matrix S = C t(C)
output: random vector x with multivariate normal(m,S) pdf
*/
{
  int i, j;
  double sum, *z;
  /* GetRNGstate(); */
  z = D(d);
  for (i = 0; i < d; i++)
    {
    z[i] = norm_rand(); /* Sample normal(0,1) */
    }
  for (i = 0; i < d; i++)
    {
    sum = 0.0;
    for (j = 0; j < d; j++)
      {
      sum = sum + (C[i][j] * z[j]); /* Shift by m and rescale by S */
      }
    x[i] = m[i] + sum;
    }
  /* PutRNGstate(); */
}

double MVN_PDF(int d, double *x, double *m, double **P)
/* 
input: dimension d, random vector x, mean vector m, precision (inverse covariance) matrix P
output: multivariate normal(m,inverse(P)) kernel on log scale 
*/
{
  int i, j;
  double log_kernel, *y;
  log_kernel = 0.0;
  y = D(d);
  for (i = 0; i < d; i++)
    {
    y[i] = x[i] - m[i]; /* Center x */
    }
  for (i = 0; i < d; i++)
    {
    for (j = 0; j < d; j++)
      {
      log_kernel = log_kernel + (y[i] * P[i][j] * y[j]); /* Quadratic form in y */
      }
    }
  log_kernel = - (log_kernel / 2); /* Log kernel */
  return log_kernel;
}

int MH_Decision(double log_ratio)
/*
input: ratio of pdfs times ratio of proposal pdfs on log scale
output: decision: accept proposal or not
*/
{
  int accept;
  double p, u;
  /* GetRNGstate(); */
  if (log_ratio < 0.0) 
    {
    u = unif_rand(); /* Sample uniform[0,1] */
    p = e(log_ratio); 
    if (u < p) accept = 1; /* Accept proposal with probability p */
    else accept = 0; /* Reject proposal with probability 1 - p */
    }
  else accept = 1; /* Accept proposal with probability 1 */
  /* PutRNGstate(); */
  return accept;
}

