/***************************************************************************/
/* Copyright 2009 Nobody                                                   */
/*                                                                         */
/* This file is part of hergm.                                             */
/*                                                                         */
/*    hergm is free software: you can redistribute it and/or modify        */
/*    it under the terms of the GNU General Public License as published by */
/*    the Free Software Foundation, either version 3 of the License, or    */
/*    (at your option) any later version.                                  */
/*                                                                         */
/*    hergm is distributed in the hope that it will be useful,             */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/*    GNU General Public License for more details.                         */
/*                                                                         */
/*    You should have received a copy of the GNU General Public License    */
/*    along with hergm.  If not, see <http://www.gnu.org/licenses/>.       */
/*                                                                         */
/***************************************************************************/

#include "h_ergm_utils.h"

double ln(double x)
{
  double y;
  if (x < epsilon) y = log(epsilon);
  else if (x > maximum) y = log(maximum);
  else y = log(x);
  return y;
}

double e(double x)
{
  double y;
  if (x < log(epsilon)) y = epsilon;
  else if (x > log(maximum)) y = maximum;
  else y = exp(x);
  return y;
}

double S(int d, double *p)
/* 
input: given discrete sample space, number of possible outcomes and probabilities of possible outcomes
output: Shannon entropy of discrete distribution on natural logarithmic scale
*/
{
  int i;
  double entropy;
  entropy = 0.0;
  for (i = 0; i < d; i++)
    {
    entropy = entropy - (p[i] * ln(p[i]));
    }
  return entropy;
}

int Sample_Discrete(double *p)
/*
input: probability vector p
output: value of discrete random variable with pmf p
*/
{
  int i;
  double sum, u;
  u = unif_rand(); /* Sample uniform[0,1] */
  i = 0;
  sum = p[0];
  while (sum < u) /* Interval [0,1] partioned into subintervals; find subinterval in which u falls */
    {
    i = i + 1;
    sum = sum + p[i];
    }
  return i;
}

void Sample_Dirichlet(int d, double alpha, double *p)
/*
input: dimension, parameter
output: probability vector
*/
{
  int i;
  double sum;
  sum = 0.0;
  for (i = 0; i < d; i++)
    { 
    p[i] = rgamma(alpha,1.0);
    sum = sum + p[i];
    }
  for (i = 0; i < d; i++)
    {
    p[i] = p[i] / sum;
    }
}

double* Sample_MVN(int d, double *m, double **C)
/* 
input: dimension d, mean vector m, Cholesky factor C of covariance matrix S = C t(C)
output: random vector x with multivariate normal(m,S) pdf
*/
{
  int i, j;
  double sum, *x, *z;
  x = (double*) calloc(d,sizeof(double));
  if (x == NULL) { Rprintf("\n\ncalloc failed: SampleMVN, x\n\n"); error("Error: out of memory"); }
  z = (double*) calloc(d,sizeof(double));
  if (z == NULL) { Rprintf("\n\ncalloc failed: SampleMVN, z\n\n"); error("Error: out of memory"); }
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
  free(z);
  return x;
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
  y = (double*) calloc(d,sizeof(double));
  if (y == NULL) { Rprintf("\n\ncalloc failed: MVN_PDF, y\n\n"); error("Error: out of memory"); }
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
  free(y);
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
  if (log_ratio < 0.0) 
    {
    u = unif_rand(); /* Sample uniform[0,1] */
    p = e(log_ratio); 
    if (u < p) accept = 1; /* Accept proposal with probability p */
    else accept = 0; /* Reject proposal with probability 1 - p */
    }
  else accept = 1; /* Accept proposal with probability 1 */
  return accept;
}

