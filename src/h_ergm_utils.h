#include <math.h>
#include <Rmath.h>
#include "h_ergm_basics.h"

int Sample_Discrete(double *p);
/*
input: probability vector p
output: value of discrete random variable with pmf p
*/

double Multinomial_PMF(int n, int d, int *x, double *p);
/*
input: sample size n, dimension d, random vector x, probability vector p
output: multinomial(n,p) pmf on log scale 
*/

void Mean_Conditional_MVN(int d_x, int d_y, double *mean_x, double *mean_y, double *x, double **b, double *conditional_mean_y);
/* 
input: dimensions of vectors x, y, mean vectors of marginal (multivariate) normal distributions of vectors x, y, vector x 
output: mean of conditional (multivariate) normal distribution of vector y given vector x
*/

void Sample_MVN(int d, double *m, double **C, double *x);
/* 
input: dimension d, mean vector m, Cholesky factor C of covariance matrix S = C t(C)
output: random vector x with multivariate normal(m,S) pdf
*/

double MVN_PDF(int d, double *x, double *m, double **P);
/* 
input: dimension d, random vector x, mean vector m, precision (inverse covariance) matrix P
output: multivariate normal(m,inverse(P)) kernel on log scale 
*/

int MH_Decision(double log_ratio); 
/*
input: ratio of pdfs times ratio of proposal pdfs on log scale
output: decision: accept proposal or not
*/

