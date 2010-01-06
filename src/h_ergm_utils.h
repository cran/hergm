#include "h_ergm_basics.h"

double epsilon, maximum;

double ln(double x);

double e(double x);

int Sample_Discrete(double *p);
/*
input: probability vector p
output: value of discrete random variable with pmf p
*/

double* Sample_MVN(int d, double *m, double **C);
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

