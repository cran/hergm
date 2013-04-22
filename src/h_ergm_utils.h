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

#include "h_ergm_basics.h"

double epsilon, maximum;

double ln(double x);

double e(double x);

double S(int d, double *p);
/* 
input: given discrete sample space, number of possible outcomes and probabilities of possible outcomes
output: Shannon entropy of discrete distribution on natural logarithmic scale
*/

int Sample_Discrete(double *p);
/*
input: probability vector p
output: value of discrete random variable with pmf p
*/

void Sample_Dirichlet(int d, double alpha, double *p);
/*
input: dimension, parameter
output: probability vector
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

