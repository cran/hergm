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
  double *mean2_mean; /* Mean of Gaussian prior of mean of baseline distribution of Dirichlet / stick-breaking prior */
  double *mean2_precision; /* Precision of Gaussian prior of mean of baseline distribution of Dirichlet / stick-breaking prior */
  double precision2_shape; /* Shape of Gamma prior of precision of Gaussian baseline distribution of Dirichlet / stick-breaking prior */
  double precision2_rate; /* Rate (inverse scale) of Gamma prior of precision of Gaussian baseline distribution of Dirichlet / stick-breaking prior */
  double *mean1; /* Mean of Gaussian prior */
  double *mean2; /* Mean of Gaussian baseline distribution of Dirichlet / stick-breaking prior */
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

priorstructure* Initialize_Prior(int d1, int d2, double *mean2_mean, double *mean2_precision, double precision2_shape, double precision2_rate, double *mean1, double *mean2, double *cf1, double *cf2, double *precision1, double *precision2);
/* 
input: number of non-hierarchical, hierarchical ergm terms, R input in the form of vectors and (by vec operator) vectorized matrices
output: prior of non-structural, structural parameters
*/

void Finalize_Prior_ls(priorstructure_ls *prior_ls);
/*
input: shape, rate (inverse scale) of Gamma prior of clustering parameter
output: prior
*/

void Finalize_Priorstructure(priorstructure *prior, int d1, int d2);
/*
input: number of non-hierchical, hierarchical ergm terms
output: prior
*/

