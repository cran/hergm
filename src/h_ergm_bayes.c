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

#include "h_ergm_bayes.h"

priorstructure_ls* Initialize_Prior_ls(double shape, double rate)
/*
input: shape, rate (inverse scale) of Gamma prior of clustering parameter
output: prior
*/
{
  priorstructure_ls *prior_ls;
  prior_ls = (priorstructure_ls*) calloc(1,sizeof(priorstructure_ls));
  if (prior_ls == NULL) { Rprintf("\n\ncalloc failed: Initialize_Prior_ls, prior_ls\n\n"); error("Error: out of memory"); }
  prior_ls->alpha_shape = shape; 
  prior_ls->alpha_rate = rate; 
  return prior_ls;
}

priorstructure* Initialize_Priorstructure(int d1, int d2)
/*
input: number of non-hierchical, hierarchical ergm terms
output: prior
*/
{
  int i;
  priorstructure *prior;
  prior = (priorstructure*) calloc(1,sizeof(priorstructure));
  if (prior == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior\n\n"); error("Error: out of memory"); }
  prior->mean2_mean = (double*) calloc(d2,sizeof(double)); 
  if (prior->mean2_mean == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior->mean2_mean\n\n"); error("Error: out of memory"); }
  prior->mean2_precision = (double*) calloc(d2,sizeof(double)); 
  if (prior->mean2_precision == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior->mean2_precision\n\n"); error("Error: out of memory"); }
  prior->mean1 = (double*) calloc(d1,sizeof(double)); 
  if (prior->mean1 == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior->mean1\n\n"); error("Error: out of memory"); }
  prior->mean2 = (double*) calloc(d2,sizeof(double)); 
  if (prior->mean2 == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior->mean2\n\n"); error("Error: out of memory"); }
  prior->cf1 = (double**) calloc(d1,sizeof(double*));
  if (prior->cf1 == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior->cf1\n\n"); error("Error: out of memory"); }
  prior->precision1 = (double**) calloc(d1,sizeof(double*));
  if (prior->precision1 == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior->precision1\n\n"); error("Error: out of memory"); }
  for (i = 0; i < d1; i++)
    {
    prior->cf1[i] = (double*) calloc(d1,sizeof(double));
    if (prior->cf1[i] == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior->cf1[%i]\n\n",i); error("Error: out of memory"); }
    prior->precision1[i] = (double*) calloc(d1,sizeof(double));
    if (prior->precision1[i] == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior->precision1[%i]\n\n",i); error("Error: out of memory"); }
    }
  prior->cf2 = (double**) calloc(d2,sizeof(double*));
  if (prior->cf2 == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior->cf2\n\n"); error("Error: out of memory"); }
  prior->precision2 = (double**) calloc(d2,sizeof(double*));
  if (prior->precision2 == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior->precision2\n\n"); error("Error: out of memory"); }
  for (i = 0; i < d2; i++)
    {  
    prior->cf2[i] = (double*) calloc(d2,sizeof(double)); 
    if (prior->cf2[i] == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior->cf2[%i]\n\n",i); error("Error: out of memory"); }
    prior->precision2[i] = (double*) calloc(d2,sizeof(double)); 
    if (prior->precision2[i] == NULL) { Rprintf("\n\ncalloc failed: Initialize_Priorstructure, prior->precision2[%i]\n\n",i); error("Error: out of memory"); }
    }
  return prior;
}

priorstructure* Initialize_Prior(int d1, int d2, double *mean2_mean, double *mean2_precision, double precision2_shape, double precision2_rate, double *mean1, double *mean2, double *cf1, double *cf2, double *precision1, double *precision2)
/* 
input: number of non-hierarchical, hierarchical ergm terms, R input in the form of vectors and (by vec operator) vectorized matrices
output: prior of non-structural, structural parameters
*/
{
  int i, j, k;
  priorstructure *prior;
  prior = Initialize_Priorstructure(d1,d2);
  Set_D_D(d2,prior->mean2_mean,mean2_mean); 
  Set_D_D(d2,prior->mean2_precision,mean2_precision); 
  prior->precision2_shape = precision2_shape;
  prior->precision2_rate = precision2_rate;
  Set_D_D(d1,prior->mean1,mean1); 
  Set_D_D(d2,prior->mean2,mean2); 
  k = 0;
  for (j = 0; j < d1; j++) /* d1 columns */
    {
    for (i = 0; i < d1; i++) /* d1 rows */
      {
      prior->cf1[i][j] = cf1[k]; 
      prior->precision1[i][j] = precision1[k]; 
      k = k + 1;
      }
    }
  k = 0;
  for (j = 0; j < d2; j++) /* d2 columns */
    {
    for (i = 0; i < d2; i++) /* d2 rows */
      {
      prior->cf2[i][j] = cf2[k]; 
      prior->precision2[i][j] = precision2[k]; 
      k = k + 1;      
      }
    }
  return prior;
}

void Finalize_Prior_ls(priorstructure_ls *prior_ls)
/*
input: shape, rate (inverse scale) of Gamma prior of clustering parameter
output: prior
*/
{
  free(prior_ls);
}

void Finalize_Priorstructure(priorstructure *prior, int d1, int d2)
/*
input: number of non-hierchical, hierarchical ergm terms
output: prior
*/
{
  int i;
  free(prior->mean2_mean);
  free(prior->mean2_precision);
  free(prior->mean1);
  free(prior->mean2);
  for (i = 0; i < d1; i++)
    {
    free(prior->cf1[i]);
    free(prior->precision1[i]);
    }
  free(prior->cf1);
  free(prior->precision1);
  for (i = 0; i < d2; i++)
    {  
    free(prior->cf2[i]);
    free(prior->precision2[i]);
    }
  free(prior->cf2);
  free(prior->precision2);
  free(prior);
}

