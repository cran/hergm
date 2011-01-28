/***************************************************************************/
/* Copyright 2009 Michael Schweinberger                                    */
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

#include "h_ergm_latent.h"

latentstructure* Initialize_Latentstructure(int number, int n, int minimum_size, int threshold, int d, int *between)
/*
input: maximum number of categories, number of nodes, minimum number of nodes so that structural parameters show up in ergm pmf, number of structural parameters, indicators of wether between-category parameters are restricted to 0
ouput: latent structure
*/
{
  int i, k;
  latentstructure *ls;
  ls = (latentstructure*) calloc(1,sizeof(latentstructure));
  if (ls == NULL) { Rprintf("\n\ncalloc failed: Initialize_Latentstructure, ls\n\n"); exit(1); }
  /* Basics: */
  ls->number = number; /* (Maximum) number of categories */
  ls->n = n; /* Number of nodes */
  /* Law generating structure: */
  ls->p = (double*) calloc(number,sizeof(double)); /* Category-bound probability */
  if (ls->p == NULL) { Rprintf("\n\ncalloc failed: Initialize_Latentstructure, ls->p\n\n"); exit(1); }
  /* Structure: */
  ls->size = (int*) calloc(number,sizeof(int)); /* Category-bound variable: number of nodes belonging to category */
  if (ls->size == NULL) { Rprintf("\n\ncalloc failed: Initialize_Latentstructure, ls->size\n\n"); exit(1); }
  ls->indicator = (int*) calloc(n,sizeof(int)); /* Node-bound variable: category to which node belongs */
  if (ls->indicator == NULL) { Rprintf("\n\ncalloc failed: Initialize_Latentstructure, ls->indicator\n\n"); exit(1); }
  /* Law generating data: */
  ls->minimum_size = minimum_size; /* Minimum number of nodes so that structural parameters show up in PMF */
  ls->threshold = threshold; /* Category-bound PMF tractable as long as number of nodes in category is smaller than threshold */
  ls->d = d; /* Number of category-bound parameters */
  ls->number_between = 0; /* Number of between-category parameters */
  for (i = 0; i < d; i++) 
    {
    ls->number_between = ls->number_between + between[i];
    }
  if (ls->number_between > 0)
    {
    ls->between = (int*) calloc(ls->number_between,sizeof(int)); /* Indicators of whether between-category parameters are restricted to 0 */
    if (ls->between == NULL) { Rprintf("\n\ncalloc failed: Initialize_Latentstructure, ls->between\n\n"); exit(1); }
    k = -1;
    for (i = 0; i < d; i++)
      {
      if (between[i] == 1) 
        {
        k = k + 1;
        ls->between[k] = i;
        } 
      }
    }
  ls->theta = (double**) calloc(d,sizeof(double*));
  if (ls->theta == NULL) { Rprintf("\n\ncalloc failed: Initialize_Latentstructure, ls->theta\n\n"); exit(1); }
  for (i = 0; i < d; i++)
    {
    ls->theta[i] = (double*) calloc(number+1,sizeof(double)); /* Category-bound parameters */
    if (ls->theta[i] == NULL) { Rprintf("\n\ncalloc failed: Initialize_Latentstructure, ls->theta[%i]\n\n",i); exit(1); }
    }
  return ls;
}

ergmstructure* Initialize_Ergm(int terms, int *hierarchical, int d, int d1, int d2, int *structural)
/* 
input: number of ergm terms, indicator of hierarchical model terms, number of parameters, number of non-structural, structural parameters, indicator of structural parameters
output: ergm structure
*/
{
  ergmstructure *ergm;
  ergm = (ergmstructure*) calloc(1,sizeof(ergmstructure));
  if (ergm == NULL) { Rprintf("\n\ncalloc failed: Initialize_Ergm, ergm\n\n"); exit(1); }
  ergm->terms = terms; /* Number of ergm terms */
  ergm->hierarchical = (int*) calloc(terms,sizeof(int)); /* Indicator of hierarchical ergm terms */
  if (ergm->hierarchical == NULL) { Rprintf("\n\ncalloc failed: Initialize_Ergm, ergm->hierarchical\n\n"); exit(1); }
  Set_I_I(terms,ergm->hierarchical,hierarchical);
  ergm->d = d; /* Number of parameters */
  ergm->d1 = d1; /* Number of non-structural parameters */
  ergm->d2 = d2; /* Number of structural parameters */
  ergm->structural = (int*) calloc(d,sizeof(int)); /* Indicator of structural parameters */
  if (ergm->structural == NULL) { Rprintf("\n\ncalloc failed: Initialize_Ergm, ergm->structural\n\n"); exit(1); }
  Set_I_I(d,ergm->structural,structural);
  ergm->theta = (double*) calloc(d1,sizeof(double)); /* Non-structural parameters */
  if (ergm->theta == NULL) { Rprintf("\n\ncalloc failed: Initialize_Ergm, ergm->theta\n\n"); exit(1); }
  return ergm;
}

void Finalize_Latentstructure(latentstructure *ls, int d)
/*
input: maximum number of categories, number of nodes, minimum number of nodes so that structural parameters show up in ergm pmf, number of structural parameters
ouput: latent structure
*/
{
  int i;
  free(ls->p);
  free(ls->size);
  free(ls->indicator);
  if (ls->number_between > 0) free(ls->between);
  for (i = 0; i < d; i++)
    {
    free(ls->theta[i]);
    }
  free(ls->theta);
  free(ls);
}

void Finalize_Ergm(ergmstructure *ergm)
/* 
input: number of ergm terms, indicator of hierarchical model terms, number of parameters, number of non-structural, structural parameters, indicator of structural parameters
output: ergm structure
*/
{
  free(ergm->hierarchical);
  free(ergm->structural);
  free(ergm->theta);
  free(ergm);
}

