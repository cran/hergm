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
Structure and structural parameters:
- structure: number of categories; each node is member of one and only one category; number of categories and category memberships are latent
- structural parameters: parameters of law generating structure and parameters of law generating data
*/
{
  /* Basics: */
  int number; /* (Maximum) number of categories */
  int n; /* Number of nodes */
  /* Law generating structure: */
  double alpha; /* Shape parameter of beta distribution */
  double *p; /* Category-bound probability */
  /* Structure: */
  int *size; /* Category-bound variable: number of nodes belonging to category */
  int *indicator; /* Node-bound variable: category to which node belongs */
  int *fixed; /* Node-bound variable: indicator of whether indicator is fixed */
  int number_fixed; /* Number of fixed indicators */
  /* Law generating data: */
  int minimum_size; /* Minimum number of nodes so that structural parameters show up in PMF */
  int threshold; /* Category-bound within-block PMF tractable as long as numbers of nodes smaller than threshold */
  int d; /* Number of category-bound parameters */
  int number_between; /* Number of unrestricted between-category parameters */
  int *between; /* Indicators of whether between-category parameters are unrestricted */
  double *scaling; /* Orders of statistics, e.g., number of edges = 2, number of triangles = 3 */
  double **theta; /* Category-bound parameters */
} 
  latentstructure; 

typedef struct
/*
Non-structural parameters of law generating data
*/
{
  int terms; /* Number of ergm terms */
  int *hierarchical; /* Indicator of hierarchical ergm terms */
  int d; /* Number of parameters */
  int d1; /* Number of non-structural parameters */
  int d2; /* Number of structural parameters */
  int *structural; /* Indicator of structural parameters */
  double *theta; /* Non-structural parameters */ 
}
  ergmstructure;

latentstructure* Initialize_Latentstructure(int number, int n, int *indicator, int minimum_size, int threshold, int d, int *between, double *scaling);
/*
input: maximum number of categories, number of nodes, minimum number of nodes so that structural parameters show up in ergm pmf, number of structural parameters, indicators of wether between-category parameters are restricted to 0
ouput: latent structure
*/

ergmstructure* Initialize_Ergm(int terms, int *hierarchical, int d, int d1, int d2, int *structural);
/* 
input: number of ergm terms, indicator of hierarchical model terms, number of parameters, number of non-structural, structural parameters, indicator of structural parameters
output: ergm structure
*/

void Finalize_Latentstructure(latentstructure *ls, int d);
/*
input: maximum number of categories, number of nodes, minimum number of nodes so that structural parameters show up in ergm pmf, number of structural parameters
ouput: latent structure
*/

void Finalize_Ergm(ergmstructure *ergm);
/* 
input: number of ergm terms, indicator of hierarchical model terms, number of parameters, number of non-structural, structural parameters, indicator of structural parameters
output: ergm structure
*/


