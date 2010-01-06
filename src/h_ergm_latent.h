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
  /* Law generating data: */
  int threshold; /* Minimum number of nodes so that structural parameters show up in ergm pmf */
  int d; /* Number of category-bound parameters */
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

latentstructure* Initialize_Latentstructure(int number, int n, int threshold, int d);
/*
input: maximum number of categories, number of nodes, minimum number of nodes so that structural parameters show up in ergm pmf, number of structural parameters
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


