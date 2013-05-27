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

#include "h_ergm_interface.h"

int Number_Input(int terms, double *input)
/*
input: number of ergm terms, input parameters
output: number of input parameters
*/
{
  int i, k, number;
  k = -1; 
  for (i = 0; i < terms; i++) 
    {                        
    k = k + 3; 
    number = trunc(input[k]); /* Element 3: total number of input parameters */
    k = k + number;
    }
  k = k + 1; /* Since input starts at 0, number of elements is not k but k + 1 */
  return k;
}

void Set_Input(int terms, int *hierarchical, int max_number, int n, int *indicator, double **theta, double *input)
/*
input: number of ergm terms, indicator of hierarchical ergm terms, number of nodes, node-bound category indicator, category-bound parameter
output: input parameters
*/
{
  int h, i, j, k, number;
  h = -1; /* Hierarchical ergm term h */
  k = -1; /* Input parameter k */
  for (i = 0; i < terms; i++) /* For given ergm term... */
    {
    if (hierarchical[i] == 0) /* ...if non-hierarchical, go to following ergm term */
      {
      k = k + 3; /* Elements 1, 2, 3 unchanged */
      number = trunc(input[k]); /* Element 3: total number of input parameters */
      k = k + number; /* If number > 0, elements 4, ..., 3 + number unchanged */ 
      }
    else /* ...if hierarchical, set input parameters */ 
      {
      h = h + 1; /* Hierarchical ergm term h */
      k = k + 1; 
      input[k] = 0.0; /* Element 1 */
      k = k + 1;
      input[k] = 1.0; /* Element 2: one change statistic */
      k = k + 1;
      input[k] = 1.0 + n + (max_number + 1.0); /* Element 3: total number of input parameters: (maximum) number of categories, n node-bound category indicators, (max_number + 1) category-bound parameters */
      k = k + 1; 
      input[k] = max_number; /* Elements 4: (maximum) number of categories */
      for (j = 0; j < n; j++) /* Elements 4 + n: category indicators */
        {
        k = k + 1;
        input[k] = indicator[j]; 
        }
      for (j = 0; j < max_number; j++) /* Elements 4 + n + max_number: within-category parameters */
        {
        k = k + 1; 
        input[k] = theta[h][j]; 
        }
      k = k + 1;
      input[k] = theta[h][max_number]; /* Element 4 + n + max_number + 1: between-category parameter */       
      }
    }
}

double* Set_Input_Block(int terms, int *hierarchical, int max_number, int n, int n_block, double *theta_block, double *input)
/*
input: number of ergm terms, indicator of hierarchical ergm terms, number of nodes, node-bound category indicator, category-bound parameter
output: input parameters
*/
{
  int k_block, h, i, j, k, l, number, max_number_block;
  double *input_block;
  input_block = NULL;
  max_number_block = 1;
  h = -1; /* Hierarchical ergm term h */
  k = -1; /* Input parameter k */
  k_block = -1;
  for (i = 0; i < terms; i++) /* For given ergm term... */
    {
    if (hierarchical[i] == 0) /* ...if non-hierarchical, go to following ergm term */
      {
      k = k + 1; /* Copy element 1 */
      k_block = k_block + 1;
      input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
      input_block[k_block] = input[k]; 
      k = k + 1; /* Copy element 2 */
      k_block = k_block + 1;
      input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
      input_block[k_block] = input[k]; 
      k = k + 1; /* Copy element 3 */
      k_block = k_block + 1;
      input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
      input_block[k_block] = input[k]; 
      number = trunc(input[k]); /* Element 3: total number of input parameters */
      for (l = 0; l < number; l++) /* If number > 0, copy elements 4, ..., 3 + number */
        {
        k = k + 1;
        k_block = k_block + 1;
        input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
        input_block[k_block] = input[k]; 
        } 
      }
    else /* ...if hierarchical, set input parameters */ 
      {
      h = h + 1; /* Hierarchical ergm term h */
      k = k + 1; 
      k_block = k_block + 1;
      input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
      input_block[k_block] = 0.0; /* Element 1 */
      k = k + 1;
      k_block = k_block + 1;
      input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
      input_block[k_block] = 1.0; /* Element 2: one change statistic */
      k = k + 1;
      k_block = k_block + 1;
      input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
      input_block[k_block] = 1.0 + n_block + (max_number_block + 1.0); /* Element 3: total number of input parameters: (maximum) number of categories, n node-bound category indicators, (max_number + 1) category-bound parameters */
      k = k + 1;
      k_block = k_block + 1;
      input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
      input_block[k_block] = max_number_block; /* Elements 4: (maximum) number of categories */
      k = k + n;
      for (l = 0; l < n_block; l++)
        {
        k_block = k_block + 1; /* Elements 4 + 1 + max_number: category indicators */
        input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
        input_block[k_block] = 0.0;
        }
      k = k + max_number;  
      k_block = k_block + max_number_block; /* Elements 4 + n + 1: within-category parameter */
      input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
      input_block[k_block] = theta_block[h]; 
      k = k + 1;
      k_block = k_block + 1; /* Elements 4 + n + 1 + 1: between-category parameter */
      input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
      input_block[k_block] = 0.0; 
      }
    }
  return input_block;
}

void Set_Input_Indicator(int terms, int *hierarchical, int max_number, int n, int node, int node_indicator, double *input)
/*
input: number of ergm terms, indicator of hierarchical ergm terms, number of nodes, node-bound category indicator, category-bound parameter
output: input parameters
*/
{
  int i, k, index, number;
  k = -1; /* Input parameter k */
  for (i = 0; i < terms; i++) /* For given ergm term... */
    {
    if (hierarchical[i] == 0) /* ...if non-hierarchical, go to following ergm term */
      {
      k = k + 3; /* Elements 1, 2, 3 unchanged */
      number = trunc(input[k]); /* Element 3: total number of input parameters */
      k = k + number; /* If number > 0, elements 4, ..., 3 + number unchanged */ 
      }
    else /* ...if hierarchical, set input parameters */ 
      {
      k = k + 4; /* Elements 1..4 contain general information */
      index = k + (node + 1); 
      input[index] = node_indicator; /* Elements 4 + (node + 1) contains indicator of node */ 
      k = k + n + (max_number + 1); /* Elements 5..4 + n contain category indicators, elements 4 + n + max_number + 1 contain parameters */
      }
    }
}

double* Get_Parameter(int d, int *structural, double *theta)
/*
input: number of ergm terms, indicator of structural parameters
output: parameter
*/
{
  int i, k;
  double *parameter;
  parameter = (double*) calloc(d,sizeof(double)); 
  if (parameter == NULL) { Rprintf("\n\ncalloc failed: Get_Parameter, parameter\n\n"); error("Error: out of memory"); }
  k = -1;
  for (i = 0; i < d; i++)
    {
    if (structural[i] == 0) /* Non-structural parameter */
      {
      k = k + 1;
      parameter[i] = theta[k];
      }
    else /* Structural parameter */
      {
      parameter[i] = 1.0; /* Structural parameters enter ergm pmf through input parameters of "change statistics" */
      }
    }
  return parameter;
}

void Set_Parameter(int d, int *structural, double *theta, double *parameter)
/*
input: number of ergm terms, indicator of structural parameters
output: parameter
*/
{
  int i, k;
  k = -1;
  for (i = 0; i < d; i++)
    {
    if (structural[i] == 0) /* Non-structural parameter */
      {
      k = k + 1;
      parameter[i] = theta[k];
      }
    else /* Structural parameter */
      {
      parameter[i] = 1.0; /* Structural parameters enter ergm pmf through input parameters of "change statistics" */
      }
    }
}

double Minus_Energy(int d, double *input, double *parameter, 
                       int *heads, int *tails, int *nedges, 
		       int *n, int *directed,  int *bipartite,
		       int *nterms, char **funnames,
		       char **sonames,
                       double *statistic)
/*
input: number of parameters, input parameters, parameters
output: statistic, inner product <parameter, statistic>
*/
{
  int i;
  double sum;
  /*
  Rprintf("\nMinus_Energy: number of edges = %i",*nedges);
  */
  for (i = 0; i < d; i++) /* Statistic must be null */
    {
    statistic[i] = 0.0;
    }
  int timings = 0, time = 0, lasttoggle = 0;
  /*
  Rprintf("\n\nh_ergm_interface.c:");
  Rprintf("\ntimings=%p",&timings);
  Rprintf("\ntimings=%i",timings);
  */
  network_stats_wrapper(tails,heads,&timings,&time,&lasttoggle,nedges,n,directed,bipartite,nterms,funnames,sonames,input,statistic); /* Compute statistic given input */
  sum = 0.0;
  for (i = 0; i < d; i++)
    {
    sum = sum + (parameter[i] * statistic[i]);
    /*
    Rprintf("\nparameter[%i] = %-8.4f, statistic[%i] = %-8.4f, product = %-8.4f", i, parameter[i], i, statistic[i], parameter[i] * statistic[i]);
    */
    }
  return sum;
} 

