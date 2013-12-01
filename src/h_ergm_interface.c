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

double* Extract_Input_Blocks(int terms, int *hierarchical, int max_number, int n, int *indicator, double *input, int *block, double **theta)
/*
input: number of terms, indicators of hierarchical terms, maximum number of blocks, number of nodes, indicators of block memberships, input, number and labels of included blocks, parameters
output: input list of included blocks, which is a subvector of input
*/
{
  int h, i, included, *indicator_block, k, k_block, l, number, number_blocks, number_nodes;
  double *input_block;
  indicator_block = (int*) calloc(n,sizeof(int));
  if (indicator_block == NULL) { Rprintf("\n\ncalloc failed: Set_Input_Blocks, indicator_block\n\n"); error("Error: out of memory"); }
  number_blocks = block[0]; /* Number of blocks included; the labels of included blocks are stored in block[1], ..., block[number_blocks] */
  number_nodes = 0; /* Number of nodes which are members of included blocks */
  for (i = 0; i < n; i++)
    {
    included = 0;
    k = 0;
    while ((k < number_blocks) && (included == 0)) /* Check whether node is member of included blocks */
      {
      k = k + 1;
      if (indicator[i] == block[k]) 
        {
        included = 1;
        number_nodes = number_nodes + 1;
        indicator_block[number_nodes-1] = block[k]; /* Conclusion: node is member of included blocks */
        }
      }
    }
  input_block = NULL;
  h = -1; /* Hierarchical ergm term h */
  k = -1; /* Input parameter k */
  k_block = -1; /* Block input parameter k_block */
  for (i = 0; i < terms; i++) /* For given ergm term... */
    {
    if (hierarchical[i] == 0) /* ...if non-hierarchical, copy input parameters */
      {
      for (l = 0; l < 3; l++) /* Copy elements 1, 2, 3 */
        {
        k = k + 1;
        k_block = k_block + 1;
        input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
        input_block[k_block] = input[k]; 
        }
      number = trunc(input[k]); /* Element 3: total number of input parameters */
      for (l = 0; l < number; l++) /* If number > 0, copy elements 4, ..., 3 + number */
        {
        k = k + 1;
        k_block = k_block + 1;
        input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
        input_block[k_block] = input[k]; 
        } 
      }
    else /* ...if hierarchical, copy input parameters and set block indicators */ 
      {
      h = h + 1; /* Hierarchical ergm term h */
      for (l = 0; l < 2; l++) /* Copy elements 1, 2 */
        {
        k = k + 1; 
        k_block = k_block + 1;
        input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
        input_block[k_block] = input[k];
        }
      k = k + 1; 
      k_block = k_block + 1;
      input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
      input_block[k_block] = 1.0 + number_nodes + (max_number + 1.0); /* Element 3: total number of input parameters: (maximum) number of categories, number_nodes block indicators, (max_number + 1) within- and between-block parameters */
      k = k + 1; 
      k_block = k_block + 1;
      input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
      input_block[k_block] = max_number; /* Element 4: (maximum) number of blocks */
      k = k + n; /* Please note: k must be incremented by n rather than number_nodes */
      for (l = 0; l < number_nodes; l++) /* Set elements 4 + 1 + number_nodes: block indicators */
        {
        k_block = k_block + 1; 
        input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
        input_block[k_block] = indicator_block[l]; /* Block indicator of member l of included blocks */
        }
      for (l = 0; l < max_number + 1; l++) /* Copy elements 4 + number_nodes + max_number + 1: within- and between-block block parameters */
        {
        k = k + 1;
        k_block = k_block + 1;
        input_block = (double*) realloc(input_block,(k_block+1)*sizeof(double));
        input_block[k_block] = theta[h][l];
        }
      }
    }
  free(indicator_block);
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

