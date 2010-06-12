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
  if (parameter == NULL) { Rprintf("\n\ncalloc failed: Get_Parameter, parameter\n\n"); exit(1); }
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

void Sample_Graph(int number, int n, int ls_d, int terms, int *hierarchical, int ergm_d, double *statistic,
                          int *heads, int *tails, int *dnedges,
                          int *maxpossibleedges,
                          int *dn, int *directed, int *bipartite, 
                          int *nterms, char **funnames,
                          char **sonames, 
                          char **MHproposaltype, char **MHproposalpackage,
                          double *input, double *theta, int *samplesize, 
                          double *sample, int *burnin, int *interval,  
                          int *newnetworkheads, 
                          int *newnetworktails, 
                          int *fVerbose, 
                          int *attribs, int *maxout, int *maxin, int *minout,
                          int *minin, int *condAllDegExact, int *attriblength, 
                          int *maxedges,
                          int *mheads, int *mtails, int *mdnedges)
/*
input: (maximum) number of categories, number of nodes, number of structural parameters, number of parameters
output: one sample from posterior predictive distribution
*/
{
  int *h, *t, i, *indicator, k, *nedges, s;
  double **parameter;
  s = 1; /* Sample one graph from posterior predictive distribution */
  for (i = 0; i < ergm_d; i++)
    {
    sample[i] = 0.0;
    }
  MCMC_wrapper(heads,tails,dnedges, /* Sample one graph from posterior predictive distribution given input and theta */
               maxpossibleedges,
               dn,directed,bipartite,
               nterms,funnames,
               sonames,
               MHproposaltype,MHproposalpackage,
               input,theta,&s,
               sample,burnin,interval, 
               newnetworkheads,
               newnetworktails,
               fVerbose,
               attribs,maxout,maxin,minout,
               minin,condAllDegExact,attriblength,
               maxedges,
               mheads,mtails,mdnedges);
  indicator = (int*) calloc(n,sizeof(int));
  if (indicator == NULL) { Rprintf("\n\ncalloc failed: Sample_Graph, indicator\n\n"); exit(1); }
  for (i = 0; i < n; i++) /* Identical indicators */
    {
    indicator[i] = 1;
    }
  parameter = (double**) calloc(ls_d,sizeof(double*));
  if (parameter == NULL) { Rprintf("\n\ncalloc failed: Sample_Graph, parameter\n\n"); exit(1); }
  for (i = 0; i < ls_d; i++)
    {
    parameter[i] = (double*) calloc(number+1,sizeof(double));
    if (parameter[i] == NULL) { Rprintf("\n\ncalloc failed: Sample_Graph, parameter[%i]\n\n",i); exit(1); }
    }
  for (i = 0; i < ls_d; i++) /* Identical structural parameters, so that structural function of graph, structural parameters reduces to corresponding non-structural function of graph */
    {
    for (k = 0; k < (number + 1); k++)
      {
      parameter[i][k] = 1.0; 
      }
    }
  Set_Input(terms,hierarchical,number,n,indicator,parameter,input);
  nedges = newnetworkheads; /* First element of newnetworkheads = newnetworkheads[0] is number of edges */
  h = (int*) calloc(*nedges,sizeof(int));
  if (h == NULL) { Rprintf("\n\ncalloc failed: Sample_Graph, h\n\n"); exit(1); }
  t = (int*) calloc(*nedges,sizeof(int));
  if (t == NULL) { Rprintf("\n\ncalloc failed: Sample_Graph, t\n\n"); exit(1); }
  for (i = 0; i < *nedges; i++) /* Since first element of newnetworkheads and newnetworktails is number of edges, heads and tails must be extracted */
    {
    h[i] = newnetworkheads[i+1];
    t[i] = newnetworktails[i+1];
    }
  network_stats_wrapper(h,t,nedges,dn,directed,bipartite,nterms,funnames,sonames,input,statistic); /* Compute non-structural function of graph */
  free(h);
  free(t);
  free(indicator);
  for (i = 0; i < ls_d; i++)
    {
    free(parameter[i]);
    }
  free(parameter);
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
  network_stats_wrapper(heads,tails,nedges,n,directed,bipartite,nterms,funnames,sonames,input,statistic); /* Compute statistic given input */
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

double Ratio_Partition_Functions(int s, int d, double sum_observed, double *statistic_generating, double *statistic, double *theta_generating, double *theta)
/*
input: sample size, dimension, difference of inner products under alternative, data-generating parameter for observed graph, value of statistic under data-generating, alternative parameter, value of data-generating, alternative parameter
output: ratio of partition functions of ergms under alternative and data-generating parameter on log scale
*/
{ 
  int i, j, k;
  double ratio_n_const, log, log_generating, log_ratio_n_const, log_ratio, moment1, moment2, sum, variance;
  /* Ratio of partition functions is expectation of exponential functions:
  - estimated by: MCMC sample average
  - computations stabilized by: log-normal approximation of log expectation = log ratio of partition constants
  */
  moment1 = 0.0;
  moment2 = 0.0;
  k = 0; 
  for (i = 0; i < s; i++) /* Sample i */
    {
    sum = 0.0; 
    for (j = 0; j < d; j++) /* Given sample i... */
      {
      log = theta[j] * statistic[k]; /* ...compute inner product; note indices j, k */
      log_generating = theta_generating[j] * statistic_generating[k]; /* ...compute inner product; note indices j, k */
      sum = sum + (log - log_generating); /* ...compute difference in inner products */
      k = k + 1; /* Get following elements of statistic, statistic_generating */
      }
    /*
    Rprintf("\nLog of exponential function (uncentered) = % -f",sum);
    */
    moment1 = moment1 + sum;
    moment2 = moment2 + (sum * sum);
    }
  moment1 = moment1 + sum_observed; /* Add observed term */
  moment2 = moment2 + (sum_observed * sum_observed); /* Add square of observed term */
  moment1 = moment1 / (s + 1); /* First moment; note: s terms plus observed term */
  moment2 = moment2 / (s + 1); /* Second moment; note: s terms plus observed term */
  variance = moment2 - (moment1 * moment1); /* Variance */
  log_ratio_n_const = moment1 + (variance / 2.0); /* log expectation = log ratio of partition constants */
  return log_ratio_n_const;
}

double PMF_Edge_Independence_Node(int node, int d, double *input, double *theta, 
                             int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames, int *n_edges, int *heads, int *tails)
/*
input: input
output: minus energy of node i on log scale, computed under the assumption of conditional edge-independence given latent structure
*/
{
  int zero = 0;
  int one = 1;
  int i, j, edge, *number_edges, *pseudo_heads, *pseudo_tails;
  double sign, change, log_p_i_k, *statistic;
  Network nw;
  number_edges = &one;
  statistic = (double*) calloc(d,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: PMF_Independence_Node, statistic\n\n"); exit(1); }
  nw = NetworkInitialize(heads,tails,(Edge)*n_edges,(Vertex)*n,(int)*directed,(Vertex)*bipartite,zero);
  if (nw.outedges == NULL) { Rprintf("\n\ncalloc failed: PMF_Independence_Node, nw\n\n"); exit(1); }
  /* 
  Note 1: if undirected graph and i < j, undirected edge (i, j) should be stored as (i, j)
  Note 2: i, j should be integers between 1 and n
  */  
  log_p_i_k = 0.0; /* Log partition function */
  i = node + 1; /* The passed argument node of PMF_Independence_Node is in 0..n-1, whereas Partition_Function_Independence_Node assumes that it is in 1..n */
  pseudo_tails = &i; /* If undirected graph and i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
  for (j = 1; j < i; j++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
    {
    pseudo_heads = &j;
    if (EdgetreeSearch((int)*pseudo_heads,(int)*pseudo_tails,nw.outedges) == 0) sign = 1.0;
    else sign = -1.0; /* Sign is +1.0 if edge outedge (= inedge) is absent and -1.0 otherwise */
    change = Minus_Energy(d,input,theta,pseudo_heads,pseudo_tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
    log_p_i_k = log_p_i_k - ln(1.0 + e(sign * change));
    /*
    Rprintf("\nsign = %-4.0f",sign);
    Rprintf("\nchange = %-8.4f",change);
    */
    }
  pseudo_heads = &i; /* If undirected graph and i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
  for (j = i + 1; j < *n + 1; j++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
    {
    pseudo_tails = &j;
    if (EdgetreeSearch((int)*pseudo_heads,(int)*pseudo_tails,nw.outedges) == 0) sign = 1.0;
    else sign = -1.0; /* Sign is +1.0 if edge outedge (= inedge) is absent and -1.0 otherwise */
    change = Minus_Energy(d,input,theta,pseudo_heads,pseudo_tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
    log_p_i_k = log_p_i_k - ln(1.0 + e(sign * change));
    }
  free(statistic);
  NetworkDestroy(&nw);
  return log_p_i_k;
}

double PMF_Dyad_Independence_Node(int node, int d, double *input, double *theta, 
                             int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames, int *n_edges, int *heads, int *tails)
/*
input: input
output: minus energy of node i on log scale, computed under the assumption of conditional dyad-independence given latent structure
*/
{
  int zero = 0;
  int one = 1;
  int two = 2;
  int i, j, energy_0, energy_1, energy_2, energy_3, dyad, edge, *number_edges, *pseudo_heads, *pseudo_tails;
  double change, log_p_i_k, *statistic;
  Network nw;
  statistic = (double*) calloc(d,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, statistic\n\n"); exit(1); }
  nw = NetworkInitialize(heads,tails,(Edge)*n_edges,(Vertex)*n,(int)*directed,(Vertex)*bipartite,zero);
  if (nw.outedges == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, nw\n\n"); exit(1); }
  /* 
  Note 1: if undirected graph and i < j, undirected edge (i, j) should be stored as (i, j)
  Note 2: i, j should be integers between 1 and n
  */  
  /*
  Rprintf("\nnode %i",node+1);
  */
  log_p_i_k = 0.0; /* Log probability */
  i = node + 1; /* The passed argument node of PMF_Independence_Node is in 0..n-1, whereas Partition_Function_Independence_Node assumes that it is in 1..n */
  for (j = 1; j < i; j++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
    {
    /*
    0: no edge
    1: edge (i, j)
    2: edge (j, i)
    3: edges (i, j) and (j, i)
    */
    /* No edge present */
    energy_0 = 0.0; 
    /* One edge present */
    number_edges = &one;
    pseudo_heads = (int*) calloc(*number_edges,sizeof(int)); 
    if (pseudo_heads == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_heads 1\n\n"); exit(1); }
    pseudo_tails = (int*) calloc(*number_edges,sizeof(int)); 
    if (pseudo_tails == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_tails 1\n\n"); exit(1); }
    pseudo_heads[0] = i; /* Edge (i, j) */
    pseudo_tails[0] = j;
    if (EdgetreeSearch((int)*pseudo_heads,(int)*pseudo_tails,nw.outedges) != 0) dyad = 1; 
    else dyad = 0;
    energy_1 = Minus_Energy(d,input,theta,pseudo_heads,pseudo_tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
    /*
    Rprintf("\nenergy_1 = %-8.4f",change);
    */
    pseudo_heads[0] = j; /* Edge (j, i) */
    pseudo_tails[0] = i;
    if (EdgetreeSearch((int)*pseudo_heads,(int)*pseudo_tails,nw.outedges) != 0) 
      {
      if (dyad == 1) dyad = 3;
      else dyad = 2;
      }
    energy_2 = Minus_Energy(d,input,theta,pseudo_heads,pseudo_tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
    /*
    Rprintf("\nchange = %-8.4f",change);
    */
    free(pseudo_heads);
    free(pseudo_tails);
    /* Two edges present */
    number_edges = &two;
    pseudo_heads = (int*) calloc(*number_edges,sizeof(int)); 
    if (pseudo_heads == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_heads 2\n\n"); exit(1); }
    pseudo_tails = (int*) calloc(*number_edges,sizeof(int)); 
    if (pseudo_tails == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_tails 2\n\n"); exit(1); }
    pseudo_heads[0] = i; /* Edge (i, j) */
    pseudo_tails[0] = j;
    pseudo_heads[1] = j; /* Edge (j, i) */
    pseudo_tails[1] = i;
    energy_3 = Minus_Energy(d,input,theta,pseudo_heads,pseudo_tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
    /*
    Rprintf("\nchange = %-8.4f",change);
    */
    free(pseudo_heads);
    free(pseudo_tails);
    switch(dyad)
      {
      case 0: log_p_i_k = log_p_i_k + (energy_0 - ln(e(energy_0) + e(energy_1) + e(energy_2) + e(energy_3)));
      break;
      case 1: log_p_i_k = log_p_i_k + (energy_1 - ln(e(energy_0) + e(energy_1) + e(energy_2) + e(energy_3)));
      break;
      case 2: log_p_i_k = log_p_i_k + (energy_2 - ln(e(energy_0) + e(energy_1) + e(energy_2) + e(energy_3)));
      break;
      case 3: log_p_i_k = log_p_i_k + (energy_3 - ln(e(energy_0) + e(energy_1) + e(energy_2) + e(energy_3)));
      break;
      }
    /*
    if (dyad > 0) Rprintf(" %i",j);
    */
    }
  for (j = i + 1; j < *n + 1; j++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
    {
    /*
    0: no edge
    1: edge (i, j)
    2: edge (j, i)
    3: edges (i, j) and (j, i)
    */
    /* No edge present */
    energy_0 = 0.0; 
    /* One edge present */
    number_edges = &one;
    pseudo_heads = (int*) calloc(*number_edges,sizeof(int)); 
    if (pseudo_heads == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_heads 1\n\n"); exit(1); }
    pseudo_tails = (int*) calloc(*number_edges,sizeof(int)); 
    if (pseudo_tails == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_tails 1\n\n"); exit(1); }
    pseudo_heads[0] = i; /* Edge (i, j) */
    pseudo_tails[0] = j;
    if (EdgetreeSearch((int)*pseudo_heads,(int)*pseudo_tails,nw.outedges) != 0) dyad = 1; 
    else dyad = 0;
    energy_1 = Minus_Energy(d,input,theta,pseudo_heads,pseudo_tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
    /*
    Rprintf("\nenergy_1 = %-8.4f",change);
    */
    pseudo_heads[0] = j; /* Edge (j, i) */
    pseudo_tails[0] = i;
    if (EdgetreeSearch((int)*pseudo_heads,(int)*pseudo_tails,nw.outedges) != 0) 
      {
      if (dyad == 1) dyad = 3;
      else dyad = 2;
      }
    energy_2 = Minus_Energy(d,input,theta,pseudo_heads,pseudo_tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
    /*
    Rprintf("\nchange = %-8.4f",change);
    */
    free(pseudo_heads);
    free(pseudo_tails);
    /* Two edges present */
    number_edges = &two;
    pseudo_heads = (int*) calloc(*number_edges,sizeof(int)); 
    if (pseudo_heads == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_heads 2\n\n"); exit(1); }
    pseudo_tails = (int*) calloc(*number_edges,sizeof(int)); 
    if (pseudo_tails == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_tails 2\n\n"); exit(1); }
    pseudo_heads[0] = i; /* Edge (i, j) */
    pseudo_tails[0] = j;
    pseudo_heads[1] = j; /* Edge (j, i) */
    pseudo_tails[1] = i;
    energy_3 = Minus_Energy(d,input,theta,pseudo_heads,pseudo_tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
    /*
    Rprintf("\nchange = %-8.4f",change);
    */
    free(pseudo_heads);
    free(pseudo_tails);
    switch(dyad)
      {
      case 0: log_p_i_k = log_p_i_k + (energy_0 - ln(e(energy_0) + e(energy_1) + e(energy_2) + e(energy_3)));
      break;
      case 1: log_p_i_k = log_p_i_k + (energy_1 - ln(e(energy_0) + e(energy_1) + e(energy_2) + e(energy_3)));
      break;
      case 2: log_p_i_k = log_p_i_k + (energy_2 - ln(e(energy_0) + e(energy_1) + e(energy_2) + e(energy_3)));
      break;
      case 3: log_p_i_k = log_p_i_k + (energy_3 - ln(e(energy_0) + e(energy_1) + e(energy_2) + e(energy_3)));
      break;
      }
    /*
    if (dyad > 0) Rprintf(" %i",j);
    */
    }
  free(statistic);
  NetworkDestroy(&nw);
  return log_p_i_k;
}

int* Degree_Sequence(int n, int directed, int n_edges, int *heads, int *tails)
/*
input: number of nodes, indicator of directed network, number of edges, heads and tails of edge list
output: degree sequence
*/
{
  int *degree, i, j, k, sum;
  degree = (int*) calloc(n,sizeof(int)); 
  if (degree == NULL) { Rprintf("\n\ncalloc failed: Degree_Sequence, degree\n\n"); exit(1); }
  for (k = 0; k < n_edges; k++)
    {
    i = heads[k] - 1;
    j = tails[k] - 1;
    degree[i] = degree[i] + 1;
    if (directed == 0) degree[j] = degree[j] + 1;
    }
  return degree;
}

int* Degree_Freq(int n, int* degree)
/*
input: number of nodes, degree sequence
output: degree frequencies
*/
{
  int *degree_freq, i, k;
  degree_freq = (int*) calloc(n,sizeof(int)); 
  if (degree_freq == NULL) { Rprintf("\n\ncalloc failed: Degree_Freq, degree_freq\n\n"); exit(1); }
  for (i = 0; i < n; i++)
    {
    k = degree[i];
    degree_freq[k] = degree_freq[k] + 1;
    }
  return degree_freq;
}

double* Block_Degree_Freq(int n, int *degree, int block_number, int *block_size, int *block_indicator)
/*
input: number of nodes, degree sequence, number of blocks, size of blocks, indicator of block membership
output: relative frequencies of degree by block
*/
{
  int i, k;
  double *block_degree_freq;
  block_degree_freq = (double*) calloc(block_number,sizeof(double)); 
  if (block_degree_freq == NULL) { Rprintf("\n\ncalloc failed: Block_Degree_Freq, block_degree_freq\n\n"); exit(1); }
  for (i = 0; i < n; i++)
    {
    k = block_indicator[i];
    block_degree_freq[k] = block_degree_freq[k] + degree[i];
    }
  for (k = 0; k < block_number; k++)
    {
    if (block_size[k] == 0) block_degree_freq[k] = 0.0;
    else block_degree_freq[k] = block_degree_freq[k] / block_size[k];
    }
  return block_degree_freq;
}

