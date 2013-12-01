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

#include "h_ergm_mcmc.h"

void P_Edge_Independence(int *number_terms, int *number_parameters, double *input, double *theta,  int *n, int *directed, int *bipartite, char **funnames, char **sonames, double *p)
/*
input: undirected graph; number of terms; number of parameters;  input vector; parameter vector; number of nodes; other variables
output: probabilities of edges between nodes i and j on log scale, computed under the assumption of conditional edge-independence given latent structure,
and ordered in accordance with i < j
*/
{
  int one = 1;
  int index, i, j, *number_edges, *heads, *tails;
  double log_odds, *statistic;
  number_edges = &one;
  statistic = (double*) calloc(*number_parameters,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: P_Edge_Independence, statistic\n\n"); error("Error: out of memory"); }
  /* 
  Note 1: if undirected graph and i < j, undirected edge (i, j) should be stored as (i, j)
  Note 2: i, j should be integers between 1 and n
  */  
  index = 0;
  for (i = 1; i < *n; i++) /* Row i */
    {
    heads = &i; 
    for (j = i + 1; j < *n + 1; j++) /* Row i, column j > i (undirected, directed graph) */
      {
      tails = &j;
      log_odds = Minus_Energy(*number_parameters,input,theta,heads,tails,number_edges,n,directed,bipartite,number_terms,funnames,sonames,statistic); /* Compute log-odds of probability of edge statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
      p[index] = -ln(1.0 + e(-log_odds));
      index = index + 1;
      }
    }
  free(statistic);
}

double Partition_Function_Edge_Independence(latentstructure *ls, ergmstructure *ergm, double *input, double *theta, 
                                         int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale, computed under the assumption of conditional edge-independence given latent structure
*/
{
  int one = 1;
  int i, j, *number_edges, *heads, *tails;
  double a, b, *statistic;
  number_edges = &one;
  statistic = (double*) calloc(ergm->d,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: Partition_Function_Edge_Independence, statistic\n\n"); error("Error: out of memory"); }
  /* 
  Note 1: if undirected graph and i < j, undirected edge (i, j) should be stored as (i, j)
  Note 2: i, j should be integers between 1 and n
  */  
  a = 0.0; /* Log partition function */
  for (i = 1; i < ls->n + 1; i++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
    {
    heads = &i; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
      {
      for (j = i + 1; j < ls->n + 1; j++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
        {
        tails = &j;
        b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
        a = a + ln(1.0 + e(b));
        }
      }
    }
  free(statistic);
  return a;
}

double Partition_Function_Dyad_Independence(latentstructure *ls, ergmstructure *ergm, double *input, double *theta, 
                                         int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale, computed under the assumption of conditional dyad-independence given latent structure
*/
{
  int one = 1;
  int two = 2;
  int i, j, *number_edges, *heads, *tails;
  double a, b, *statistic, sum;
  statistic = (double*) calloc(ergm->d,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: Partition_Function_Dyad_Independence, statistic\n\n"); error("Error: out of memory"); }
  /* 
  Note 1: if undirected graph and i < j, undirected edge (i, j) should be stored as (i, j)
  Note 2: i, j should be integers between 1 and n
  */  
  a = 0.0; /* Log partition function */
  for (i = 1; i < ls->n + 1; i++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
    {
    for (j = i + 1; j < ls->n + 1; j++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
      {
      /* No edge present */
      sum = 1.0; 
      /* One edge present */
      number_edges = &one;
      heads = (int*) calloc(*number_edges,sizeof(int)); 
      if (heads == NULL) { Rprintf("\n\ncalloc failed: Partition_Function_Dyad_Independence, heads 1\n\n"); error("Error: out of memory"); }
      tails = (int*) calloc(*number_edges,sizeof(int)); 
      if (tails == NULL) { Rprintf("\n\ncalloc failed: Partition_Function_Dyad_Independence, tails 1\n\n"); error("Error: out of memory"); }
      heads[0] = i; /* Edge (i, j) */
      tails[0] = j;
      b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
      sum = sum + e(b);
      heads[0] = j; /* Edge (j, i) */
      tails[0] = i;
      b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
      sum = sum + e(b);
      free(heads);
      free(tails);
      /* Two edges present */
      number_edges = &two;
      heads = (int*) calloc(*number_edges,sizeof(int)); 
      if (heads == NULL) { Rprintf("\n\ncalloc failed: Partition_Function_Dyad_Independence, heads 2\n\n"); error("Error: out of memory"); }
      tails = (int*) calloc(*number_edges,sizeof(int)); 
      if (tails == NULL) { Rprintf("\n\ncalloc failed: Partition_Function_Dyad_Independence, tails 2\n\n"); error("Error: out of memory"); }
      heads[0] = i; /* Edge (i, j) */
      tails[0] = j;
      heads[1] = j; /* Edge (j, i) */
      tails[1] = i;
      b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
      sum = sum + e(b);
      free(heads);
      free(tails);
      /* Take log */
      a = a + ln(sum);
      }
    }
  free(statistic);
  return a;
}

double PMF_Independence(latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input, double *theta, 
                        int *n_edges, int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: probability mass on log scale, computed under the assumption of dyad-dependence
*/
{
  double a, log_p, *statistic, u;
  statistic = (double*) calloc(ergm->d,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: PMF_Independence, statistic\n\n"); error("Error: out of memory"); }
  u = Minus_Energy(ergm->d,input,theta,heads,tails,n_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> on log scale */
  /*
  Rprintf("\nPMF_Independence: minus potential energy function = %f",- u);
  */
  if (*directed == 0) a = Partition_Function_Edge_Independence(ls,ergm,input,theta,n,directed,bipartite,nterms,funnames,sonames); /* Log partition function: undirected case */
  else a = Partition_Function_Dyad_Independence(ls,ergm,input,theta,n,directed,bipartite,nterms,funnames,sonames); /* Log partition function: directed case */
  /*
  Rprintf("\nPMF_Independence: log partition function = %f",a);
  */
  log_p = u - a; /* Probability mass */
  free(statistic);
  return log_p;
}

double PMF_i_k_Node(int i, int l, latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input_proposal, 
                    int *n_edges, int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: node i, catogory l, latent structure, ergm structure
output: conditional PMF of graph given latent structure 
*/
{
  int k;
  double log_p_i_k, *theta;
  k = ls->indicator[i]; /* Store indicator */
  ls->indicator[i] = l; /* Set indicator */
  /*
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_proposal); /* Set input given indicator
  */
  Set_Input_Indicator(ergm->terms,ergm->hierarchical,ls->number,ls->n,i,l,input_proposal); /* Set input given indicator; reset in Gibbs_Indicators_Independence */
  theta = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter */
  if (*directed == 0) log_p_i_k = PMF_Edge_Independence_Node(i,ergm->d,input_proposal,theta,n,directed,bipartite,nterms,funnames,sonames,n_edges,heads,tails); /* Probability mass under given indicator */
  else log_p_i_k = PMF_Dyad_Independence_Node(i,ergm->d,input_proposal,theta,n,directed,bipartite,nterms,funnames,sonames,n_edges,heads,tails); /* Probability mass under given indicator */
  ls->indicator[i] = k; /* Reset indicator */
  /*
  Rprintf("\nPMF_i_k: %f",log_p_i_k);
  */
  free(theta);
  return log_p_i_k;
}

double PMF_Edge_Independence_Node(int node, int d, double *input, double *theta, 
                             int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames, int *n_edges, int *heads, int *tails)
/*
input: input
output: minus energy of node i on log scale, computed under the assumption of conditional edge-independence given latent structure
*/
{
  int one = 1;
  int i, j, edge, *lasttoggle, *number_edges, *pseudo_heads, *pseudo_tails;
  double sign, change, log_p_i_k, *statistic;
  Network nw;
  number_edges = &one;
  statistic = (double*) calloc(d,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: PMF_Independence_Node, statistic\n\n"); error("Error: out of memory"); }
  nw = NetworkInitialize(tails,heads,(Edge)*n_edges,(Vertex)*n,(int)*directed,(Vertex)*bipartite, 0, 0, NULL);
  if (nw.outedges == NULL) { Rprintf("\n\ncalloc failed: PMF_Independence_Node, nw\n\n"); error("Error: out of memory"); }
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
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, statistic\n\n"); error("Error: out of memory"); }
  nw = NetworkInitialize(tails,heads,(Edge)*n_edges,(Vertex)*n,(int)*directed,(Vertex)*bipartite, 0, 0, NULL);
  if (nw.outedges == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, nw\n\n"); error("Error: out of memory"); }
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
    if (pseudo_heads == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_heads 1\n\n"); error("Error: out of memory"); }
    pseudo_tails = (int*) calloc(*number_edges,sizeof(int)); 
    if (pseudo_tails == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_tails 1\n\n"); error("Error: out of memory"); }
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
    if (pseudo_heads == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_heads 2\n\n"); error("Error: out of memory"); }
    pseudo_tails = (int*) calloc(*number_edges,sizeof(int)); 
    if (pseudo_tails == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_tails 2\n\n"); error("Error: out of memory"); }
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
    if (pseudo_heads == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_heads 1\n\n"); error("Error: out of memory"); }
    pseudo_tails = (int*) calloc(*number_edges,sizeof(int)); 
    if (pseudo_tails == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_tails 1\n\n"); error("Error: out of memory"); }
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
    if (pseudo_heads == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_heads 2\n\n"); error("Error: out of memory"); }
    pseudo_tails = (int*) calloc(*number_edges,sizeof(int)); 
    if (pseudo_tails == NULL) { Rprintf("\n\ncalloc failed: PMF_Dyad_Independence_Node, pseudo_tails 2\n\n"); error("Error: out of memory"); }
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

void Gibbs_Indicators_Independence(latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input_proposal, 
                       int *n_edges, int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames, double *q_i)
/*
input: latent structure, ergm structure
output: indicators
*/
{
  int i, k, *sample, sample_size;
  double center, log_p_i_k, p_i_k, *p_i, sum, u;
  p_i = (double*) calloc(ls->number,sizeof(double));
  if (p_i == NULL) { Rprintf("\n\ncalloc failed: Gibbs_Indicators_Independence, p_i\n\n"); error("Error: out of memory"); }
  sample = (int*) calloc(ls->n,sizeof(int));
  if (sample == NULL) { Rprintf("\n\ncalloc failed: Gibbs_Indicators_Independence, sample\n\n"); error("Error: out of memory"); }
  for (k = 0; k < ls->number; k++) /* Reset size */
    {
    ls->size[k] = 0;
    }
  sample_size = trunc(ls->n / 10.0); 
  if (sample_size < 10) sample_size = 10;
  for (k = 0; k < sample_size; k++)
    {
    i = Sample_Discrete(q_i);
    sample[i] = 1;
    }
  for (i = 0; i < ls->n; i++) /* Node i */
    {
    if (sample[i] == 1) /* Indicator of node i updated: y/n */ /* 333 */
      {
      /*
      Rprintf("\nnode %-3i",i);
      */
      sum = 0.0;
      for (k = 0; k < ls->number; k++) /* Category k */
        {
        log_p_i_k = PMF_i_k_Node(i,k,ls,ergm,heads,tails,input_proposal,n_edges,n,directed,bipartite,nterms,funnames,sonames);
        if (k == 0)
          {
          center = log_p_i_k;
          log_p_i_k = 0;
          }
        else log_p_i_k = log_p_i_k - center;
        p_i_k = e(log_p_i_k);
        if ((ls->p[k] * p_i_k) < epsilon) p_i[k] = epsilon; /* Mass */
        else p_i[k] = ls->p[k] * p_i_k;
        sum = sum + p_i[k];
        }   
      for (k = 0; k < ls->number; k++) /* Probability mass */
        {
        p_i[k] = p_i[k] / sum; /* Underflow impossible since sum is at least ls->number * epsilon: see above */
        /*
        Rprintf(" %f (%f)",p_i[k], ls->theta[0][k]);
        */
        }
      k = Sample_Discrete(p_i); /* Sample full conditional of category indicators */
      ls->indicator[i] = k; /* Update indicator */ 
      Set_Input_Indicator(ergm->terms,ergm->hierarchical,ls->number,ls->n,i,k,input_proposal); /* Since the full conditionals of the category indicators of nodes i+1..n may depend on the category indicator of node i, input_proposal must be updated */
      }
    else k = ls->indicator[i];
    ls->size[k] = ls->size[k] + 1; /* Update size */
    }   
  free(p_i);
  free(sample);
}

int** Edge_List_Blocks(latentstructure *ls, int *block, int *total_number_edges, int *total_heads, int *total_tails)
/*
input: latent structure, block, number of edges and edge list in terms of heads and tails, number of labels of included blocks
output: number of edges and edge_list of members of included blocks
*/
{
  int number_nodes, **edge_list, head, i, k, tail, included, number_blocks, number_edges, *label;
  label = (int*) calloc(ls->n,sizeof(int));
  if (label == NULL) { Rprintf("\n\ncalloc failed: Edge_List_Blocks, label\n\n"); error("Error: out of memory"); }
  number_blocks = block[0]; /* Number of blocks included; the labels of included blocks are stored in block[1], ..., block[number_blocks] */
  number_nodes = 0; /* Number of nodes which are members of included blocks */
  for (i = 0; i < ls->n; i++)
    {
    included = 0; /* Indicator of whether node is member of included blocks */
    k = 0; 
    while ((k < number_blocks) && (included == 0)) /* Check whether node is member of included blocks */
      { 
      k = k + 1;
      if (ls->indicator[i] == block[k]) included = 1; /* Conclusion: node is member of included blocks */ 
      }
    if (included == 1) /* If node is member of included blocks, then... */
      {
      number_nodes = number_nodes + 1; 
      label[i] = number_nodes; /* ...change its label so that members of included blocks have labels 1..number_nodes */
      }
    }
  edge_list = (int**) calloc(3,sizeof(int*)); /* Edge list of included blocks */
  edge_list[0] = (int*) calloc(1,sizeof(int)); /* Stores number of included blocks */
  edge_list[1] = NULL; /* Stores heads of edge_list of included blocks */
  edge_list[2] = NULL; /* Stores tails of edge_list of included blocks */
  number_edges = 0;
  for (i = 0; i < *total_number_edges; i++) /* Go through edge_list and extract edges between members of included blocks */
    {
    head = total_heads[i];
    tail = total_tails[i];
    if ((label[head-1] > 0) && (label[tail-1] > 0)) /* Labels of members of included blocks are positive, others are by default 0 */
      {
      number_edges = number_edges + 1; /* Number of edges between members of included blocks */
      edge_list[1] = realloc(edge_list[1],number_edges*sizeof(int));
      edge_list[2] = realloc(edge_list[2],number_edges*sizeof(int));
      edge_list[1][number_edges-1] = label[head-1]; /* Heads of edges between members of included blocks */
      edge_list[2][number_edges-1] = label[tail-1]; /* Tails of edges between members of included blocks */
      }
    }
  edge_list[0][0] = number_edges; /* Number of edges between members of included blocks */
  free(label);
  return edge_list;
}

int* Degree_Sequence(int n, int directed, int n_edges, int *heads, int *tails)
/*
input: number of nodes, indicator of directed network, number of edges, heads and tails of edge list
output: degree sequence
*/
{
  int *degree, i, j, k, sum;
  degree = (int*) calloc(n,sizeof(int)); 
  if (degree == NULL) { Rprintf("\n\ncalloc failed: Degree_Sequence, degree\n\n"); error("Error: out of memory"); }
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
  if (degree_freq == NULL) { Rprintf("\n\ncalloc failed: Degree_Freq, degree_freq\n\n"); error("Error: out of memory"); }
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
  if (block_degree_freq == NULL) { Rprintf("\n\ncalloc failed: Block_Degree_Freq, block_degree_freq\n\n"); error("Error: out of memory"); }
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

double Within_Block_Partition_Function_2_Graph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale given n = 2 nodes and undirected graphs
*/
{
  int count, edges, *heads, i, j, n_k, *number_edges, *tails, **y;
  double a, b, *statistic, sum;
  n_k = 2;
  y = (int**) calloc(n_k,sizeof(int*));
  for (i = 0; i < n_k; i++)
    {
    y[i] = (int*) calloc(n_k,sizeof(int));
    }
  statistic = (double*) calloc(ergm->d,sizeof(double));
  sum = 0.0;
  for (y[0][1] = 0; y[0][1] < 2; y[0][1]++)
    { 
    edges = y[0][1];
    number_edges = &edges;
    heads = (int*) calloc(*number_edges,sizeof(int));
    tails = (int*) calloc(*number_edges,sizeof(int));
    count = -1;
    for (i = 0; i < n_k - 1; i++)
      {
      for (j = i + 1; j < n_k; j++)
        {
        if (y[i][j] == 1) 
          {
          count = count + 1;
          heads[count] = block_members[i] + 1; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
          tails[count] = block_members[j] + 1;
          }
        }
      }
    b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); 
    sum = sum + e(b);
    free(heads);
    free(tails);
    }
  a = ln(sum);
  for (i = 0; i < n_k; i++)
    {
    free(y[i]);
    }
  free(y);
  free(statistic);
  return a;
}

double Within_Block_Partition_Function_3_Graph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale given n = 3 nodes and undirected graphs
*/
{
  int count, edges, *heads, i, j, n_k, *number_edges, *tails, **y;
  double a, b, *statistic, sum;
  n_k = 3;
  y = (int**) calloc(n_k,sizeof(int*));
  for (i = 0; i < n_k; i++)
    {
    y[i] = (int*) calloc(n_k,sizeof(int));
    }
  statistic = (double*) calloc(ergm->d,sizeof(double));
  sum = 0.0;
  for (y[0][1] = 0; y[0][1] < 2; y[0][1]++)
    {
    for (y[0][2] = 0; y[0][2] < 2; y[0][2]++)
      {
      for (y[1][2] = 0; y[1][2] < 2; y[1][2]++)
        {
        edges = y[0][1] + y[0][2] + y[1][2];
        number_edges = &edges;
        heads = (int*) calloc(*number_edges,sizeof(int));
        tails = (int*) calloc(*number_edges,sizeof(int));
        count = -1;
        for (i = 0; i < n_k - 1; i++)
          {
          for (j = i + 1; j < n_k; j++)
            {
            if (y[i][j] == 1) 
              {
              count = count + 1;
              heads[count] = block_members[i] + 1; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
              tails[count] = block_members[j] + 1;
              }
            }
          }
        b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); 
        sum = sum + e(b);
        free(heads);
        free(tails);
        }
      } 
    }
  a = ln(sum);
  for (i = 0; i < n_k; i++)
    {
    free(y[i]);
    }
  free(y);
  free(statistic);
  return a;
}

double Within_Block_Partition_Function_4_Graph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale given n = 4 nodes and undirected graphs
*/
{
  int count, edges, *heads, i, j, n_k, *number_edges, *tails, **y;
  double a, b, *statistic, sum;
  n_k = 4;
  y = (int**) calloc(n_k,sizeof(int*));
  for (i = 0; i < n_k; i++)
    {
    y[i] = (int*) calloc(n_k,sizeof(int));
    }
  statistic = (double*) calloc(ergm->d,sizeof(double));
  sum = 0.0;
  for (y[0][1] = 0; y[0][1] < 2; y[0][1]++)
    {
    for (y[0][2] = 0; y[0][2] < 2; y[0][2]++)
      {
      for (y[0][3] = 0; y[0][3] < 2; y[0][3]++)
        {
        for (y[1][2] = 0; y[1][2] < 2; y[1][2]++)
          {
          for (y[1][3] = 0; y[1][3] < 2; y[1][3]++)
            {
            for (y[2][3] = 0; y[2][3] < 2; y[2][3]++)
              {
              edges = y[0][1] + y[0][2] + y[0][3] + y[1][2] + y[1][3] + y[2][3];
              number_edges = &edges;
              heads = (int*) calloc(*number_edges,sizeof(int));
              tails = (int*) calloc(*number_edges,sizeof(int));
              count = -1;
              for (i = 0; i < n_k - 1; i++)
                {
                for (j = i + 1; j < n_k; j++)
                  {
                  if (y[i][j] == 1) 
                    {
                    count = count + 1;
                    heads[count] = block_members[i] + 1; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
                    tails[count] = block_members[j] + 1;
                    }
                  }
                }
              b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic);
              sum = sum + e(b);
              free(heads);
              free(tails);
              }
            }    
          }  
        }
      } 
    }
  a = ln(sum);
  for (i = 0; i < n_k; i++)
    {
    free(y[i]);
    }
  free(y);
  free(statistic);
  return a;
}

double Within_Block_Partition_Function_5_Graph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale given n = 5 nodes and undirected graphs
*/
{
  int count, edges, *heads, i, j, n_k, *number_edges, *tails, **y;
  double a, b, *statistic, sum;
  n_k = 5;
  y = (int**) calloc(n_k,sizeof(int*));
  for (i = 0; i < n_k; i++)
    {
    y[i] = (int*) calloc(n_k,sizeof(int));
    }
  statistic = (double*) calloc(ergm->d,sizeof(double));
  sum = 0.0;
  for (y[0][1] = 0; y[0][1] < 2; y[0][1]++)
    {
    for (y[0][2] = 0; y[0][2] < 2; y[0][2]++)
      {
      for (y[0][3] = 0; y[0][3] < 2; y[0][3]++)
        {
        for (y[0][4] = 0; y[0][4] < 2; y[0][4]++)
          {
          for (y[1][2] = 0; y[1][2] < 2; y[1][2]++)
            {
            for (y[1][3] = 0; y[1][3] < 2; y[1][3]++)
              {
              for (y[1][4] = 0; y[1][4] < 2; y[1][4]++)
                {
                for (y[2][3] = 0; y[2][3] < 2; y[2][3]++)
                  { 
                  for (y[2][4] = 0; y[2][4] < 2; y[2][4]++)
                    { 
                    for (y[3][4] = 0; y[3][4] < 2; y[3][4]++)
                      { 
                      edges = y[0][1] + y[0][2] + y[0][3] + y[0][4] + y[1][2] + y[1][3] + y[1][4] + y[2][3] + y[2][4] + y[3][4];
                      number_edges = &edges;
                      heads = (int*) calloc(*number_edges,sizeof(int));
                      tails = (int*) calloc(*number_edges,sizeof(int));
                      count = -1;
                      for (i = 0; i < n_k - 1; i++)
                        {
                        for (j = i + 1; j < n_k; j++)
                          {
                          if (y[i][j] == 1)         
                            {
                            count = count + 1;
                            heads[count] = block_members[i] + 1; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
                            tails[count] = block_members[j] + 1;
                            }
                          }
                        }
                      b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic);
                      sum = sum + e(b);
                      free(heads);
                      free(tails);
                      }
                    }    
                  }  
                }
              } 
            }
          }
        }
      }
    }  
  a = ln(sum);
  for (i = 0; i < n_k; i++)
    {
    free(y[i]);
    }
  free(y);
  free(statistic);
  return a;
}

double Within_Block_Partition_Function_6_Graph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale given n = 6 nodes and undirected graphs
*/
{
  int count, edges, *heads, i, j, n_k, *number_edges, *tails, **y;
  double a, b, *statistic, sum;
  n_k = 6;
  y = (int**) calloc(n_k,sizeof(int*));
  for (i = 0; i < n_k; i++)
    {
    y[i] = (int*) calloc(n_k,sizeof(int));
    }
  statistic = (double*) calloc(ergm->d,sizeof(double));
  sum = 0.0;
  for (y[0][1] = 0; y[0][1] < 2; y[0][1]++)
    {
    for (y[0][2] = 0; y[0][2] < 2; y[0][2]++)
      {
      for (y[0][3] = 0; y[0][3] < 2; y[0][3]++)
        {
        for (y[0][4] = 0; y[0][4] < 2; y[0][4]++)
          {
          for (y[0][5] = 0; y[0][5] < 2; y[0][5]++)
            {
            for (y[1][2] = 0; y[1][2] < 2; y[1][2]++)
              {
              for (y[1][3] = 0; y[1][3] < 2; y[1][3]++)
                {
                for (y[1][4] = 0; y[1][4] < 2; y[1][4]++)
                  {
                  for (y[1][5] = 0; y[1][5] < 2; y[1][5]++)
                    {
                    for (y[2][3] = 0; y[2][3] < 2; y[2][3]++)
                      { 
                      for (y[2][4] = 0; y[2][4] < 2; y[2][4]++)
                        { 
                        for (y[2][5] = 0; y[2][5] < 2; y[2][5]++)
                          { 
                          for (y[3][4] = 0; y[3][4] < 2; y[3][4]++)
                            { 
                            for (y[3][5] = 0; y[3][5] < 2; y[3][5]++)
                              { 
                              for (y[4][5] = 0; y[4][5] < 2; y[4][5]++)
                                { 
                                edges = y[0][1] + y[0][2] + y[0][3] + y[0][4] + y[0][5] + y[1][2] + y[1][3] + y[1][4] + y[1][5] + y[2][3] + y[2][4] + y[2][5] + y[3][4] + y[3][5] + y[4][5];
                                number_edges = &edges;
                                heads = (int*) calloc(*number_edges,sizeof(int));
                                tails = (int*) calloc(*number_edges,sizeof(int));
                                count = -1;
                                for (i = 0; i < n_k - 1; i++)
                                  {
                                  for (j = i + 1; j < n_k; j++)
                                    {
                                    if (y[i][j] == 1)         
                                      {
                                      count = count + 1;
                                      heads[count] = block_members[i] + 1; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
                                      tails[count] = block_members[j] + 1;
                                      }
                                    }
                                  }
                                b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic);
                                sum = sum + e(b);
                                free(heads);
                                free(tails);
                                }
                              }    
                            }
                          }
                        }
                      }
                    }        
                  }  
                }
              } 
            }
          }
        }
      }
    }  
  a = ln(sum);
  for (i = 0; i < n_k; i++)
    {
    free(y[i]);
    }
  free(y);
  free(statistic);
  return a;
}

double Within_Block_Partition_Function_7_Graph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale given n = 7 nodes and undirected graphs
*/
{
  int count, edges, *heads, i, j, n_k, *number_edges, *tails, **y;
  double a, b, *statistic, sum;
  n_k = 7;
  y = (int**) calloc(n_k,sizeof(int*));
  for (i = 0; i < n_k; i++)
    {
    y[i] = (int*) calloc(n_k,sizeof(int));
    }
  statistic = (double*) calloc(ergm->d,sizeof(double));
  sum = 0.0;
  for (y[0][1] = 0; y[0][1] < 2; y[0][1]++)
    {
    for (y[0][2] = 0; y[0][2] < 2; y[0][2]++)
      {
      for (y[0][3] = 0; y[0][3] < 2; y[0][3]++)
        {
        for (y[0][4] = 0; y[0][4] < 2; y[0][4]++)
          {
          for (y[0][5] = 0; y[0][5] < 2; y[0][5]++)
            {
            for (y[0][6] = 0; y[0][6] < 2; y[0][6]++)
              {
              for (y[1][2] = 0; y[1][2] < 2; y[1][2]++)
                {
                for (y[1][3] = 0; y[1][3] < 2; y[1][3]++)
                  {
                  for (y[1][4] = 0; y[1][4] < 2; y[1][4]++)
                    {
                    for (y[1][5] = 0; y[1][5] < 2; y[1][5]++)
                      {
                      for (y[1][6] = 0; y[1][6] < 2; y[1][6]++)
                        {
                        for (y[2][3] = 0; y[2][3] < 2; y[2][3]++)
                          { 
                          for (y[2][4] = 0; y[2][4] < 2; y[2][4]++)
                            { 
                            for (y[2][5] = 0; y[2][5] < 2; y[2][5]++)
                              { 
                              for (y[2][6] = 0; y[2][6] < 2; y[2][6]++)
                                { 
                                for (y[3][4] = 0; y[3][4] < 2; y[3][4]++)
                                  { 
                                  for (y[3][5] = 0; y[3][5] < 2; y[3][5]++)
                                    { 
                                    for (y[3][6] = 0; y[3][6] < 2; y[3][6]++)
                                      { 
                                      for (y[4][5] = 0; y[4][5] < 2; y[4][5]++)
                                        { 
                                        for (y[4][6] = 0; y[4][6] < 2; y[4][6]++)
                                          { 
                                          for (y[5][6] = 0; y[5][6] < 2; y[5][6]++)
                                            { 
                                            edges = y[0][1] + y[0][2] + y[0][3] + y[0][4] + y[0][5] + y[0][6] + y[1][2] + y[1][3] + y[1][4] + y[1][5] + y[1][6] + y[2][3] + y[2][4] + y[2][5] + y[2][6] + y[3][4] + y[3][5] + y[3][6] + y[4][5] + y[4][6] + y[5][6];
                                            number_edges = &edges;
                                            heads = (int*) calloc(*number_edges,sizeof(int));
                                            tails = (int*) calloc(*number_edges,sizeof(int));
                                            count = -1;
                                            for (i = 0; i < n_k - 1; i++)
                                              {
                                              for (j = i + 1; j < n_k; j++)
                                                {
                                                if (y[i][j] == 1)         
                                                  {
                                                  count = count + 1;
                                                  heads[count] = block_members[i] + 1; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
                                                  tails[count] = block_members[j] + 1;
                                                  }
                                                }
                                              }
                                            b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic);
                                            sum = sum + e(b);
                                            free(heads);
                                            free(tails);
                                            }
                                          }  
                                        }
                                      }  
                                    }
                                  }
                                }
                              }    
                            }
                          }
                        }
                      }
                    }        
                  }  
                }
              } 
            }
          }
        }
      }
    }  
  a = ln(sum);
  for (i = 0; i < n_k; i++)
    {
    free(y[i]);
    }
  free(y);
  free(statistic);
  return a;
}

double Within_Block_Partition_Function_2_Digraph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale given n = 2 nodes and directed graphs
*/
{
  int count, edges, *heads, i, j, n_k, *number_edges, *tails, **y;
  double a, b, *statistic, sum;
  n_k = 2;
  y = (int**) calloc(n_k,sizeof(int*));
  for (i = 0; i < n_k; i++)
    {
    y[i] = (int*) calloc(n_k,sizeof(int));
    }
  statistic = (double*) calloc(ergm->d,sizeof(double));
  sum = 0.0;
  for (y[0][1] = 0; y[0][1] < 2; y[0][1]++)
    {
    for (y[1][0] = 0; y[1][0] < 2; y[1][0]++)
      {
      edges = y[0][1] + y[1][0];
      number_edges = &edges;
      heads = (int*) calloc(*number_edges,sizeof(int));
      tails = (int*) calloc(*number_edges,sizeof(int));
      count = -1;
      for (i = 0; i < n_k; i++)
        {
        for (j = 0; j < n_k; j++)
          {
          if (y[i][j] == 1) /* If directed, y[i][i] == 0 */ 
            {
            count = count + 1;
            heads[count] = block_members[i] + 1; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
            tails[count] = block_members[j] + 1;
            }
          }
        }
      b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); 
      sum = sum + e(b);
      free(heads);
      free(tails);
      } 
    }
  a = ln(sum);
  for (i = 0; i < n_k; i++)
    {
    free(y[i]);
    }
  free(y);
  free(statistic);
  return a;
}

double Within_Block_Partition_Function_3_Digraph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale given n = 3 nodes and directed graphs
*/
{
  int count, edges, *heads, i, j, n_k, *number_edges, *tails, **y;
  double a, b, *statistic, sum;
  n_k = 3;
  y = (int**) calloc(n_k,sizeof(int*));
  for (i = 0; i < n_k; i++)
    {
    y[i] = (int*) calloc(n_k,sizeof(int));
    }
  statistic = (double*) calloc(ergm->d,sizeof(double));
  sum = 0.0;
  for (y[0][1] = 0; y[0][1] < 2; y[0][1]++)
    {
    for (y[1][0] = 0; y[1][0] < 2; y[1][0]++)
      {
      for (y[0][2] = 0; y[0][2] < 2; y[0][2]++)
        {
        for (y[2][0] = 0; y[2][0] < 2; y[2][0]++)
          {
          for (y[1][2] = 0; y[1][2] < 2; y[1][2]++)
            {
            for (y[2][1] = 0; y[2][1] < 2; y[2][1]++)
              {
              edges = y[0][1] + y[1][0] + y[0][2] + y[2][0] + y[1][2] + y[2][1];
              number_edges = &edges;
              heads = (int*) calloc(*number_edges,sizeof(int));
              tails = (int*) calloc(*number_edges,sizeof(int));
              count = -1;
              for (i = 0; i < n_k; i++)
                {
                for (j = 0; j < n_k; j++)
                  {
                  if (y[i][j] == 1) /* If directed, y[i][i] == 0 */ 
                    {
                    count = count + 1;
                    heads[count] = block_members[i] + 1; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
                    tails[count] = block_members[j] + 1;
                    }
                  }
                }
              b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); 
              sum = sum + e(b);
              free(heads);
              free(tails);
              }
            }
          }    
        }
      } 
    }
  a = ln(sum);
  for (i = 0; i < n_k; i++)
    {
    free(y[i]);
    }
  free(y);
  free(statistic);
  return a;
}

double Within_Block_Partition_Function_4_Digraph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale given n = 4 nodes and directed graphs
*/
{
  int count, edges, *heads, i, j, n_k, *number_edges, *tails, **y;
  double a, b, *statistic, sum;
  n_k = 4;
  y = (int**) calloc(n_k,sizeof(int*));
  for (i = 0; i < n_k; i++)
    {
    y[i] = (int*) calloc(n_k,sizeof(int));
    }
  statistic = (double*) calloc(ergm->d,sizeof(double));
  sum = 0.0;
  for (y[0][1] = 0; y[0][1] < 2; y[0][1]++)
    {
    for (y[1][0] = 0; y[1][0] < 2; y[1][0]++)
      {
      for (y[0][2] = 0; y[0][2] < 2; y[0][2]++)
        {
        for (y[2][0] = 0; y[2][0] < 2; y[2][0]++)
          {
          for (y[0][3] = 0; y[0][3] < 2; y[0][3]++)
            {
            for (y[3][0] = 0; y[3][0] < 2; y[3][0]++)
              {
              for (y[1][2] = 0; y[1][2] < 2; y[1][2]++)
                {
                for (y[2][1] = 0; y[2][1] < 2; y[2][1]++)
                  {
                  for (y[1][3] = 0; y[1][3] < 2; y[1][3]++)
                    {
                    for (y[3][1] = 0; y[3][1] < 2; y[3][1]++)
                      {
                      for (y[2][3] = 0; y[2][3] < 2; y[2][3]++)
                        {
                        for (y[3][2] = 0; y[3][2] < 2; y[3][2]++)
                          {
                          edges = y[0][1] + y[1][0] + y[0][2] + y[2][0] + y[0][3] + y[3][0] + y[1][2] + y[2][1] + y[1][3] + y[3][1] + y[2][3] + y[3][2];
                          number_edges = &edges;
                          heads = (int*) calloc(*number_edges,sizeof(int));
                          tails = (int*) calloc(*number_edges,sizeof(int));
                          count = -1;
                          for (i = 0; i < n_k; i++)
                            {
                            for (j = 0; j < n_k; j++)
                              {
                              if (y[i][j] == 1) /* If directed, y[i][i] == 0 */ 
                                {
                                count = count + 1;
                                heads[count] = block_members[i] + 1; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
                                tails[count] = block_members[j] + 1;
                                }
                              }
                            }
                          b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); 
                          sum = sum + e(b);
                          free(heads);
                          free(tails);
                          }
                        }   
                      }
                    }
                  }
                }    
              }  
            }
          }    
        }
      } 
    }
  a = ln(sum);
  for (i = 0; i < n_k; i++)
    {
    free(y[i]);
    }
  free(y);
  free(statistic);
  return a;
}

double Within_Block_Partition_Function_5_Digraph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale given n = 5 nodes and directed graphs
*/
{
  int count, edges, *heads, i, j, n_k, *number_edges, *tails, **y;
  double a, b, *statistic, sum;
  n_k = 5;
  y = (int**) calloc(n_k,sizeof(int*));
  for (i = 0; i < n_k; i++)
    {
    y[i] = (int*) calloc(n_k,sizeof(int));
    }
  statistic = (double*) calloc(ergm->d,sizeof(double));
  sum = 0.0;
  for (y[0][1] = 0; y[0][1] < 2; y[0][1]++)
    {
    for (y[1][0] = 0; y[1][0] < 2; y[1][0]++)
      {
      for (y[0][2] = 0; y[0][2] < 2; y[0][2]++)
        {
        for (y[2][0] = 0; y[2][0] < 2; y[2][0]++)
          {
          for (y[0][3] = 0; y[0][3] < 2; y[0][3]++)
            {
            for (y[3][0] = 0; y[3][0] < 2; y[3][0]++)
              {
              for (y[0][4] = 0; y[0][4] < 2; y[0][4]++)
                {
                for (y[4][0] = 0; y[4][0] < 2; y[4][0]++)
                  {
                  for (y[1][2] = 0; y[1][2] < 2; y[1][2]++)
                    {
                    for (y[2][1] = 0; y[2][1] < 2; y[2][1]++)
                      {
                      for (y[1][3] = 0; y[1][3] < 2; y[1][3]++)
                        {
                        for (y[3][1] = 0; y[3][1] < 2; y[3][1]++)
                          {
                          for (y[1][4] = 0; y[1][4] < 2; y[1][4]++)
                            {
                            for (y[4][1] = 0; y[4][1] < 2; y[4][1]++)
                              {
                              for (y[2][3] = 0; y[2][3] < 2; y[2][3]++)
                                {
                                for (y[3][2] = 0; y[3][2] < 2; y[3][2]++)
                                  {
                                  for (y[2][4] = 0; y[2][4] < 2; y[2][4]++)
                                    {
                                    for (y[4][2] = 0; y[4][2] < 2; y[4][2]++)
                                      {
                                      for (y[3][4] = 0; y[3][4] < 2; y[3][4]++)
                                        {
                                        for (y[4][3] = 0; y[4][3] < 2; y[4][3]++)
                                          {
                                          edges = y[0][1] + y[1][0] + y[0][2] + y[2][0] + y[0][3] + y[3][0] + y[0][4] + y[4][0] + y[1][2] + y[2][1] + y[1][3] + y[3][1] + y[1][4] + y[4][1] + y[2][3] + y[3][2] + y[2][4] + y[4][2] + y[3][4] + y[4][3];
                                          number_edges = &edges;
                                          heads = (int*) calloc(*number_edges,sizeof(int));
                                          tails = (int*) calloc(*number_edges,sizeof(int));
                                          count = -1;
                                          for (i = 0; i < n_k; i++)
                                            {
                                            for (j = 0; j < n_k; j++)
                                              {
                                              if (y[i][j] == 1) /* If directed, y[i][i] == 0 */ 
                                                {
                                                count = count + 1;
                                                heads[count] = block_members[i] + 1; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
                                                tails[count] = block_members[j] + 1;
                                                }
                                              }
                                            }
                                          b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); 
                                          sum = sum + e(b);
                                          free(heads);
                                          free(tails);
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }            
                        }   
                      }
                    }
                  }
                }    
              }  
            }
          }    
        }
      } 
    }
  a = ln(sum);
  for (i = 0; i < n_k; i++)
    {
    free(y[i]);
    }
  free(y);
  free(statistic);
  return a;
}

double Within_Block_Partition_Function(int model, latentstructure *ls, int k, ergmstructure *ergm, int *heads, int *tails, double *input, int *n, int *directed, int *number_terms, char **funnames, char **sonames)
/*
input: node i, latent structure, ergm structure
output: within-block partition function on log scale, evaluated either by complete enumeration or lower bounded by variational methods / mean-field methods
*/
{
  int zero = 0;
  int i, j, *bipartite, *block_members, count, variational_verbose = 0;
  double *eta, lower_bound_k, *theta;
  bipartite = &zero;
  if (ls->size[k] < 2) lower_bound_k = 0; /* No possible edge */
  else /* Possible edges */
    {       
    eta = (double*) calloc(2,sizeof(double));
    if (eta == NULL) { Rprintf("\n\ncalloc failed: Within_Block_Partition_Function, eta\n\n"); error("Error: out of memory"); }
    if ((ergm->d1 == 0) && (ergm->d2 == 2))
      {
      eta[0] = ls->theta[0][k]; /* Within-block edge parameter */
      eta[1] = ls->theta[1][k]; /* Within-block parameter */
      }  
    else if ((ergm->d1 == 1) && (ergm->d2 == 1))
      {
      eta[0] = ergm->theta[0]; /* Between- and within-block edge parameter */
      eta[1] = ls->theta[0][k]; /* Within-block parameter */
      }
    else if ((ergm->d1 == 1) && (ergm->d2 == 2))
      {
      eta[0] = ergm->theta[0] + ls->theta[0][k]; /* Between- and within-block edge parameter */
      eta[1] = ls->theta[1][k]; /* Within-block parameter */
      }
    Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input); /* Set input given ls->theta */
    theta = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter */
    /*
    Rprintf("\nn = %i eta = (%6.3f,%6.3f)",ls->size[k],eta[0],eta[1]);     
    */
    if (ls->size[k] < ls->threshold) /* Computation of within-block partition function by complete enumeration feasible */ 
      {
      block_members = NULL;
      count = -1;
      for (i = 0; i < ls->n; i++)
        {
        if (ls->indicator[i] == k) 
          {
          count = count + 1;            
          block_members = (int*) realloc(block_members,(count+1)*sizeof(int));
          block_members[count] = i;
          }
        }
      if (*directed == 0)
        {
        if (ls->size[k] == 2) lower_bound_k = Within_Block_Partition_Function_2_Graph(ls,block_members,ergm,input,theta,n,directed,bipartite,number_terms,funnames,sonames);
        else if (ls->size[k] == 3) lower_bound_k = Within_Block_Partition_Function_3_Graph(ls,block_members,ergm,input,theta,n,directed,bipartite,number_terms,funnames,sonames);
        else if (ls->size[k] == 4) lower_bound_k = Within_Block_Partition_Function_4_Graph(ls,block_members,ergm,input,theta,n,directed,bipartite,number_terms,funnames,sonames);
        else if (ls->size[k] == 5) lower_bound_k = Within_Block_Partition_Function_5_Graph(ls,block_members,ergm,input,theta,n,directed,bipartite,number_terms,funnames,sonames);
        /*
        else if (ls->size[k] == 6) lower_bound_k = Within_Block_Partition_Function_6_Graph(ls,block_members,ergm,input,theta,n,directed,bipartite,number_terms,funnames,sonames);
        else if (ls->size[k] == 7) lower_bound_k = Within_Block_Partition_Function_7_Graph(ls,block_members,ergm,input,theta,n,directed,bipartite,number_terms,funnames,sonames);
        */
        }
      else 
        {
        if (ls->size[k] == 2) lower_bound_k = Within_Block_Partition_Function_2_Digraph(ls,block_members,ergm,input,theta,n,directed,bipartite,number_terms,funnames,sonames);
        else if (ls->size[k] == 3) lower_bound_k = Within_Block_Partition_Function_3_Digraph(ls,block_members,ergm,input,theta,n,directed,bipartite,number_terms,funnames,sonames);
        else if (ls->size[k] == 4) lower_bound_k = Within_Block_Partition_Function_4_Digraph(ls,block_members,ergm,input,theta,n,directed,bipartite,number_terms,funnames,sonames);
        else if (ls->size[k] == 5) lower_bound_k = Within_Block_Partition_Function_5_Digraph(ls,block_members,ergm,input,theta,n,directed,bipartite,number_terms,funnames,sonames);
        }
      /*
      Rprintf("\nExact A: %8.4f",lower_bound_k);
      */
      free(block_members);
      }
    else /* Computation of within-block partition function by complete enumeration infeasible: lower bound within-block partition function by variational methods */
      {
      lower_bound_k = Variational_EM(ls->size[k],model,eta,*directed,variational_verbose);
      /*
      Rprintf("\nLower bound of A: %8.4f",lower_bound_k);
      */
      }
    free(eta);
    free(theta);
    } 
  return lower_bound_k;
}

double Between_Block_Partition_Function(latentstructure *ls, ergmstructure *ergm, double *input, double *theta, 
                                         int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: between-block partition function on log scale, computed under the assumption of conditional edge-independence given latent structure
*/
{
  int one = 1;
  int two = 2;
  int i, j, *number_edges, *heads, *tails;
  double a, b, *statistic, sum;
  statistic = (double*) calloc(ergm->d,sizeof(double));
  /* 
  Note 1: if undirected graph and i < j, undirected edge (i, j) should be stored as (i, j)
  Note 2: i, j should be integers between 1 and n
  */  
  a = 0.0; /* Log partition function */
  for (i = 1; i < ls->n; i++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
    {
    for (j = i + 1; j < ls->n + 1; j++) 
      {  
      if (ls->indicator[i-1] != ls->indicator[j-1]) /* Edges (i, j) and (j, i) are between blocks, because i and j are not in the same block */
        {
        if (*directed == 0)
          {
          number_edges = &one;
          heads = &i; 
          tails = &j;
          b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
          a = a + ln(1.0 + e(b));
          }
        else
          {
          sum = 1.0; /* Both edges (i, j) and (j, i) are absent */
          number_edges = &one;
          heads = &i; /* Edge (i, j) present while edge (j, i) absent */
          tails = &j;
          b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
          sum = sum + e(b);
          number_edges = &one;
          heads = &j; /* Edge (j, i) present while edge (i, j) absent */
          tails = &i;
          b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
          sum = sum + e(b);
          number_edges = &two;
          heads = (int*) calloc(2,sizeof(int));
          tails = (int*) calloc(2,sizeof(int));
          heads[0] = i; /* Both edges (i, j) and (j, i) present */
          tails[0] = j;
          heads[1] = j;
          tails[1] = i; 
          b = Minus_Energy(ergm->d,input,theta,heads,tails,number_edges,n,directed,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
          sum = sum + e(b);
          a = a + ln(sum);          
          free(heads);
          free(tails);        
          }
        }
      }
    }
  free(statistic);
  return a;
}

double* Candidate_Generating_Distribution_Indicators_Dependence(int node, int model, latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input, int *n_edges, int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: node i, latent structure, ergm structure
output: candidate-generating distribution
*/
{
  int k, block, indicator;
  double a_i, energy, lower_bound, lower_bound_k, *lower_bound_k_present, maximum, *q_i, *statistic, sum, *theta;
  lower_bound_k_present = (double*) calloc(ls->number,sizeof(double));
  if (lower_bound_k_present == NULL) { Rprintf("\n\ncalloc failed: Candidate_Generating_Distribution_Indicators_Dependence, lower_bound_k_present\n\n"); error("Error: out of memory"); }
  q_i = (double*) calloc(ls->number,sizeof(double));
  if (q_i == NULL) { Rprintf("\n\ncalloc failed: Candidate_Generating_Distribution_Indicators_Dependence, q_i\n\n"); error("Error: out of memory"); }
  statistic = (double*) calloc(ls->number,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: Candidate_Generating_Distribution_Indicators_Dependence, statistic\n\n"); error("Error: out of memory"); }
  indicator = ls->indicator[node];
  ls->size[indicator] = ls->size[indicator] - 1;
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input); /* Set input given ls->theta */
  theta = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter */
  for (k = 0; k < ls->number; k++)
    {
    lower_bound_k_present[k] = Within_Block_Partition_Function(model,ls,k,ergm,heads,tails,input,n,directed,nterms,funnames,sonames)
;
    }
  maximum = -DBL_MAX;
  for (block = 0; block < ls->number; block++) /* Block */
    {
    ls->indicator[node] = block;
    ls->size[block] = ls->size[block] + 1;
    Set_Input_Indicator(ergm->terms,ergm->hierarchical,ls->number,ls->n,node,block,input); /* Set input given indicator; reset in Gibbs_Indicators_Independence */
    energy = Minus_Energy(ergm->d,input,theta,heads,tails,n_edges,n,directed,bipartite,nterms,funnames,sonames,statistic);
    /* Lower bound on log partition function */
    lower_bound = 0.0;
    /* Within-block log partition functions */
    /*
    Rprintf("\n- block %i of size %i",block,ls->size[block]);
    */
    for (k = 0; k < ls->number; k++)
      {
      if (k == block) lower_bound_k = Within_Block_Partition_Function(model,ls,k,ergm,heads,tails,input,n,directed,nterms,funnames,sonames)
;
      else lower_bound_k = lower_bound_k_present[k];
      lower_bound = lower_bound + lower_bound_k;
      /*
      Rprintf(" %8.4f",lower_bound_k);
      */
      }
    /* Between-block log partition functions */
    lower_bound = lower_bound + Between_Block_Partition_Function(ls,ergm,input,theta,n,directed,bipartite,nterms,funnames,sonames);
    q_i[block] = energy - lower_bound;
    if (q_i[block] > maximum) maximum = q_i[block];
    ls->size[block] = ls->size[block] - 1; /* Reset size */
    }   
  sum = 0.0;
  for (k = 0; k < ls->number; k++) /* Translate mass to forgo underflow */
    {   
    q_i[k] = ln(ls->p[k]) + q_i[k] - maximum; 
    /*
    Rprintf(" %8.4f",q_i[k]);
    */
    sum = sum + e(q_i[k]);
    }
  a_i = ln(sum);
  /*
  Rprintf("\n- full conditional of node %i:",node+1);
  */
  for (k = 0; k < ls->number; k++) /* Probability mass */
    {
    q_i[k] = e(q_i[k] - a_i);
    /*  
    Rprintf(" %8.4f",q_i[k]);
    */
    }
  ls->indicator[node] = indicator; /* Reset indicator */
  ls->size[indicator] = ls->size[indicator] + 1; /* Reset size */
  free(lower_bound_k_present);
  free(statistic);
  free(theta);
  return q_i;
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
                          int *mheads, int *mtails, int *mdnedges, int *status)
/*
input: (maximum) number of categories, number of nodes, number of structural parameters, number of parameters
output: one sample from posterior predictive distribution
*/
{
  int number_networks, *h, *t, i, *indicator, k, *nedges, s;
  double **parameter;
  s = 1; /* Sample one graph from posterior predictive distribution */
  for (i = 0; i < ergm_d; i++)
    {
    sample[i] = 0.0;
    }
  number_networks = 1;
  MCMC_wrapper(&number_networks,dnedges,tails,heads, /* Sample one graph from posterior predictive distribution given input and theta */
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
               status);
  indicator = (int*) calloc(n,sizeof(int));
  if (indicator == NULL) { Rprintf("\n\ncalloc failed: Sample_Graph, indicator\n\n"); error("Error: out of memory"); }
  for (i = 0; i < n; i++) /* Identical indicators */
    {
    indicator[i] = 1;
    }
  parameter = (double**) calloc(ls_d,sizeof(double*));
  if (parameter == NULL) { Rprintf("\n\ncalloc failed: Sample_Graph, parameter\n\n"); error("Error: out of memory"); }
  for (i = 0; i < ls_d; i++)
    {
    parameter[i] = (double*) calloc(number+1,sizeof(double));
    if (parameter[i] == NULL) { Rprintf("\n\ncalloc failed: Sample_Graph, parameter[%i]\n\n",i); error("Error: out of memory"); }
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
  if (h == NULL) { Rprintf("\n\ncalloc failed: Sample_Graph, h\n\n"); error("Error: out of memory"); }
  t = (int*) calloc(*nedges,sizeof(int));
  if (t == NULL) { Rprintf("\n\ncalloc failed: Sample_Graph, t\n\n"); error("Error: out of memory"); }
  for (i = 0; i < *nedges; i++) /* Since first element of newnetworkheads and newnetworktails is number of edges, heads and tails must be extracted */
    {
    h[i] = newnetworkheads[i+1];
    t[i] = newnetworktails[i+1];
    }
  int timings = 0, time = 0, lasttoggle = 0;
  /*
  Rprintf("\n\nh_ergm_mcmc.c:");
  Rprintf("\ntimings=%p",&timings);
  Rprintf("\ntimings=%i",timings);
  */
  network_stats_wrapper(t,h,&timings,&time,&lasttoggle,nedges,dn,directed,bipartite,nterms,funnames,sonames,input,statistic); /* Compute non-structural function of graph */
  free(h);
  free(t);
  free(indicator);
  for (i = 0; i < ls_d; i++)
    {
    free(parameter[i]);
    }
  free(parameter);
}

