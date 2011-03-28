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
#include "h_ergm_bayes.h"
#include "h_ergm_variational.h"
#include "h_ergm_interface.h"

void P_Edge_Independence(int *number_terms, int *number_parameters, double *input, double *theta,  int *n, int *directed, int *bipartite, char **funnames, char **sonames, double *p);
/*
input: directed graph; number of terms; number of parameters;  input vector; parameter vector; number of nodes; other variables
output: probabilities of edges between nodes i and j on log scale, computed under the assumption of conditional edge-independence given latent structure,
and ordered in accordance with i < j
*/

double Partition_Function_Edge_Independence(latentstructure *ls, ergmstructure *ergm, double *input, double *theta, 
                                         int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale, computed under the assumption of conditional edge-independence given latent structure
*/

double Partition_Function_Dyad_Independence(latentstructure *ls, ergmstructure *ergm, double *input, double *theta, 
                                         int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale, computed under the assumption of conditional dyad-independence given latent structure
*/

double PMF_Independence(latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input, double *theta, 
                        int *n_edges, int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: probability mass on log scale, computed under the assumption of dyad-dependence
*/

double PMF_i_k_Node(int i, int l, latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input_proposal, 
                       int *n_edges, int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: node i, catogory l, latent structure, ergm structure
output: conditional PMF of graph given latent structure 
*/

double PMF_Edge_Independence_Node(int node, int d, double *input, double *theta, 
                             int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames, int *n_edges, int *heads, int *tails);
/*
input: input
output: minus energy of node i on log scale, computed under the assumption of conditional dyad-independence given latent structure
*/

double PMF_Dyad_Independence_Node(int node, int d, double *input, double *theta, 
                             int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames, int *n_edges, int *heads, int *tails);
/*
input: input
output: minus energy of node i on log scale, computed under the assumption of conditional dyad-independence given latent structure
*/

int* Degree_Sequence(int n, int directed, int n_edges, int *heads, int *tails);
/*
input: number of nodes, indicator of directed network, number of edges, heads and tails of edge list
output: degree sequence
*/

int* Degree_Freq(int n, int* degree);
/*
input: number of nodes, degree sequence
output: degree frequencies
*/

double* Block_Degree_Freq(int n, int *degree, int block_number, int *block_size, int *block_indicator);
/*
input: number of nodes, degree sequence, number of blocks, size of blocks, indicator of block membership
output: relative frequencies of degree by block
*/

double Within_Block_Partition_Function_2_Graph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale given n = 2 nodes and undirected graphs
*/

double Within_Block_Partition_Function_3_Graph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale given n = 3 nodes and undirected graphs
*/

double Within_Block_Partition_Function_4_Graph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale given n = 4 nodes and undirected graphs
*/

double Within_Block_Partition_Function_5_Graph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale given n = 5 nodes and undirected graphs
*/

double Within_Block_Partition_Function_6_Graph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale given n = 6 nodes and undirected graphs
*/

double Within_Block_Partition_Function_7_Graph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale given n = 7 nodes and undirected graphs
*/

double Within_Block_Partition_Function_2_Digraph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale given n = 2 nodes and directed graphs
*/

double Within_Block_Partition_Function_3_Digraph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale given n = 3 nodes and directed graphs
*/

double Within_Block_Partition_Function_4_Digraph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale given n = 4 nodes and directed graphs
*/

double Within_Block_Partition_Function_5_Digraph(latentstructure *ls, int *block_members, ergmstructure *ergm, double *input, double *theta, 
           int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale given n = 5 nodes and directed graphs
*/

double Within_Block_Partition_Function(int model, latentstructure *ls, int k, ergmstructure *ergm, int *heads, int *tails, double *input, int *n, int *directed, int *number_terms, char **funnames, char **sonames);
/*
input: node i, latent structure, ergm structure
output: within-block partition function on log scale, evaluated either by complete enumeration or lower bounded by variational methods / mean-field methods
*/


double Between_Block_Partition_Function(latentstructure *ls, ergmstructure *ergm, double *input, double *theta, 
                                         int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: between-block partition function on log scale, computed under the assumption of conditional edge-independence given latent structure
*/

double* Candidate_Generating_Distribution_Indicators_Dependence(int node, int model, latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input, int *n_edges, int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: node i, latent structure, ergm structure
output: candidate-generating distribution
*/

double Ratio_Partition_Functions(int s, int d, double sum_observed, double *statistic_generating, double *statistic, double *theta_generating, double *theta);
/*
input: sample size, dimension, difference of inner products under alternative, data-generating parameter for observed graph, value of statistic under data-generating, alternative parameter, value of data-generating, alternative parameter
output: ratio of partition functions of ergms under alternative and data-generating parameter on log scale
*/

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
                          int *mheads, int *mtails, int *mdnedges);
/*
input: (maximum) number of categories, number of nodes, number of structural parameters, number of parameters
output: one sample from posterior predictive distribution
*/

