#include "MCMC.h"
#include "netstats.h"
#include "h_ergm_utils.h"

int Number_Input(int terms, double *input);
/*
input: number of ergm terms, input parameters
output: number of input parameters
*/

void Set_Input(int terms, int *hierarchical, int max_number, int n, int *indicator, double **theta, double *input);
/*
input: number of ergm terms, indicator of hierarchical ergm terms, number of nodes, node-bound category indicator, category-bound parameter
output: input parameters
*/

void Set_Input_Indicator(int terms, int *hierarchical, int max_number, int n, int node, int node_indicator, double *input);
/*
input: number of ergm terms, indicator of hierarchical ergm terms, number of nodes, node-bound category indicator, category-bound parameter
output: input parameters
*/

double* Get_Parameter(int d, int *structural, double *theta);
/*
input: number of ergm terms, indicator of structural parameters
output: parameter
*/

void Set_Parameter(int d, int *structural, double *theta, double *parameter);
/*
input: number of ergm terms, indicator of structural parameters
output: parameter
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

double Minus_Energy(int d, double *input, double *parameter, 
                       int *heads, int *tails, int *nedges, 
		       int *n, int *directed,  int *bipartite,
		       int *nterms, char **funnames,
		       char **sonames,
                       double *statistic);
/*
input: number of parameters, input parameters, parameters
output: statistic, inner product <parameter, statistic>
*/

double Ratio_Partition_Functions(int s, int d, double sum_observed, double *statistic_generating, double *statistic, double *theta_generating, double *theta);
/*
input: sample size, dimension, difference of inner products under alternative, data-generating parameter for observed graph, value of statistic under data-generating, alternative parameter, value of data-generating, alternative parameter
output: ratio of partition functions of ergms under alternative and data-generating parameter on log scale
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

