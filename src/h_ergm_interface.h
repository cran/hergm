#include "MCMC.h"
#include "netstats.h"
#include "h_ergm_basics.h"

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

void Set_Parameter(int d, int *structural, double *theta, double *parameter);
/*
input: number of ergm terms, indicator of structural parameters
output: parameter
*/

void Get_Parameter(int d, int max_number, int *structural, double *parameter, double *ergm_theta, double **ls_theta);
/*
input: number of ergm terms, number of categories, indicator of structural parameters, initial parameter specified by user
output: non-structural, structural parameters
*/

void MCMC_wrapper_h (int *heads, int *tails, int *dnedges,
                   int *maxpossibleedges,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, int *samplesize, 
                   double *sample, int *burnin, int *interval,  
                   int *newnetworkheads, 
                   int *newnetworktails, 
                   int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   int *maxedges,
                   int *mheads, int *mtails, int *mdnedges,
                   double *inputs_h, double *sample_h);

void MCMCSample_h (char *MHproposaltype, char *MHproposalpackage,
  double *theta, double *networkstatistics, 
  long int samplesize, long int burnin, 
  long int interval, int hammingterm, int fVerbose,
  Network *nwp, Model *m, DegreeBound *bd,
  double *input,
  int nmax,
  int *newnetworkheads, int *newnetworktails, int *n_nodes, int *dflag, int *bipartite,
  int *n_terms, char **funnames, char **sonames, double *input_h, double *networkstatistics_h);

void MCMCSample_h_2 (char *MHproposaltype, char *MHproposalpackage,
  double *theta, double *networkstatistics, 
  long int samplesize, long int burnin, 
  long int interval, int hammingterm, int fVerbose,
  Network *nwp, Model *m, DegreeBound *bd,
  double *input,
  int nmax,
  int *newnetworkheads, int *newnetworktails, int *n_nodes, int *dflag, int *bipartite,
  int *n_terms, char **funnames, char **sonames, int n_number, int n_input_h, double **inputs, double **networkstatistics_h);

void MCMC_wrapper_h_2 (int *heads, int *tails, int *dnedges,
                   int *maxpossibleedges,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   double *inputs, double *theta0, int *samplesize, 
                   double *sample, int *burnin, int *interval,  
                   int *newnetworkheads, 
                   int *newnetworktails, 
                   int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   int *maxedges,
                   int *mheads, int *mtails, int *mdnedges,
                   int n_number, int n_inputs_h, double **inputs_h, double **sample_h);

void Sample_Graph(int number, int n, int ls_d, int terms, int *hierarchical, int ergm_d, double *statistic,
                          int *heads, int *tails, int *dnedges,
                          int *maxpossibleedges,
                          int *dn, int *dflag, int *bipartite, 
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
		       int *n, int *dflag,  int *bipartite,
		       int *nterms, char **funnames,
		       char **sonames,
                       double *statistic);
/*
input: number of parameters, input parameters, parameters
output: statistic, inner product <parameter, statistic>
*/

double Ratio_Partition_Functions_1(int s, int d, double sum_observed, double *statistic_generating, double *statistic, double *theta_generating, double *theta);
/*
input: sample size, dimension, difference of inner products under alternative, data-generating parameter for observed graph, value of statistic under data-generating, alternative parameter, value of data-generating, alternative parameter
output: ratio of normalizing functions of ergms under alternative and data-generating parameter on log scale
*/

double Ratio_Partition_Functions_2(int s, int d, double sum_observed, double *statistic_generating, double *statistic, double *theta_generating, double *theta);
/*
input: sample size, dimension, difference of inner products under alternative, data-generating parameter for observed graph, value of statistic under data-generating, alternative parameter, value of data-generating, alternative parameter
output: ratio of partition functions of ergms under alternative and data-generating parameter on log scale
*/

double Ratio_Ergm_Pmfs(int *heads, int *tails, int *dnedges,
                   int *maxpossibleedges,
                   int *dn, int *dflag, int *bipartite, 
                   int *nterms, char **funnames,
                   char **sonames, 
                   char **MHproposaltype, char **MHproposalpackage,
                   int *samplesize, 
                   int *burnin, int *interval,  
                   int *newnetworkheads, 
                   int *newnetworktails, 
                   int *fVerbose, 
                   int *attribs, int *maxout, int *maxin, int *minout,
                   int *minin, int *condAllDegExact, int *attriblength, 
                   int *maxedges,
                   int *mheads, int *mtails, int *mdnedges,
                   double *input_proposal, double *input_present, int d, double *theta_proposal, double *theta_present);
/*
input: two values of input parameters, two values of parameter
output: ratio of ergm pmfs on log scale
*/

void Full_Conditional_Indicator(int *heads, int *tails, int *dnedges,
                        int *maxpossibleedges,
                        int *dn, int *dflag, int *bipartite, 
                        int *nterms, char **funnames,
                        char **sonames, 
                        char **MHproposaltype, char **MHproposalpackage,
                        int *samplesize, 
                        int *burnin, int *interval,  
                        int *newnetworkheads, 
                        int *newnetworktails, 
                        int *fVerbose, 
                        int *attribs, int *maxout, int *maxin, int *minout,
                        int *minin, int *condAllDegExact, int *attriblength, 
                        int *maxedges,
                        int *mheads, int *mtails, int *mdnedges,
                        int number, int n_input, double **input, int d, double *theta, int generating, double *p);
/*
input: number of categories, number of input parameters, dimension of parameter, category indicator of given node to be used to generate importance sample 
output: full conditional of category indicator of given node (given node is specified by calling function and input is set accordingly) 
*/

