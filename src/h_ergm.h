#include "h_ergm_latent.h"
#include "h_ergm_bayes.h"
#include "h_ergm_interface.h"
#include "h_ergm_initialize.c"

double Sample_Alpha(priorstructure_ls *prior_ls, latentstructure *ls);
/*
input: prior, latent structure
output: clustering parameter
*/

double* Stick_Breaking(double *shape1, double *shape2, latentstructure *ls);
/*
input: shape parameters of Beta distribution, latent structure
output: category probability vector
*/

double* Sample_P(latentstructure *ls);
/*
input: latent structure
output: category probability vector
*/

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

void Gibbs_Indicators_Independence(latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input_proposal, 
                       int *n_edges, int *n, int *directed, int *bipartite, int *nterms, char **funnames, char **sonames, double *q_i);
/*
input: latent structure, ergm structure
output: indicators
note: function more efficient than sister function Gibbs_Indicators_Independence
*/

int Sample_Ergm_Theta_Independence(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges, int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames, char **sonames, 
                        double *input, int print, int n_between, double *scale_factor);
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/

int Sample_Ls_Theta_Independence(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges, int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames, char **sonames, double *input_proposal, double *input_present, int print, int n_between, double *scale_factor);
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/

int Sample_Parameters_Independence(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges, int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames, char **sonames, 
                        double *input_proposal, double *input_present, int print, int n_between, double *scale_factor);
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/

int Sample_Parameters_Dependence(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges,
                        int *maxpossibleedges,
                        int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames,
                        char **sonames, 
                        char **MHproposaltype, char **MHproposalpackage,
                        double *sample,
                        int *burnin, int *interval,  
                        int *verbose, 
                        int *attribs, int *maxout, int *maxin, int *minout,
                        int *minin, int *condAllDegExact, int *attriblength, 
                        int *maxedges,
                        int *mheads, int *mtails, int *mdnedges,
                        double *input_present, int print,
                        int *newnetworkheads, int *newnetworktails, int n_between, double *scale_factor, double *q_i);
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/

void Gibbs_Parameters(ergmstructure *ergm, latentstructure *ls, priorstructure *prior);
/*
input: ergm structure, latent structure, prior
output: non-structural parameters not showing up in the ergm pmf
*/

double* Gibbs_Parameters_Means(priorstructure *prior, latentstructure *ls);
/*
input: prior structure, latent structure
output: means of parameters
*/

double* Gibbs_Parameters_Precisions(priorstructure *prior, latentstructure *ls);
/*
input: prior structure, latent structure
output: precisions of parameters
*/

void Initial_State(int *parallel, double *alpha, int *indicator, priorstructure_ls *prior_ls, priorstructure *prior, latentstructure *ls, ergmstructure *ergm, double *theta, double *scale_factor);
/* 
input: clustering parameter, priors, latent structure, ergm structure, user-specified initial value of non-structural parameters
*/

int Sample_CRP(latentstructure *ls);
/*
input: latent structure ls
output: partition of set of nodes drawn from Chinese restaurant process with scaling parameter ls->alpha
*/

int Sample_Graph_Edge_Independence(int *directed, latentstructure *ls, double *ln_p, int *heads, int *tails);
/*
input: latent structure; probability of edge between nodes i and j on log scale
output: graph sampled from PMF p and number of edges
*/

void Simulation(int *dyaddependence,
             int *hierarchical,
             int *d, 
             int *d1, 
             int *d2,
             int *structural,
             int *min_size,
             int *max_number,
             double *alpha,
             double *alpha_shape,
             double *alpha_rate,
             double *m1,
             double *m2,
             double *b,
             double *cf1,
             double *cf2,
             double *p1,
             double *p2,
             double *m2_mean,
             double *m2_precision,
             double *p2_shape,
             double *p2_rate,
             double *eta,
             int *indicator,
             int *heads, int *tails, int *dnedges,
             int *maxpossibleedges,
             int *dn, int *directed, int *bipartite, 
             int *nterms, char **funnames,
             char **sonames, 
             char **MHproposaltype, char **MHproposalpackage,
             double *inputs, double *inputs_h, double *theta, int *samplesize, 
             double *sample,
             int *burnin, int *interval,  
             int *v, 
             int *attribs, int *maxout, int *maxin, int *minout,
             int *minin, int *condAllDegExact, int *attriblength, 
             int *maxedges,
             int *max_iterations, int *n_between_block_parameters, int *output, double *mcmc, int *sample_heads, int *sample_tails, int *call_RNGstate, int *hyperprior);
/*
input: R input
output: simulated graph
*/

void Inference(int *dyaddependence,
             int *hierarchical,
             int *d, 
             int *d1, 
             int *d2,
             int *structural,
             int *min_size,
             int *max_number,
             double *alpha,
             double *alpha_shape,
             double *alpha_rate,
             double *m1,
             double *m2,
             double *b,
             double *cf1,
             double *cf2,
             double *p1,
             double *p2,
             double *m2_mean,
             double *m2_precision,
             double *p2_shape,
             double *p2_rate,
             int *indicator,
             int *heads, int *tails, int *dnedges,
             int *maxpossibleedges,
             int *dn, int *directed, int *bipartite, 
             int *nterms, char **funnames,
             char **sonames, 
             char **MHproposaltype, char **MHproposalpackage,
             double *inputs, double *inputs_h, double *theta, int *samplesize, 
             double *sample,
             int *burnin, int *interval,  
             int *newnetworkheads, 
             int *newnetworktails, 
             int *v, 
             int *attribs, int *maxout, int *maxin, int *minout,
             int *minin, int *condAllDegExact, int *attriblength, 
             int *maxedges,
             int *max_iterations, int *n_between_block_parameters, int *output, double *mcmc, double *scalefactor, double *mh_accept, double *q_i, int *call_RNGstate, int *parallel, int *hyperprior);
/*
input: R input
output: MCMC sample of unknowns from posterior
*/

