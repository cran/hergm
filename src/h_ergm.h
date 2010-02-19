#include "h_ergm_utils.h"
#include "h_ergm_latent.h"
#include "h_ergm_bayes.h"
#include "h_ergm_interface.h"

void Sample_Alpha(priorstructure_ls *prior_ls, latentstructure *ls);
/*
input: prior, latent structure
output: clustering parameter
*/

void Stick_Breaking(double *shape1, double *shape2, latentstructure *ls);
/*
input: shape parameters of Beta distribution, latent structure
output: category probability vector
*/

void Sample_P(latentstructure *ls);
/*
input: latent structure
output: category probability vector
*/

void Set_Input_i(int i, double *input_present, int *index, int n_input, double **input, ergmstructure *ergm, latentstructure *ls);
/*
input: node, input (vector), starting index of row of input (matrix), number of columns of input (matrix), input (matrix), ergm structure, latent structure
output: input where indicator of node is set to all possible values 
*/

void Gibbs_Indicator_i(int i, double *p, latentstructure *ls);
/* 
input: node, full conditional, latent structure
output: category indicator
*/

void Sample_Indicators_1(latentstructure *ls, ergmstructure *ergm,
                       int *heads, int *tails, int *dnedges,
                       int *maxpossibleedges,
                       int *dn, int *dflag, int *bipartite, 
                       int *nterms, char **funnames,
                       char **sonames, 
                       char **MHproposaltype, char **MHproposalpackage,
                       int *samplesize, 
                       int *burnin, int *interval,  
                       int *newnetworkheads, 
                       int *newnetworktails, 
                       int *verbose, 
                       int *attribs, int *maxout, int *maxin, int *minout,
                       int *minin, int *condAllDegExact, int *attriblength, 
                       int *maxedges,
                       int *mheads, int *mtails, int *mdnedges,
                       double *input_present, double *theta, int print);
/*
input: latent structure
output: category indicators
*/

void Gibbs_Indicators_1(int h, int i, int j, latentstructure *ls, ergmstructure *ergm,
                       int *heads, int *tails, int *dnedges,
                       int *maxpossibleedges,
                       int *dn, int *dflag, int *bipartite, 
                       int *nterms, char **funnames,
                       char **sonames, 
                       char **MHproposaltype, char **MHproposalpackage,
                       int *samplesize, 
                       int *burnin, int *interval,  
                       int *newnetworkheads, 
                       int *newnetworktails, 
                       int *verbose, 
                       int *attribs, int *maxout, int *maxin, int *minout,
                       int *minin, int *condAllDegExact, int *attriblength, 
                       int *maxedges,
                       int *mheads, int *mtails, int *mdnedges,
                       double *input_present, double *theta);
/*
input: node, latent structure
output: category indicators
*/

void Gibbs_Indicators_2(int i, int j, latentstructure *ls, ergmstructure *ergm,
                       int *heads, int *tails, int *dnedges,
                       int *maxpossibleedges,
                       int *dn, int *dflag, int *bipartite, 
                       int *nterms, char **funnames,
                       char **sonames, 
                       char **MHproposaltype, char **MHproposalpackage,
                       int *samplesize, 
                       int *burnin, int *interval,  
                       int *newnetworkheads, 
                       int *newnetworktails, 
                       int *verbose, 
                       int *attribs, int *maxout, int *maxin, int *minout,
                       int *minin, int *condAllDegExact, int *attriblength, 
                       int *maxedges,
                       int *mheads, int *mtails, int *mdnedges,
                       double *input_present, double *theta);
/*
input: two nodes, latent structure
output: category indicators
*/

void Gibbs_Indicators_3(int h, int i, int j, latentstructure *ls, ergmstructure *ergm,
                       int *heads, int *tails, int *dnedges,
                       int *maxpossibleedges,
                       int *dn, int *dflag, int *bipartite, 
                       int *nterms, char **funnames,
                       char **sonames, 
                       char **MHproposaltype, char **MHproposalpackage,
                       int *samplesize, 
                       int *burnin, int *interval,  
                       int *newnetworkheads, 
                       int *newnetworktails, 
                       int *verbose, 
                       int *attribs, int *maxout, int *maxin, int *minout,
                       int *minin, int *condAllDegExact, int *attriblength, 
                       int *maxedges,
                       int *mheads, int *mtails, int *mdnedges,
                       double *input_present, double *theta);
/*
input: three nodes, latent structure
output: category indicators
*/

int Sample_Indicators_2(latentstructure *ls, ergmstructure *ergm,
                       int *heads, int *tails, int *dnedges,
                       int *maxpossibleedges,
                       int *dn, int *dflag, int *bipartite, 
                       int *nterms, char **funnames,
                       char **sonames, 
                       char **MHproposaltype, char **MHproposalpackage,
                       int *samplesize, 
                       int *burnin, int *interval,  
                       int *newnetworkheads, 
                       int *newnetworktails, 
                       int *verbose, 
                       int *attribs, int *maxout, int *maxin, int *minout,
                       int *minin, int *condAllDegExact, int *attriblength, 
                       int *maxedges,
                       int *mheads, int *mtails, int *mdnedges,
                       double *input_proposal, double *input_present, double *theta, int print);
/*
input: latent structure
output: category indicators
*/

void P_Independence(int *number_terms, int *number_parameters, double *input, double *theta,  int *n, int *flag, int *bipartite, char **funnames, char **sonames, double *p);
/*
input: number of terms; number of parameters;  input vector; parameter vector; number of nodes; other variables
output: probabilities of edges between nodes i and j on log scale, computed under the assumption of conditional dyad-independence given latent structure,
and ordered in accordance with i < j
*/

double Partition_Function_Independence(latentstructure *ls, ergmstructure *ergm, double *input, double *theta, 
                                         int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function on log scale, computed under the assumption of conditional dyad-independence given latent structure
*/

double Partition_Function_Independence_Node(int node, latentstructure *ls, ergmstructure *ergm, double *input, double *theta, 
                                         int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: partition function of node i on log scale, computed under the assumption of conditional dyad-independence given latent structure
*/

double PMF_Independence(latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input, double *theta, 
                        int *n_edges, int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: probability mass on log scale, computed under the assumption of dyad-dependence
*/

double PMF_Independence_Node(int i, latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input, double *theta, 
                        int *n_edges, int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: input
output: mass of node i on log scale, computed under the assumption of dyad-dependence
*/

double PMF_i_k(int i, int l, latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input_proposal, 
                       int *n_edges, int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: node i, catogory l, latent structure, ergm structure
output: conditional PMF of graph given latent structure 
*/

double PMF_i_k_Node(int i, int l, latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input_proposal, 
                       int *n_edges, int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: node i, catogory l, latent structure, ergm structure
output: conditional PMF of graph given latent structure 
*/

void Gibbs_Indicators_Independence(latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input_proposal, 
                       int *n_edges, int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames);
/*
input: latent structure, ergm structure
output: indicators
note: function more efficient than sister function Gibbs_Indicators_Independence
*/

int Sample_Parameters_1(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges,
                        int *maxpossibleedges,
                        int *dn, int *dflag, int *bipartite, 
                        int *nterms, char **funnames,
                        char **sonames, 
                        char **MHproposaltype, char **MHproposalpackage,
                        int *samplesize, 
                        int *burnin, int *interval,  
                        int *newnetworkheads, 
                        int *newnetworktails, 
                        int *verbose, 
                        int *attribs, int *maxout, int *maxin, int *minout,
                        int *minin, int *condAllDegExact, int *attriblength, 
                        int *maxedges,
                        int *mheads, int *mtails, int *mdnedges,
                        double *input_proposal, double *input_present, int print, double scale_factor);
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/

void Sample_Parameters_2(ergmstructure *ergm, latentstructure *ls, priorstructure *prior);
/*
input: ergm structure, latent structure, prior
output: non-structural parameters not showing up in the ergm pmf
*/

int Sample_Parameters_3(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges,
                        int *maxpossibleedges,
                        int *dn, int *dflag, int *bipartite, 
                        int *nterms, char **funnames,
                        char **sonames, 
                        char **MHproposaltype, char **MHproposalpackage,
                        int *samplesize, 
                        int *burnin, int *interval,  
                        int *newnetworkheads, 
                        int *newnetworktails, 
                        int *verbose, 
                        int *attribs, int *maxout, int *maxin, int *minout,
                        int *minin, int *condAllDegExact, int *attriblength, 
                        int *maxedges,
                        int *mheads, int *mtails, int *mdnedges,
                        double *input_proposal, double *input_present, int print, int n_between, double scale_factor);
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/

int Sample_Parameters_Indicators(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges,
                        int *maxpossibleedges,
                        int *dn, int *dflag, int *bipartite, 
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
                        int *newnetworkheads, int *newnetworktails, int n_between, double scale_factor);
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/

void Initial_State(double *alpha, int *indicator, priorstructure_ls *prior_ls, priorstructure *prior, latentstructure *ls, ergmstructure *ergm, double *theta);
/* 
input: clustering parameter, priors, latent structure, ergm structure, user-specified initial value of non-structural parameters
*/

int Sample_Graph_Independence(latentstructure *ls, double *ln_p, int *heads, int *tails);
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
             int *indicator,
             int *heads, int *tails, int *dnedges,
             int *maxpossibleedges,
             int *dn, int *dflag, int *bipartite, 
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
             int *max_iterations, int *n_between_block_parameters, int *output, double *mcmc, int *sample_heads, int *sample_tails);
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
             int *indicator,
             int *heads, int *tails, int *dnedges,
             int *maxpossibleedges,
             int *dn, int *dflag, int *bipartite, 
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
             int *max_iterations, int *n_between_block_parameters, int *output, double *mcmc, double *scalefactor, double *mh_accept);
/*
input: R input
output: MCMC sample of unknowns from posterior
*/

