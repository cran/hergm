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
#include "signal.h"

/*
Optimization algorithms based on evaluations of the objective function or its gradient:
#include <math.h>
#include <nlopt.h>
*/
 
#define Min(a,b) ((a)<(b) ? (a) : (b))

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

double* Stick_Breaking_External(double *shape1, double *shape2, int number, int n);
/*
input: shape parameters of Beta distribution, latent structure
output: category probability vector
*/

double* Sample_P(latentstructure *ls);
/*
input: latent structure
output: category probability vector
*/

int Sample_Ergm_Theta_Independence(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges, int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames, char **sonames, 
                        double *input, int print, double *scale_factor);
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/

int Sample_Ls_Theta_Independence(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges, int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames, char **sonames, double *input_proposal, double *input_present, int print, double *scale_factor);
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/

int Sample_Parameters_Independence(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges, int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames, char **sonames, 
                        double *input_proposal, double *input_present, int print, double *scale_factor);
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/

int Sample_Indicators_Dependence(int model, ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
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
                        int *newnetworkheads, int *newnetworktails, double *scale_factor, int update_node, double *temperature, int *status);
/*
input: ergm structure, latent structure, prior
output: indicators
*/

int Sample_Ergm_Theta_Dependence(int model, ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
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
                        int *newnetworkheads, int *newnetworktails, double *scale_factor, int *status);
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/

int Sample_Ls_Theta_Dependence(int model, ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
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
                        int *newnetworkheads, int *newnetworktails, double *scale_factor, int update_block, int *status);
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/

int Sample_Ls_Theta_Between(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges,
                        int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames, char **sonames, 
                        double *input_present, int print,
                        double *scale_factor);
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

double* Gibbs_Parameters_Means_Conditional(priorstructure *prior, latentstructure *ls);
/*
input: prior structure, latent structure
output: means of parameters
*/

double* Gibbs_Parameters_Precisions(priorstructure *prior, latentstructure *ls);
/*
input: prior structure, latent structure
output: precisions of parameters
*/

double* Gibbs_Parameters_Precisions_Marginal(priorstructure *prior, latentstructure *ls);
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

int Sample_Graph_Edge_Independence(latentstructure *ls, double *ln_p, int *heads, int *tails);
/*
input: latent structure; probability of edge between nodes i and j on log scale
output: graph sampled from PMF p and number of edges
*/

void Dirichlet(int *n,
               int *number,
               double *alpha,
               double *eta_mean,
               double *eta_sd,
               int *indicator,
               double *eta);
/*
input: R input
output: draw from truncated Dirichlet process prior:
- indicator of block membership
- parameters
*/

void Simulation(int *dyaddependence,
             int *hierarchical,
             double *scaling,
             int *d, 
             int *d1, 
             int *d2,
             int *structural,
             int *min_size,
             int *max_number,
             int *null_alpha,
             int *null_eta,
             int *null_indicator,
             double *alpha,
             double *alpha_shape,
             double *alpha_rate,
             double *m1,
             double *m2,
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
             int *max_iterations, int *between, double *mean_between, int *output, double *mcmc, int *sample_heads, int *sample_tails, int *prior_assumption, int *status);
/*
input: R input
output: simulated graph
*/

void Inference(int *model_type,
             int *dyaddependence,
             int *hierarchical,
             double *scaling,
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
             int *max_iterations, int *between, int *output, double *mcmc, double *scalefactor, double *mh_accept, double *q_i, int *parallel, double *temperature, int *prior_assumptions, int *status);
/*
input: R input
output: MCMC sample of unknowns from posterior
*/

