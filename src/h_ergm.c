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

#include "h_ergm.h"

double Sample_Alpha(priorstructure_ls *prior_ls, latentstructure *ls)
/*
input: prior, latent structure
output: clustering parameter
*/
{
  double alpha, rate, shape;
  shape = prior_ls->alpha_shape + (ls->number - 1.0); /* Shape of full conditional Gamma */
  rate = prior_ls->alpha_rate - ln(ls->p[ls->number-1]); /* Rate (inverse scale) of full conditional Gamma */
  alpha = rgamma(shape,1.0/rate); 
  return alpha;
}

double* Stick_Breaking(double *shape1, double *shape2, latentstructure *ls)
/*
input: shape parameters of Beta distribution, latent structure
output: category probability vector
*/
{
  int i;
  double *b, c, *p;
  b = (double*) calloc(ls->number,sizeof(double));
  if (b == NULL) { Rprintf("\n\ncalloc failed: Stick_Breaking, b\n\n"); error("Error: out of memory"); }
  p = (double*) calloc(ls->number,sizeof(double));
  if (p == NULL) { Rprintf("\n\ncalloc failed: Stick_Breaking, p\n\n"); error("Error: out of memory"); }
  /* Sample beta random variates: */
  /*
  Rprintf("\nStick_Breaking");
  */
  for (i = 0; i < (ls->number - 1); i++)
    {
    b[i] = rbeta(shape1[i],shape2[i]); 
    /*
    Rprintf("\nshape1[%i] = %f shape2[%i] = %f",i,shape1[i],i,shape2[i]);
    */
    }
  b[ls->number-1] = 1.0; /* Ensure that probabilities sum to one */
  /* p as function of b: */
  p[0] = b[0];
  /*
  Rprintf("\nb[%i] = %f, p[%i] = %f",0,b[0],0,p[0]);
  */
  c = 1.0;
  for (i = 1; i < ls->number; i++)
    {
    c = c * (1.0 - b[i-1]);
    p[i] = b[i] * c; 
    /*
    Rprintf("\nb[%i] = %f c = %f p[%i] = %f",i,b[i],c,i,p[i]);
    */
    }
  /*
  c = 0;
  for (i = 0; i < ls->number; i++)
    {
    c = c + p[i];
    }
  Rprintf("\nlength of p = %f",c);
  */
  free(b);
  return p;
}

double* Stick_Breaking_External(double *shape1, double *shape2, int number, int n)
/*
input: shape parameters of Beta distribution, latent structure
output: category probability vector
*/
{
  int i;
  double *b, c, *p;
  b = (double*) calloc(number,sizeof(double));
  if (b == NULL) { Rprintf("\n\ncalloc failed: Stick_Breaking, b\n\n"); error("Error: out of memory"); }
  p = (double*) calloc(number,sizeof(double));
  if (p == NULL) { Rprintf("\n\ncalloc failed: Stick_Breaking, p\n\n"); error("Error: out of memory"); }
  /* Sample beta random variates: */
  /*
  Rprintf("\nStick_Breaking");
  */
  for (i = 0; i < (number - 1); i++)
    {
    b[i] = rbeta(shape1[i],shape2[i]); 
    /*
    Rprintf("\nshape1[%i] = %f shape2[%i] = %f",i,shape1[i],i,shape2[i]);
    */
    }
  b[number-1] = 1.0; /* Ensure that probabilities sum to one */
  /* p as function of b: */
  p[0] = b[0];
  /*
  Rprintf("\nb[%i] = %f, p[%i] = %f",0,b[0],0,p[0]);
  */
  c = 1.0;
  for (i = 1; i < number; i++)
    {
    c = c * (1.0 - b[i-1]);
    p[i] = b[i] * c; 
    /*
    Rprintf("\nb[%i] = %f c = %f p[%i] = %f",i,b[i],c,i,p[i]);
    */
    }
  /*
  c = 0;
  for (i = 0; i < number; i++)
    {
    c = c + p[i];
    }
  Rprintf("\nlength of p = %f",c);
  */
  free(b);
  return p;
}

double* Sample_P(latentstructure *ls)
/*
input: latent structure
output: category probability vector
*/
{
  int i, rest;
  double *p, *shape1, *shape2;
  shape1 = (double*) calloc(ls->number-1,sizeof(double)); /* Element 0..ls->number-2 required */
  if (shape1 == NULL) { Rprintf("\n\ncalloc failed: Sample_P, sample1\n\n"); error("Error: out of memory"); }
  shape2 = (double*) calloc(ls->number-1,sizeof(double)); /* Element 0..ls->number-2 required */
  if (shape2 == NULL) { Rprintf("\n\ncalloc failed: Sample_P, sample2\n\n"); error("Error: out of memory"); }
  rest = ls->n;
  for (i = 0; i < (ls->number - 1); i++)
    {
    rest = rest - ls->size[i]; /* Number of nodes in category i + 1, ..., ls->number */
    shape1[i] = 1.0 + ls->size[i]; /* First shape parameter of Beta distribution */
    shape2[i] = ls->alpha + rest; /* Second shape parameter of Beta distribution */
    }
  p = Stick_Breaking(shape1,shape2,ls); /* Construct category probability vector by stick-breaking */
  free(shape1);
  free(shape2);
  return p;
}

int Sample_Ergm_Theta_Independence(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges, int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames, char **sonames, 
                        double *input, int print, double *scale_factor)
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/
{
  int accept, i;
  double **cf, *ergm_theta, log_present, log_proposal, log_ratio, *theta_present, *theta_proposal;
  log_ratio = 0.0;
  /* Propose ergm->theta: random walk (ratio of proposal pdfs cancels): */
  cf = Scale(ergm->d1,ergm->d1,prior->cf1,scale_factor[0]); /* Rescale Cholesky factor of Gaussian prior */        
  ergm_theta = Sample_MVN(ergm->d1,ergm->theta,cf); /* Random walk Metropolis-Hastings */
  log_proposal = MVN_PDF(ergm->d1,ergm_theta,prior->mean1,prior->precision1); /* Prior pdf: proposal */
  log_present = MVN_PDF(ergm->d1,ergm->theta,prior->mean1,prior->precision1); /* Prior pdf: present */
  log_ratio = log_ratio + (log_proposal - log_present);
  /*
  Rprintf("\n- log ratio (parameters) = %8.4f",log_ratio);  
  */
  /* Decide: */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input); /* Set input given ls->theta */
  theta_proposal = Get_Parameter(ergm->d,ergm->structural,ergm_theta); /* Set parameter */
  theta_present = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter */  
  log_proposal = PMF_Independence(ls,ergm,heads,tails,input,theta_proposal,dnedges,dn,directed,bipartite,nterms,funnames,sonames); /* Probability mass under proposed parameter */
  log_present = PMF_Independence(ls,ergm,heads,tails,input,theta_present,dnedges,dn,directed,bipartite,nterms,funnames,sonames); /* Probability mass under present parameters */
  log_ratio = log_ratio + (log_proposal - log_present);
  accept = MH_Decision(log_ratio);
  if (accept == 1) /* Proposal accepted: set ergm->theta to ergm_theta */
    {
    Set_D_D(ergm->d1,ergm->theta,ergm_theta);
    }
  if (print >= 1)
    {
    Rprintf("\nSample parameters:");
    Rprintf("\n- log ratio: %8.4f",log_ratio);  
    Rprintf("\n- decision: %i",accept);
    }
  free(ergm_theta);
  free(theta_present);
  free(theta_proposal);
  for (i = 0; i < ergm->d1; i++)
    { 
    free(cf[i]);
    }
  free(cf);
  return accept;
}

int Sample_Ls_Theta_Independence(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges, int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames, char **sonames, double *input_proposal, double *input_present, int print, double *scale_factor)
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/
{
  int accept, i;
  double **cf, *present, log_present, log_proposal, log_ratio, **ls_theta, *proposal, *theta;
  /* Proposal:
  note 1: all ls->theta such that ls->size >= ls->minimum_size and all ergm->theta are updated by random walk Metropolis-Hastings algorithm
  note 2: ratio of proposal pdfs cancels under random walk Metropolis-Hastings algorithm */
  log_ratio = 0.0;
  /* Propose ls->theta: */
  ls_theta = (double**) calloc(ls->d,sizeof(double*)); /* Remark: since memory is allocated to ls_theta by calloc, between-block parameters are 0 */
  if (ls_theta == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Independence, ls_theta\n\n"); error("Error: out of memory"); }
  for (i = 0; i < ls->d; i++)
    {
    ls_theta[i] = (double*) calloc(ls->number+1,sizeof(double));
    if (ls_theta[i] == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Independence, ls_theta[%i]\n\n",i); error("Error: out of memory"); }
    }
  present = (double*) calloc(ls->d,sizeof(double));  
  if (present == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Independence, present\n\n"); error("Error: out of memory"); }
  cf = Scale(ls->d,ls->d,prior->cf2,scale_factor[1]); /* Rescale Cholesky factor of Gaussian prior */ 
  for (i = 0; i < ls->number; i++) 
    {
    Get_Column(ls->d,present,ls->theta,i); /* Set mean to ls->theta[][i] */
    if (ls->size[i] < ls->minimum_size) Set_Column(ls->d,ls_theta,i,present); /* Set proposal = present */ 
    else 
      {
      /* Generate candidate: */
      proposal = Sample_MVN(ls->d,present,cf); /* Random walk Metropolis-Hastings algorithm */
      Set_Column(ls->d,ls_theta,i,proposal); /* Set ls_theta[][i] to proposal */
      /* Add ratio of prior pdf: */
      log_proposal = MVN_PDF(ls->d,proposal,prior->mean2,prior->precision2); /* Prior pdf of proposal */
      log_present = MVN_PDF(ls->d,present,prior->mean2,prior->precision2); /* Prior pdf of present */
      log_ratio = log_ratio + (log_proposal - log_present);
      free(proposal);
      }
    }
  for (i = 0; i < ls->d; i++) /* Set between-block parameters */
    {
    ls_theta[i][ls->number] = ls->theta[i][ls->number];
    }
  /* Decide: */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls_theta,input_proposal); /* Set input given ls_theta */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input given ls->theta */
  theta = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter; note: if ergm->d1 == 0, ergm_theta is not used */  
  log_proposal = PMF_Independence(ls,ergm,heads,tails,input_proposal,theta,dnedges,dn,directed,bipartite,nterms,funnames,sonames); /* Probability mass under proposed parameter */
  log_present = PMF_Independence(ls,ergm,heads,tails,input_present,theta,dnedges,dn,directed,bipartite,nterms,funnames,sonames); /* Probability mass under present parameters */
  log_ratio = log_ratio + (log_proposal - log_present);
  accept = MH_Decision(log_ratio);
  if (accept == 1) /* Proposal accepted: set ergm->theta and ls->theta to proposal */
    {
    Set_DD_DD(ls->d,ls->number+1,ls->theta,ls_theta);
    }
  if (print >= 1)
    {
    Rprintf("\nSample block parameters:");
    Rprintf("\n- M-H acceptance probability: %8.4f",Min(e(log_ratio),1.0));  
    Rprintf("\n- decision: %i",accept);
    }
  free(theta);
  free(present);
  for (i = 0; i < ls->d; i++)
    {
    free(cf[i]);
    free(ls_theta[i]);
    }
  free(cf);
  free(ls_theta);
  return accept;
}

int Sample_Parameters_Independence(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges, int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames, char **sonames, 
                        double *input_proposal, double *input_present, int print, double *scale_factor)
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/
{
  int accept, i;
  double **cf, *present, *ergm_theta, log_present, log_proposal, log_ratio, **ls_theta, *proposal, *theta_present, *theta_proposal;
  /* Proposal:
  note 1: all ls->theta such that ls->size >= ls->minimum_size and all ergm->theta are updated by random walk Metropolis-Hastings algorithm
  note 2: ratio of proposal pdfs cancels under random walk Metropolis-Hastings algorithm */
  log_ratio = 0.0;
  /* Propose ergm->theta: */
  if (ergm->d1 > 0)
    {
    cf = Scale(ergm->d1,ergm->d1,prior->cf1,scale_factor[0]); /* Rescale Cholesky factor of Gaussian prior */        
    ergm_theta = Sample_MVN(ergm->d1,ergm->theta,cf); /* Random walk Metropolis-Hastings */
    log_proposal = MVN_PDF(ergm->d1,ergm_theta,prior->mean1,prior->precision1); /* Prior pdf: proposal */
    log_present = MVN_PDF(ergm->d1,ergm->theta,prior->mean1,prior->precision1); /* Prior pdf: present */
    log_ratio = log_ratio + (log_proposal - log_present);
    for (i = 0; i < ergm->d1; i++)
      { 
      free(cf[i]);
      }
    free(cf);
    }
  /* Propose ls->theta: */
  ls_theta = (double**) calloc(ls->d,sizeof(double*)); /* Remark: since memory is allocated to ls_theta by calloc, between-block parameters are 0 */
  if (ls_theta == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Independence, ls_theta\n\n"); error("Error: out of memory"); }
  for (i = 0; i < ls->d; i++)
    {
    ls_theta[i] = (double*) calloc(ls->number+1,sizeof(double));
    if (ls_theta[i] == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Independence, ls_theta[%i]\n\n",i); error("Error: out of memory"); }
    }
  present = (double*) calloc(ls->d,sizeof(double));  
  if (present == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Independence, present\n\n"); error("Error: out of memory"); }
  cf = Scale(ls->d,ls->d,prior->cf2,scale_factor[1]); /* Rescale Cholesky factor of Gaussian prior */ 
  for (i = 0; i < ls->number; i++) 
    {
    Get_Column(ls->d,present,ls->theta,i); /* Set mean to ls->theta[][i] */
    if (ls->size[i] < ls->minimum_size) Set_Column(ls->d,ls_theta,i,present); /* Set proposal = present */ 
    else 
      {
      /* Generate candidate: */
      proposal = Sample_MVN(ls->d,present,cf); /* Random walk Metropolis-Hastings algorithm */
      Set_Column(ls->d,ls_theta,i,proposal); /* Set ls_theta[][i] to proposal */
      /* Add ratio of prior pdf: */
      log_proposal = MVN_PDF(ls->d,proposal,prior->mean2,prior->precision2); /* Prior pdf of proposal */
      log_present = MVN_PDF(ls->d,present,prior->mean2,prior->precision2); /* Prior pdf of present */
      log_ratio = log_ratio + (log_proposal - log_present);
      free(proposal);
      }
    }
  for (i = 0; i < ls->d; i++) /* Set between-block parameters */
    {
    ls_theta[i][ls->number] = ls->theta[i][ls->number];
    }
  /* Decide: */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls_theta,input_proposal); /* Set input given ls_theta */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input given ls->theta */
  theta_proposal = Get_Parameter(ergm->d,ergm->structural,ergm_theta); /* Set parameter; note: if ergm_d1 == 0, ergm_theta is not used */
  theta_present = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter; note: if ergm->d1 == 0, ergm_theta is not used */  
  log_proposal = PMF_Independence(ls,ergm,heads,tails,input_proposal,theta_proposal,dnedges,dn,directed,bipartite,nterms,funnames,sonames); /* Probability mass under proposed parameter */
  log_present = PMF_Independence(ls,ergm,heads,tails,input_present,theta_present,dnedges,dn,directed,bipartite,nterms,funnames,sonames); /* Probability mass under present parameters */
  log_ratio = log_ratio + (log_proposal - log_present);
  accept = MH_Decision(log_ratio);
  if (accept == 1) /* Proposal accepted: set ergm->theta and ls->theta to proposal */
    {
    if (ergm->d1 > 0) Set_D_D(ergm->d1,ergm->theta,ergm_theta);
    Set_DD_DD(ls->d,ls->number+1,ls->theta,ls_theta);
    }
  if (print >= 1)
    {
    Rprintf("\nSample parameters:");
    Rprintf("\n- log ratio: %8.4f",log_ratio);  
    Rprintf("\n- decision: %i",accept);
    }
  if (ergm->d1 > 0) free(ergm_theta);
  free(theta_present);
  free(theta_proposal);
  free(present);
  for (i = 0; i < ls->d; i++)
    {
    free(cf[i]);
    free(ls_theta[i]);
    }
  free(cf);
  free(ls_theta);
  return accept;
}

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
                        int *newnetworkheads, int *newnetworktails, double *scale_factor, int update_node, double *temperature, int *status)
/*
input: ergm structure, latent structure, prior
output: indicators
*/
{
  int burn_in, maxpossibleedges_block, *block, **edge_list_block, number_edges_block, *heads_block, *tails_block, size, number_networks, accept, mdnedges_block, *mheads_block, *mtails_block, auxiliary, i, k, large, *ls_indicator, *ls_size, n_input, present_block, proposal_block, proposal_distribution, proposal_n_edges, *proposal_heads, *proposal_tails, sample_size;
  double a_i, entropy, *input_proposal, log_denominator, log_numerator, log_present, log_proposal, log_ratio, *input_present_block, *input_proposal_block, *p, *q_i, present_a, proposal_a, present_energy, proposal_energy, sum, t, *theta, *statistic;
  n_input = Number_Input(ergm->terms,input_present);
  input_proposal = (double*) calloc(n_input,sizeof(double));
  if (input_proposal == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, input_proposal\n\n"); error("Error: out of memory"); }
  ls_indicator = (int*) calloc(ls->n,sizeof(int));
  if (ls_indicator == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, ls_indicator\n\n"); error("Error: out of memory"); }
  ls_size = (int*) calloc(ls->number,sizeof(int));
  if (ls_size == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, ls_size\n\n"); error("Error: out of memory"); }
  for (i = 0; i < ls->n; i++) 
    {
    ls_indicator[i] = ls->indicator[i];
    }
  for (k = 0; k < ls->number; k++)
    {
    ls_size[k] = ls->size[k];
    }
  /* Construct proposal distribution:
  - in general, proposal distribution is given by ls->p
  - in special cases, proposal distribution is approximated by full conditional distribution by using mean-field methods:
  there are two mean-field methods, one fast and one slow; both give rise to exact results when ls->size[k] == 2 and work well as long as ls->size[k] is not too large, but when ls->size[k] > 5, the slow method is much more accurate than the slow method */
  if (model == 0) proposal_distribution = 0; /* 666 */
  else proposal_distribution = 1;
  if (proposal_distribution == 0) 
    {
    q_i = (double*) calloc(ls->number,sizeof(double));
    if (q_i == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, q_i\n\n"); error("Error: out of memory"); }
    for (k = 0; k < ls->number; k++) 
      {
      q_i[k] = 1.0 / ls->number;
      }
    }
  else q_i = Candidate_Generating_Distribution_Indicators_Dependence(update_node,model,ls,ergm,heads,tails,input_present,dnedges,dn,directed,bipartite,nterms,funnames,sonames); /* Special case:  proposal distribution: exact or approximate full conditional distribution */
  /* 666 
  proposal_distribution = 0; 
  q_i = (double*) calloc(ls->number,sizeof(double));
  if (q_i == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, q_i\n\n"); error("Error: out of memory"); }
  for (k = 0; k < ls->number; k++) 
    {
    q_i[k] = 1.0 / ls->number;
    }
  666 */
  entropy = S(ls->number,q_i); /* Entropy of ls->p as indicator of how much the nodes are spread out across blocks */
  entropy = entropy / ln(ls->number);
  /*
  Rprintf("\nEntropy of ls->p: %8.4f",entropy);
  */
  /* Temperature: note that the entropy of the full conditional distribution may be strongly peaked, and high temperatures are required to unfreeze the algorithm */
  if ((entropy >= epsilon) && (entropy <= maximum)) t = 1.0 / entropy;
  else t = 1.0;
  if (t < temperature[0]) t = temperature[0]; 
  else if (t > temperature[1]) t = temperature[1];
  if (t != 1.0) /* Melt down the full conditional distribution, since the dependence of the indicators implies that the full conditional distribution may have low entropy */
    {
    p = (double*) calloc(ls->number,sizeof(double));
    if (p == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, p\n\n"); error("Error: out of memory"); }
    sum = 0.0;
    for (k = 0; k < ls->number; k++)
      {
      p[k] = q_i[k]; /* Store q_i[k], since it may be needed to compute the log acceptance ratio */
      q_i[k] = ln(q_i[k]) / t;
      sum = sum + e(q_i[k]);  
      }
    a_i = ln(sum);
    for (k = 0; k < ls->number; k++)
      {
      q_i[k] = e(q_i[k] - a_i);
      }
    entropy = S(ls->number,q_i); /* Entropy of ls->p as indicator of how much the nodes are spread out across blocks */
    entropy = entropy / ln(ls->number);
    }
  /* Generate proposal: */
  present_block = ls->indicator[update_node];
  proposal_block = Sample_Discrete(q_i); 
  ls_indicator[update_node] = proposal_block;
  ls_size[present_block] = ls_size[present_block] - 1;
  ls_size[proposal_block] = ls_size[proposal_block] + 1;
  /* Compute log acceptance ratio : */    
  log_ratio = 0.0;
  if ((ls->size[present_block] < 5) && (ls->size[proposal_block] < 5)) large = 0; /* The "smaller than" implies that both present and proposed sizes of blocks will be at most 5, which is sufficiently small */
  else large = 1;
  if (proposal_block == present_block) /* Special case: proposed == present block: log acceptance ratio vanishes */
    {
    auxiliary = 0; /* Metropolis-Hastings */ 
    log_ratio = 0.0;
    }
  else if (large == 0) /* Small blocks */
    {
    /* General remarks:

    Log likelihood ratio + log prior ratio
    --------------------------------------
    ln(p[proposal_block]) - ln(p[present_block]),
    where 
    - p[k] = ls->p[k] * e(energy_k - a_k),  
    - e(energy_k - a_k) is the PMF of the observed graph given that the node is member of block k,
    - energy_k is the energy of the observed graph,
    - a_k is the log partition function of the PMF of the observed graph and is given by a_k = sum_block a_within_block + a_between;
    if the proposed and present block are small before as well as after the proposed move, 
    a_within_proposed_block and a_within_present_block before and after the proposed move can be computed by complete enumeration,
    and because all other within-block partition functions cancel in the log likelihood ratio,
    the log acceptance ratio can computed exactly;
    note that the between-block partition functions can be computed as well

    Log likelihood ratio + log prior ratio + log proposal ratio
    -----------------------------------------------------------
    ln(p[proposal_block]) - ln(p[present_block]) + ln(q_i[present_block]) - ln(q_i[proposal_ratio])
    where 
    - p[k] is defined as above and q_i[k] is the propobability that block k is proposed,
    - q_i[k] may be given by
      * q_i[k] = p[k]: the log acceptance ratio vanishes;
      * q_i[k] = function(p[k], temperature t): the log acceptance ratio does not vanish, but can readily and exactly be computed, since p[k] has already been computed
      * q_i[k] = ls->p[k]: the log acceptance ratio does not vanish, but can exactly be computed

    */
    auxiliary = 0; /* Metropolis-Hastings */
    if (proposal_distribution == 1) /* Special case: proposed != present block, both are small before as well as after the proposed move, proposal distribution is given by full conditional distribution, computed either exactly by complete enumeration or approximately by mean-field methods */
      {
      if (t == 1.0) log_ratio = 0.0; /* Proposal distribution: full conditional distribution */
      else /* Proposal distribution: full conditional distribution at temperature t */
        {
        log_ratio = 0.0;
        log_ratio = log_ratio + ln(q_i[present_block]) - ln(q_i[proposal_block]); /* Log proposal ratio */
        /*
        Rprintf("\n- log ratio (log proposal ratio) = %8.4f",log_ratio);  
        */
        log_ratio = log_ratio + ln(ls->p[proposal_block]) - ln(ls->p[present_block]); /* Log prior ratio */
        /*
        Rprintf("\n- log ratio (log prior ratio) = %8.4f",log_ratio);  
        */
        log_ratio = log_ratio + ln(p[proposal_block]) - ln(p[present_block]); /* Log likelihood ratio */
        /*
        Rprintf("\n- log ratio (log likelihood ratio) = %8.4f",log_ratio);  
        */
        }
      }
    else /* Special case: proposed != present block, both are small before as well as after the proposed move, proposal distribution is given by ls->p */
      {
      log_ratio = 0.0;
      log_ratio = log_ratio + ln(q_i[present_block]) - ln(q_i[proposal_block]); /* Log proposal ratio */
      /*
      Rprintf("\n- log ratio (log proposal ratio) = %8.4f",log_ratio);  
      */
      log_ratio = log_ratio + ln(ls->p[proposal_block]) - ln(ls->p[present_block]); /* Log prior ratio */
      /*
      Rprintf("\n- log ratio (log prior ratio) = %8.4f",log_ratio);  
      */
      for (i = 0; i < n_input; i++) 
        { 
        input_proposal[i] = input_present[i];
        }
      Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls_indicator,ls->theta,input_proposal); /* Set input given ls_indicator */
      Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input given ls->indicator */
      theta = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter; note: if ergm->d1 == 0, ergm->theta is not used */
      statistic = (double*) calloc(ergm->d,sizeof(double));
      if (statistic == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, statistic\n\n"); error("Error: out of memory"); }
      ls->indicator[update_node] = proposal_block; /* Set ls->indicator to proposed block */
      ls->size[present_block] = ls->size[present_block] - 1;
      ls->size[proposal_block] = ls->size[proposal_block] + 1;
      proposal_energy = Minus_Energy(ergm->d,input_proposal,theta,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,statistic);
      proposal_a = 0.0;
      proposal_a = proposal_a + Within_Block_Partition_Function(model,ls,present_block,ergm,heads,tails,input_proposal,dn,directed,nterms,funnames,sonames);
      proposal_a = proposal_a + Within_Block_Partition_Function(model,ls,proposal_block,ergm,heads,tails,input_proposal,dn,directed,nterms,funnames,sonames);
      proposal_a = proposal_a + Between_Block_Partition_Function(ls,ergm,input_proposal,theta,dn,directed,bipartite,nterms,funnames,sonames);
      log_proposal = proposal_energy - proposal_a;
      ls->indicator[update_node] = present_block; /* Reset ls->indicator to present block */
      ls->size[proposal_block] = ls->size[proposal_block] - 1;
      ls->size[present_block] = ls->size[present_block] + 1;
      present_energy = Minus_Energy(ergm->d,input_present,theta,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,statistic);
      present_a = 0.0;
      present_a = present_a + Within_Block_Partition_Function(model,ls,present_block,ergm,heads,tails,input_present,dn,directed,nterms,funnames,sonames);
      present_a = present_a + Within_Block_Partition_Function(model,ls,proposal_block,ergm,heads,tails,input_present,dn,directed,nterms,funnames,sonames);
      present_a = present_a + Between_Block_Partition_Function(ls,ergm,input_present,theta,dn,directed,bipartite,nterms,funnames,sonames);
      log_present = present_energy - present_a;
      log_ratio = log_ratio + (log_proposal - log_present);
      free(theta);
      free(statistic);
      /*
      Rprintf("\n- node %i: block %i > %i (block size %i > %i): log_ratio: %8.4f",update_node+1,ls->indicator[update_node]+1,ls_indicator[update_node]+1,ls->size[present_block],ls_size[proposal_block],log_ratio);
      */
      }
    }
  else if (model > 0) /* Computation of PMF by complete enumeration infeasible: without covariates, generate local sample */
    {
    auxiliary = 1; /* Auxiliary-variable Metropolis-Hastings */
    log_ratio = 0.0;
    log_ratio = log_ratio + ln(q_i[present_block]) - ln(q_i[proposal_block]); /* Log proposal ratio */
    /*
    Rprintf("\n- log ratio (log proposal ratio) = %8.4f",log_ratio);  
    */
    log_ratio = log_ratio + ln(ls->p[proposal_block]) - ln(ls->p[present_block]); /* Log prior ratio */
    /*
    Rprintf("\n- log ratio (log prior ratio) = %8.4f",log_ratio);  
    */
    for (i = 0; i < n_input; i++)
      {
      input_proposal[i] = input_present[i];
      }
    Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls_indicator,ls->theta,input_proposal); /* Set input given ls_indicator */
    Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input given ls->indicator */
    theta = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter; note: if ergm->d1 == 0, ergm->theta is not used */
    statistic = (double*) calloc(ergm->d,sizeof(double));
    if (statistic == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, statistic\n\n"); error("Error: out of memory"); }
    size = ls->size[present_block] + ls->size[proposal_block];
    if (*directed == 0) maxpossibleedges_block = size * (size - 1) / 2;
    else maxpossibleedges_block = size * (size - 1);
    /* Extract subgraph corresponding to nodes which are members of the block from observed graph: */
    block = (int*) calloc(3,sizeof(int));
    if (block == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, block\n\n"); error("Error: out of memory"); }
    if (present_block == proposal_block) 
      {
      block[0] = 1; /* Number of blocks included; the labels of included blocks are stored in block[1], ..., block[number_blocks] */
      block[1] = present_block;
      }
    else
      {
      block[0] = 2; /* Number of blocks included; the labels of included blocks are stored in block[1], ..., block[number_blocks] */
      block[1] = present_block;
      block[2] = proposal_block;
      }
    edge_list_block = Edge_List_Blocks(ls,block,dnedges,heads,tails);
    number_edges_block = edge_list_block[0][0];
    heads_block = &edge_list_block[1][0];
    tails_block = &edge_list_block[2][0];
    /* 
    Rprintf("\nedge list of blocks %i and %i with %i and %i nodes, respectively, and %i edges",present_block,proposal_block,ls->size[present_block],ls->size[proposal_block],edge_list_block[0][0]);
    Print_I(number_edges_block,heads_block);
    Print_I(number_edges_block,tails_block);
    */
    input_present_block = Extract_Input_Blocks(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,input_present,block,ls->theta);
    /*
    Rprintf("\nUpdating blocks %i and %i with %i and %i nodes, respectively: present.",present_block,proposal_block,ls->size[present_block],ls->size[proposal_block]);
    n_input = Number_Input(ergm->terms,input_present_block);
    Print_D(n_input,input_present_block);
    */
    input_proposal_block = Extract_Input_Blocks(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls_indicator,input_proposal,block,ls->theta);
    /* 
    Rprintf("\nUpdating blocks %i and %i with %i and %i nodes, respectively: proposal.",present_block,proposal_block,ls->size[present_block],ls->size[proposal_block]);
    n_input = Number_Input(ergm->terms,input_proposal_block);
    Print_D(n_input,input_proposal_block);
    */
    free(block);
    mdnedges_block = 0;
    mheads_block = NULL;
    mtails_block = NULL;
    number_networks = 1;
    sample_size = 1;
    if (20 * maxpossibleedges_block < 10000) burn_in = 10000;
    else burn_in = 20 * maxpossibleedges_block;
    MCMC_wrapper(&number_networks,&number_edges_block,tails_block,heads_block,  /* Sample one subgraph from posterior predictive distribution given input and theta */
                    &size,directed,bipartite, /* Number of nodes of the subgraph */ 
                    nterms,funnames,
                    sonames, 
                    MHproposaltype,MHproposalpackage,
                    input_proposal_block,theta,&sample_size,
                    sample,&burn_in,interval,  
                    newnetworkheads, 
                    newnetworktails, 
                    verbose, 
                    attribs,maxout,maxin,minout,
                    minin,condAllDegExact,attriblength, 
                    maxedges,
                    status);
    if (print >= 0)
      {
      if (*status == 1) Rprintf("\nWARNING: Sample_Indicators_Dependence: number of edges %i is outside of (1, %i).",newnetworkheads[0],*maxedges-1);
      else if (*status == 2) Rprintf("\nWARNING: M-H proposal failed.");
      }
    /*
    Rprintf("\n- maximum number of edges: %i",maxpossibleedges_block);
    Rprintf("\n- number of edges of observed graph: %i",number_edges_block);
    Rprintf("\n- number of edges of auxiliary graph: %i",newnetworkheads[0]);
    */
    proposal_n_edges = newnetworkheads[0]; /* Number of simulated edges */
    proposal_heads = (int*) calloc(proposal_n_edges,sizeof(int)); /* Proposed heads for auxiliary variable */
    if (proposal_heads == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, proposal_heads\n\n"); error("Error: out of memory"); }
    proposal_tails = (int*) calloc(proposal_n_edges,sizeof(int)); /* Proposed tails for auxiliary variable */
    if (proposal_tails == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, proposal_tails\n\n"); error("Error: out of memory"); }
    for (i = 0; i < proposal_n_edges; i++)  
      {
      proposal_heads[i] = newnetworkheads[i+1]; /* Note: while heads corresponds to the list of observed heads, newnetworkheads contains the number of   simulated edges as well as the list of simulated heads: to use auxiliary->heads here, one must not store the number of simulated edges */
      proposal_tails[i] = newnetworktails[i+1]; /* Note: while tails corresponds to the list of observed tails, newnetworktails contains the number of   simulated edges as well as the list of simulated tails: to use auxiliary->tails here, one must not store the number of simulated edges */
      }
    /*
    Rprintf("\n- heads and tails of auxiliary graph:\n");
    for (i = 0; i < proposal_n_edges; i++)  
      {
      Rprintf(" %i",proposal_heads[i]); 
      }
    Rprintf("\n");
    for (i = 0; i < proposal_n_edges; i++)  
      {
      Rprintf(" %i",proposal_tails[i]); 
      }
    Rprintf("\n");
    */
    /* Ratio of proposal pmfs of auxiliary graph under proposal / present */
    log_numerator = Minus_Energy(ergm->d,input_present_block,theta,
                    proposal_heads,proposal_tails,&proposal_n_edges,&size,directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n");
    Rprintf("\n- function 1: log numerator = %8.4f",log_numerator);  
    Rprintf("\n- statistic of auxiliary graph:");
    Print_D(ergm->d,statistic);
    */
    log_denominator = Minus_Energy(ergm->d,input_proposal_block,theta,
                      proposal_heads,proposal_tails,&proposal_n_edges,&size,directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n- function 2: log denominator = %8.4f",log_denominator);  
    Rprintf("\n- statistic of auxiliary graph:");
    Print_D(ergm->d,statistic);
    */
    log_ratio = log_ratio + (log_numerator - log_denominator);
    /*
    Rprintf("\n- log ratio (auxiliary graph) = %8.4f",log_ratio);  
    */
    /* Ratio of mass of observed graph under proposal / present */
    log_present = Minus_Energy(ergm->d,input_present_block,theta,heads_block,tails_block,&number_edges_block,&size,directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n");
    Rprintf("\n- function 3: log present = %8.4f",log_present);  
    Rprintf("\n- statistic of observed graph:");
    Print_D(ergm->d,statistic);
    */
    log_proposal = Minus_Energy(ergm->d,input_proposal_block,theta,heads_block,tails_block,&number_edges_block,
                   &size,directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n- function 4: log proposal = %8.4f",log_proposal);  
    Rprintf("\n- statistic of observed graph:");
    Print_D(ergm->d,statistic);
    */
    log_ratio = log_ratio + (log_proposal - log_present);
    /*
    Rprintf("\n- log ratio (observed graph) = %8.4f",log_ratio);  
    */
    for (i = 0; i < 3; i++) 
      {
      free(edge_list_block[i]);
      }
    free(edge_list_block);
    free(input_present_block);
    free(input_proposal_block);
    }
  else /* Large blocks: introduce auxiliary variables */
    {
    auxiliary = 1; /* Auxiliary-variable Metropolis-Hastings */
    log_ratio = 0.0;
    log_ratio = log_ratio + ln(q_i[present_block]) - ln(q_i[proposal_block]); /* Log proposal ratio */
    /*
    Rprintf("\n- log ratio (log proposal ratio) = %8.4f",log_ratio);  
    */
    log_ratio = log_ratio + ln(ls->p[proposal_block]) - ln(ls->p[present_block]); /* Log prior ratio */
    /*
    Rprintf("\n- log ratio (log prior ratio) = %8.4f",log_ratio);  
    */
    for (i = 0; i < n_input; i++) 
      { 
      input_proposal[i] = input_present[i];
      }
    Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls_indicator,ls->theta,input_proposal); /* Set input given ls_indicator */
    Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input given ls->indicator */
    theta = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter; note: if ergm->d1 == 0, ergm->theta is not used */
    statistic = (double*) calloc(ergm->d,sizeof(double));
    if (statistic == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, statistic\n\n"); error("Error: out of memory"); }
    sample_size = 1; /* One sample point is all that is required */
    number_networks = 1;
    /*
    Rprintf("\n\nbefore: newnetworkheads[0]=%i",newnetworkheads[0]);
    */
    MCMC_wrapper(&number_networks,dnedges,tails,heads,  /* Sample one graph from posterior predictive distribution given input and theta */
                    dn,directed,bipartite, 
                    nterms,funnames,
                    sonames, 
                    MHproposaltype,MHproposalpackage,
                    input_proposal,theta,&sample_size,
                    sample,burnin,interval,  
                    newnetworkheads, 
                    newnetworktails, 
                    verbose, 
                    attribs,maxout,maxin,minout,
                    minin,condAllDegExact,attriblength, 
                    maxedges,
                    status);
    if (print >= 0)
      {
      if (*status == 1) Rprintf("\nWARNING: Sample_Indicators_Dependence: number of edges %i is outside of (1, %i).",newnetworkheads[0],*maxedges-1);
      else if (*status == 2) Rprintf("\nWARNING: M-H proposal failed.");
      }
    /*
    Rprintf("\n\nafter: newnetworkheads[0]=%i",newnetworkheads[0]);
    */
    proposal_n_edges = newnetworkheads[0]; /* Number of simulated edges */
    proposal_heads = (int*) calloc(proposal_n_edges,sizeof(int)); /* Proposed heads for auxiliary variable */
    if (proposal_heads == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, proposal_heads\n\n"); error("Error: out of memory"); }
    proposal_tails = (int*) calloc(proposal_n_edges,sizeof(int)); /* Proposed tails for auxiliary variable */
    if (proposal_tails == NULL) { Rprintf("\n\ncalloc failed: Sample_Indicators_Dependence, proposal_tails\n\n"); error("Error: out of memory"); }
    for (i = 0; i < proposal_n_edges; i++)  
      {
      proposal_heads[i] = newnetworkheads[i+1]; /* Note: while heads corresponds to the list of observed heads, newnetworkheads contains the number of   simulated edges as well as the list of simulated heads: to use auxiliary->heads here, one must not store the number of simulated edges */
      proposal_tails[i] = newnetworktails[i+1]; /* Note: while tails corresponds to the list of observed tails, newnetworktails contains the number of   simulated edges as well as the list of simulated tails: to use auxiliary->tails here, one must not store the number of simulated edges */
      }
    /* Ratio of proposal pmfs of auxiliary graph under proposal / present */
    log_numerator = Minus_Energy(ergm->d,input_present,theta,
    proposal_heads,proposal_tails,&proposal_n_edges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n- function 1: log numerator = %8.4f",log_numerator);  
    */
    log_denominator = Minus_Energy(ergm->d,input_proposal,theta,
    proposal_heads,proposal_tails,&proposal_n_edges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n- function 2: log denominator = %8.4f",log_denominator);  
    */
    log_ratio = log_ratio + (log_numerator - log_denominator);
    /*
    Rprintf("\n- log ratio (auxiliary graph) = %8.4f",log_ratio);  
    */
    /* Ratio of mass of observed graph under proposal / present */
    log_present = Minus_Energy(ergm->d,input_present,theta,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n- function 3: log present = %8.4f",log_present);  
    */
    log_proposal = Minus_Energy(ergm->d,input_proposal,theta,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n- function 4: log proposal = %8.4f",log_proposal);  
    */
    log_ratio = log_ratio + (log_proposal - log_present);
    /*
    Rprintf("\n- log ratio (observed graph) = %8.4f",log_ratio);  
    */
    /*
    Rprintf("\n- log preference of present indicator wrt data: %8.4f",log_numerator-log_present);  
    Rprintf("\n- log preference of propposed indicator wrt data: %8.4f",log_denominator-log_proposal);  
    Rprintf("\n- log exchange ratio: %8.4f",log_numerator-log_denominator+log_proposal-log_present);  
    */
    }
  accept = MH_Decision(log_ratio);
  /* Update ls->indicator and ls-size */
  if ((*status == 0) && (accept == 1)) /* Proposal accepted */
    {
    ls->indicator[update_node] = proposal_block;
    ls->size[present_block] = ls->size[present_block] - 1;
    ls->size[proposal_block] = ls->size[proposal_block] + 1;
    }
  /* Console output: */
  if (print >= 1)
    {
    Rprintf("\nSample indicator of node %i:",update_node+1);
    entropy = S(ls->number,q_i);
    entropy = entropy / (ln(ls->number));
    Rprintf("\n- entropy of proposal distribution: %2.2f%%",entropy*100);
    Rprintf("\n- proposal distribution:");
    for (k = 0; k < ls->number; k++)
      {
      Rprintf("%7.4f",q_i[k]);
      }
    Rprintf("\n- proposal: %i > %i",present_block+1,proposal_block+1);      
    if (auxiliary == 0) Rprintf("\n- M-H acceptance probability: %8.4f",Min(e(log_ratio),1.0));  
    else Rprintf("\n- auxiliary-variable M-H acceptance probability: %8.4f",Min(e(log_ratio),1.0));  
    Rprintf("\n- decision: %i",accept);
    }
  /* Free memory: */
  free(input_proposal);
  free(ls_indicator);
  free(ls_size);
  free(q_i); /* 666 */
  if (t != 1.0) free(p);
  if (auxiliary > 0)
    {
    free(proposal_heads);
    free(proposal_tails);
    free(statistic);
    free(theta);
    }
  return accept;
}

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
                        int *newnetworkheads, int *newnetworktails, double *scale_factor, int *status)
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/
{
  int number_networks, accept, i, n_input, proposal_n_edges, *proposal_heads, *proposal_tails, sample_size;
  double **cf, *ergm_theta, *input_proposal, log_denominator, log_numerator, log_present, log_proposal, log_ratio, *theta_present, *theta_proposal, *statistic;
  n_input = Number_Input(ergm->terms,input_present);
  input_proposal = (double*) calloc(n_input,sizeof(double));
  if (input_proposal == NULL) { Rprintf("\n\ncalloc failed: Sample_Ergm_Theta_Dependence, input_proposal\n\n"); error("Error: out of memory"); }
  log_ratio = 0.0;
  /* Propose ergm->theta: */
  cf = Scale(ergm->d1,ergm->d1,prior->cf1,scale_factor[0]); /* Rescale Cholesky factor of Gaussian prior */        
  ergm_theta = Sample_MVN(ergm->d1,ergm->theta,cf); /* Random walk Metropolis-Hastings */
  log_proposal = MVN_PDF(ergm->d1,ergm_theta,prior->mean1,prior->precision1); /* Prior pdf: proposal */
  log_present = MVN_PDF(ergm->d1,ergm->theta,prior->mean1,prior->precision1); /* Prior pdf: present */
  log_ratio = log_ratio + (log_proposal - log_present);
  /*
  Rprintf("\n- log ratio (parameters) = %8.4f",log_ratio);  
  */
  /* Propose auxiliary variable: */
  for (i = 0; i < n_input; i++) 
    { 
    input_proposal[i] = input_present[i];
    }
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_proposal); /* Set input given ls->theta */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input given ls->theta */
  theta_proposal = Get_Parameter(ergm->d,ergm->structural,ergm_theta); /* Set parameter; note: if ergm_d1 == 0, ergm_theta is not used */
  theta_present = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter; note: if ergm->d1 == 0, ergm_theta is not used */
  sample_size = 1; /* One sample point is all that is required */
  number_networks = 1;
  /*
  Rprintf("\nburnin=%i interval=%i maxout=%i maxin=%i minout=%i minin=%i",*burnin,*interval,*maxout,*maxin,*minout,*minin);
  Print_D(n_input,input_present);
  Print_D(n_input,input_proposal);
  Print_DD(ls->d,ls->number,ls->theta);
  Print_I(ls->n,ls->indicator);
  Print_D(ergm->d1,ergm->theta);
  Print_D(ergm->d1,ergm_theta);
  Print_D(ergm->d,theta_present);
  Print_D(ergm->d,theta_proposal);
  Print_I(*dnedges,heads);
  Print_I(*dnedges,tails);
  */
  MCMC_wrapper(&number_networks,dnedges,tails,heads,  /* Sample one graph from posterior predictive distribution given input and theta */
                  dn,directed,bipartite, 
                  nterms,funnames,
                  sonames, 
                  MHproposaltype,MHproposalpackage,
                  input_proposal,theta_proposal,&sample_size,
                  sample,burnin,interval,  
                  newnetworkheads, 
                  newnetworktails, 
                  verbose, 
                  attribs,maxout,maxin,minout,
                  minin,condAllDegExact,attriblength, 
                  maxedges,
                  status);
  if (print >= 0)
    {
    if (*status == 1) Rprintf("\nWARNING: Sample_Ergm_Theta_Dependence: number of edges %i is outside of (1, %i).",newnetworkheads[0],*maxedges-1);
    else if (*status == 2) Rprintf("\nWARNING: M-H proposal failed.");
    }
  proposal_n_edges = newnetworkheads[0]; /* Number of simulated edges */
  proposal_heads = (int*) calloc(proposal_n_edges,sizeof(int)); /* Proposed heads for auxiliary variable */
  if (proposal_heads == NULL) { Rprintf("\n\ncalloc failed: Sample_Ergm_Theta_Dependence, proposal_heads\n\n"); error("Error: out of memory"); }
  proposal_tails = (int*) calloc(proposal_n_edges,sizeof(int)); /* Proposed tails for auxiliary variable */
  if (proposal_tails == NULL) { Rprintf("\n\ncalloc failed: Sample_Ergm_Theta_Dependence, proposal_tails\n\n"); error("Error: out of memory"); }
  for (i = 0; i < proposal_n_edges; i++)  
    {
    proposal_heads[i] = newnetworkheads[i+1]; /* Note: while heads corresponds to the list of observed heads, newnetworkheads contains the number of   simulated edges as well as the list of simulated heads: to use auxiliary->heads here, one must not store the number of simulated edges */
    proposal_tails[i] = newnetworktails[i+1]; /* Note: while tails corresponds to the list of observed tails, newnetworktails contains the number of   simulated edges as well as the list of simulated tails: to use auxiliary->tails here, one must not store the number of simulated edges */
    }
  /* Ratio of proposal pmfs of auxiliary graph under proposal / present */
  statistic = (double*) calloc(ergm->d,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: Sample_Ergm_Theta_Dependence, statistic\n\n"); error("Error: out of memory"); }
  log_numerator = Minus_Energy(ergm->d,input_present,theta_present,
  proposal_heads,proposal_tails,&proposal_n_edges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
  /*
  Rprintf("\n- function 1: log numerator = %8.4f",log_numerator);  
  */
  log_denominator = Minus_Energy(ergm->d,input_proposal,theta_proposal,
  proposal_heads,proposal_tails,&proposal_n_edges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
  /*
  Rprintf("\n- function 2: log denominator = %8.4f",log_denominator);  
  */
  log_ratio = log_ratio + (log_numerator - log_denominator);
  /*
  Rprintf("\n- log ratio (auxiliary graph) = %8.4f",log_ratio);  
  */
  /* Ratio of mass of observed graph under proposal / present */
  log_present = Minus_Energy(ergm->d,input_present,theta_present,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
  /*
  Rprintf("\n- function 3: log present = %8.4f",log_present);  
  */
  log_proposal = Minus_Energy(ergm->d,input_proposal,theta_proposal,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
  /*
  Rprintf("\n- function 4: log proposal = %8.4f",log_proposal);  
  */
  log_ratio = log_ratio + (log_proposal - log_present);
  /*
  Rprintf("\n- log ratio (observed graph) = %8.4f",log_ratio);  
  */
  accept = MH_Decision(log_ratio);
  if ((*status == 0) && (accept == 1)) /* Proposal accepted */
    {
    if (ergm->d1 > 0) Set_D_D(ergm->d1,ergm->theta,ergm_theta);
    Set_DD_DD(ls->d,ls->number+1,ls->theta,ls->theta);
    }
  if (print >= 1)
    {
    Rprintf("\nSample parameters:");
    Rprintf("\n- auxiliary-variable M-H acceptance probability: %8.4f",Min(e(log_ratio),1.0)); 
    Rprintf("\n- decision: %i",accept);
    }
  for (i = 0; i < ergm->d1; i++)
    { 
    free(cf[i]);
    }
  free(cf);
  free(ergm_theta);
  free(proposal_heads);
  free(proposal_tails);
  free(statistic);
  free(theta_present);
  free(theta_proposal);
  return accept;
}

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
                        int *newnetworkheads, int *newnetworktails, double *scale_factor, int update_block, int *status)
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/
{
  int burn_in, *block, number_networks, accept, auxiliary, **edge_list_block, i, k, *heads_block, *tails_block, mdnedges_block, *mheads_block, *mtails_block, number_edges_block, maxpossibleedges_block, n_input, proposal_n_edges, *proposal_heads, *proposal_tails, sample_size, update_size;
  double **cf, *present, *input_proposal, *input_present_block, *input_proposal_block, log_denominator, log_numerator, log_present, log_proposal, log_ratio, **ls_theta, *proposal, *theta_present, *theta_proposal, *statistic, present_a, proposal_a, present_energy, proposal_energy;
  n_input = Number_Input(ergm->terms,input_present);
  input_proposal = (double*) calloc(n_input,sizeof(double));
  if (input_proposal == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Dependence, input_proposal\n\n"); error("Error: out of memory"); }
  ls_theta = (double**) calloc(ls->d,sizeof(double*));
  if (ls_theta == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Dependence, ls_theta\n\n"); error("Error: out of memory"); }
  for (i = 0; i < ls->d; i++)
    {
    ls_theta[i] = (double*) calloc(ls->number+1,sizeof(double));
    if (ls_theta[i] == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Dependence, ls_theta[%i]\n\n",i); error("Error: out of memory"); }
    }
  for (i = 0; i < ls->d; i++)
    { 
    for (k = 0; k < ls->number + 1; k++)
      {
      ls_theta[i][k] = ls->theta[i][k];
      }
    }
  /* Propose ls->theta: */
  present = (double*) calloc(ls->d,sizeof(double));  
  if (present == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Dependence, present\n\n"); error("Error: out of memory"); }
  statistic = (double*) calloc(ergm->d,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Dependence, statistic\n\n"); error("Error: out of memory"); }
  update_size = ls->size[update_block];
  cf = Scale(ls->d,ls->d,prior->cf2,scale_factor[update_size-2]); /* Rescale Cholesky factor of Gaussian prior */ 
  Get_Column(ls->d,present,ls->theta,update_block); /* Set present to ls->theta[][i] */
  /* Generate candidate: */
  proposal = Sample_MVN(ls->d,present,cf); /* Random walk Metropolis-Hastings algorithm */
  Set_Column(ls->d,ls_theta,update_block,proposal); /* Set ls_theta[][i] to proposal */
  /* Log acceptance ratio: */
  log_ratio = 0.0;
  log_proposal = MVN_PDF(ls->d,proposal,prior->mean2,prior->precision2); 
  log_present = MVN_PDF(ls->d,present,prior->mean2,prior->precision2); 
  log_ratio = log_ratio + (log_proposal - log_present); /* Log prior ratio */
  /*
  Rprintf("\n- log ratio: %8.4f",log_ratio);  
  */
  for (i = 0; i < n_input; i++) 
    { 
    input_proposal[i] = input_present[i];
    }
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls_theta,input_proposal); /* Set input given ls_theta */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input given ls->theta */
  theta_proposal = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter; note: if ergm_d1 == 0, ergm_theta is not used */
  theta_present = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter; note: if ergm->d1 == 0, ergm_theta is not used */
  if (ls->size[update_block] < ls->threshold) /* Computation of PMF by complete enumeration feasible */
    {
    auxiliary = 0; /* Metropolis-Hastings */
    proposal_energy = Minus_Energy(ergm->d,input_proposal,theta_proposal,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,statistic);
    Set_Column(ls->d,ls->theta,update_block,proposal); /* Set ls->theta to proposal */
    proposal_a = Within_Block_Partition_Function(model,ls,update_block,ergm,heads,tails,input_proposal,dn,directed,nterms,funnames,sonames);
    Set_Column(ls->d,ls->theta,update_block,present); /* Reset ls->theta to present */
    log_proposal = proposal_energy - proposal_a;
    present_energy = Minus_Energy(ergm->d,input_present,theta_present,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,statistic);
    present_a = Within_Block_Partition_Function(model,ls,update_block,ergm,heads,tails,input_present,dn,directed,nterms,funnames,sonames);
    log_present = present_energy - present_a;
    log_ratio = log_ratio + (log_proposal - log_present);
    /* 
    Rprintf("\nlog_proposal %8.4f log_present %8.4f log_ratio %8.4f",log_proposal,log_present,log_ratio);
    */
    }
  else if (model > 0) /* Computation of PMF by complete enumeration infeasible: without covariates, generate local within-block sample */
    {
    auxiliary = 1; /* Auxiliary-variable Metropolis-Hastings */
    if (*directed == 0) maxpossibleedges_block = ls->size[update_block] * (ls->size[update_block] - 1) / 2;
    else maxpossibleedges_block = ls->size[update_block] * (ls->size[update_block] - 1);
    /* Extract subgraph corresponding to nodes which are members of the block from observed graph: */
    block = (int*) calloc(2,sizeof(int));
    if (block == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Dependence, block\n\n"); error("Error: out of memory"); }
    block[0] = 1; /* Number of blocks included; the labels of included blocks are stored in block[1], ..., block[number_blocks] */
    block[1] = update_block;
    edge_list_block = Edge_List_Blocks(ls,block,dnedges,heads,tails);
    number_edges_block = edge_list_block[0][0];
    heads_block = &edge_list_block[1][0];
    tails_block = &edge_list_block[2][0];
    /*
    Rprintf("\nedge list with %i nodes and %i edges",ls->n,*dnedges);
    Print_I(*dnedges,heads);
    Print_I(*dnedges,tails);
    Rprintf("\nedge list of block %i with %i nodes and %i edges",update_block,ls->size[update_block],edge_list_block[0][0]);
    Print_I(number_edges_block,heads_block);
    Print_I(number_edges_block,tails_block);
    */
    input_present_block = Extract_Input_Blocks(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,input_present,block,ls->theta);
    /*    
    Rprintf("\nupdating block %i with %i nodes.",update_block,ls->size[update_block]);
    Rprintf("\n- present ergm->theta:");
    Print_D(ergm->d,theta_present);
    Rprintf("\n- present ls->theta:");
    Print_D(ls->d,present);
    Rprintf("\n- present input 2:");
    n_input = Number_Input(ergm->terms,input_present_block);
    Print_D(n_input,input_present_block);
    */
    input_proposal_block = Extract_Input_Blocks(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,input_proposal,block,ls_theta);
    /* 
    Rprintf("\nUpdating block %i with %i nodes.",update_block,ls->size[update_block]);
    Rprintf("\n- proposed ergm->theta:");
    Print_D(ergm->d,theta_proposal);
    Rprintf("\n- proposed ls->theta:");
    Print_D(ls->d,proposal);
    Rprintf("\n- proposed input 2:");
    n_input = Number_Input(ergm->terms,input_proposal_block);
    Print_D(n_input,input_proposal_block);
    */
    free(block);
    mdnedges_block = 0;
    mheads_block = NULL;
    mtails_block = NULL;
    number_networks = 1;
    sample_size = 1;
    if (20 * maxpossibleedges_block < 10000) burn_in = 10000;
    else burn_in = 20 * maxpossibleedges_block;
    /*
    if (print >= 1) *verbose = 5; 1 short, >5 long
    */
    MCMC_wrapper(&number_networks,&number_edges_block,tails_block,heads_block,  /* Sample one graph from posterior predictive distribution given input and theta */
                    &ls->size[update_block],directed,bipartite, 
                    nterms,funnames,
                    sonames, 
                    MHproposaltype,MHproposalpackage,
                    input_proposal_block,theta_proposal,&sample_size,
                    sample,&burn_in,interval,  
                    newnetworkheads, 
                    newnetworktails, 
                    verbose, 
                    attribs,maxout,maxin,minout,
                    minin,condAllDegExact,attriblength, 
                    maxedges,
                    status);
    if (print >= 0)
      {
      if (*status == 1) Rprintf("\nWARNING: Sample_Ls_Theta_Dependence block %i of size %i: number of edges %i is outside of (1, %i).",update_block,ls->size[update_block]+1,newnetworkheads[0],*maxedges-1);
      else if (*status == 2) Rprintf("\nWARNING: M-H proposal failed.");
      }
    /*
    if (print >= 1) *verbose = 0;
    */
    proposal_n_edges = newnetworkheads[0]; /* Number of simulated edges */
    proposal_heads = (int*) calloc(proposal_n_edges,sizeof(int)); /* Proposed heads for auxiliary variable */
    if (proposal_heads == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Dependence, proposal_heads\n\n"); error("Error: out of memory"); }
    proposal_tails = (int*) calloc(proposal_n_edges,sizeof(int)); /* Proposed tails for auxiliary variable */
    if (proposal_tails == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Dependence, proposal_tails\n\n"); error("Error: out of memory"); }
    for (i = 0; i < proposal_n_edges; i++)  
      {
      proposal_heads[i] = newnetworkheads[i+1]; /* Note: while heads corresponds to the list of observed heads, newnetworkheads contains the number of   simulated edges as well as the list of simulated heads: to use auxiliary->heads here, one must not store the number of simulated edges */
      proposal_tails[i] = newnetworktails[i+1]; /* Note: while tails corresponds to the list of observed tails, newnetworktails contains the number of   simulated edges as well as the list of simulated tails: to use auxiliary->tails here, one must not store the number of simulated edges */
      }
    /*
    if (print >= 1)
      {
      Rprintf("\n- number of edges of observed graph: %i",number_edges_block);
      Rprintf("\n- number of edges of auxiliary graph: %i",proposal_n_edges);
      }
    Rprintf("\n- heads and tails of auxiliary graph:\n");
    for (i = 0; i < proposal_n_edges; i++)  
      {
      Rprintf(" %i",proposal_heads[i]); 
      }
    Rprintf("\n");
    for (i = 0; i < proposal_n_edges; i++)  
      {
      Rprintf(" %i",proposal_tails[i]); 
      }
    Rprintf("\n");
    */
    /* Ratio of proposal pmfs of auxiliary graph under proposal / present */
    log_numerator = Minus_Energy(ergm->d,input_present_block,theta_present,
                    proposal_heads,proposal_tails,&proposal_n_edges,&ls->size[update_block],directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n");
    Rprintf("\n- function 1: log numerator = %8.4f",log_numerator);  
    Rprintf("\n- statistic of auxiliary graph:");
    Print_D(ergm->d,statistic);
    */
    log_denominator = Minus_Energy(ergm->d,input_proposal_block,theta_proposal,
                      proposal_heads,proposal_tails,&proposal_n_edges,&ls->size[update_block],directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n- function 2: log denominator = %8.4f",log_denominator);  
    Rprintf("\n- statistic of auxiliary graph:");
    Print_D(ergm->d,statistic);
    */
    log_ratio = log_ratio + (log_numerator - log_denominator);
    /* 
    Rprintf("\n- log ratio (auxiliary graph) = %8.4f",log_ratio);  
    */
    /* Ratio of mass of observed graph under proposal / present */
    log_present = Minus_Energy(ergm->d,input_present_block,theta_present,heads_block,tails_block,&number_edges_block,
                  &ls->size[update_block],directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n");
    Rprintf("\n- function 3: log present = %8.4f",log_present);  
    Rprintf("\n- statistic of observed graph:");
    Print_D(ergm->d,statistic);
    */
    log_proposal = Minus_Energy(ergm->d,input_proposal_block,theta_proposal,heads_block,tails_block,&number_edges_block,
                   &ls->size[update_block],directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n- function 4: log proposal = %8.4f",log_proposal);  
    Rprintf("\n- statistic of observed graph:");
    Print_D(ergm->d,statistic);
    */
    log_ratio = log_ratio + (log_proposal - log_present);
    /*
    Rprintf("\n- log ratio (observed graph) = %8.4f",log_ratio);  
    */
    for (i = 0; i < 3; i++) 
      {
      free(edge_list_block[i]);
      }
    free(edge_list_block);
    free(input_present_block);
    free(input_proposal_block);
    free(proposal_heads);
    free(proposal_tails);
    }
  else /* Computation of PMF by complete enumeration infeasible: with covariates, generate global sample */
    {
    auxiliary = 1; /* Auxiliary-variable Metropolis-Hastings */
    sample_size = 1; /* One sample point is all that is required */
    number_networks = 1;
    MCMC_wrapper(&number_networks,dnedges,tails,heads,  /* Sample one graph from posterior predictive distribution given input and theta */
                    dn,directed,bipartite, 
                    nterms,funnames,
                    sonames, 
                    MHproposaltype,MHproposalpackage,
                    input_proposal,theta_proposal,&sample_size,
                    sample,burnin,interval,  
                    newnetworkheads, 
                    newnetworktails, 
                    verbose, 
                    attribs,maxout,maxin,minout,
                    minin,condAllDegExact,attriblength, 
                    maxedges,
                    status);
    if (print >= 0)
      {
      if (*status == 1) Rprintf("\nWARNING: Sample_Ls_Theta_Dependence: number of edges %i is outside of (1, %i).",newnetworkheads[0],*maxedges-1);
      else if (*status == 2) Rprintf("\nWARNING: M-H proposal failed.");
      }
    proposal_n_edges = newnetworkheads[0]; /* Number of simulated edges */
    proposal_heads = (int*) calloc(proposal_n_edges,sizeof(int)); /* Proposed heads for auxiliary variable */
    if (proposal_heads == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Dependence, proposal_heads\n\n"); error("Error: out of memory"); }
    proposal_tails = (int*) calloc(proposal_n_edges,sizeof(int)); /* Proposed tails for auxiliary variable */
    if (proposal_tails == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Dependence, proposal_tails\n\n"); error("Error: out of memory"); }
    for (i = 0; i < proposal_n_edges; i++)  
      {
      proposal_heads[i] = newnetworkheads[i+1]; /* Note: while heads corresponds to the list of observed heads, newnetworkheads contains the number of   simulated edges as well as the list of simulated heads: to use auxiliary->heads here, one must not store the number of simulated edges */
      proposal_tails[i] = newnetworktails[i+1]; /* Note: while tails corresponds to the list of observed tails, newnetworktails contains the number of   simulated edges as well as the list of simulated tails: to use auxiliary->tails here, one must not store the number of simulated edges */
      }
    /* Ratio of proposal pmfs of auxiliary graph under proposal / present */
    log_numerator = Minus_Energy(ergm->d,input_present,theta_present,
    proposal_heads,proposal_tails,&proposal_n_edges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n- function 1: log numerator = %8.4f",log_numerator);  
    */
    log_denominator = Minus_Energy(ergm->d,input_proposal,theta_proposal,
    proposal_heads,proposal_tails,&proposal_n_edges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n- function 2: log denominator = %8.4f",log_denominator);  
    */
    log_ratio = log_ratio + (log_numerator - log_denominator);
    /*
    Rprintf("\n- log ratio (auxiliary graph) = %8.4f",log_ratio);  
    */
    /* Ratio of mass of observed graph under proposal / present */
    log_present = Minus_Energy(ergm->d,input_present,theta_present,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n- function 3: log present = %8.4f",log_present);  
    */
    log_proposal = Minus_Energy(ergm->d,input_proposal,theta_proposal,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
    /*
    Rprintf("\n- function 4: log proposal = %8.4f",log_proposal);  
    */
    log_ratio = log_ratio + (log_proposal - log_present);
    /*
    Rprintf("\n- log ratio (observed graph) = %8.4f",log_ratio);  
    */
    free(proposal_heads);
    free(proposal_tails);
    }
  accept = MH_Decision(log_ratio);
  if ((*status == 0) && (accept == 1)) /* Proposal accepted */
    {
    Set_DD_DD(ls->d,ls->number+1,ls->theta,ls_theta);
    }
  if (print >= 1)
    {
    Rprintf("\nSample parameters of block %i of size %i:",update_block+1,ls->size[update_block]);
    if (auxiliary == 0) Rprintf("\n- M-H acceptance probability: %8.4f",Min(e(log_ratio),1.0));  
    else Rprintf("\n- auxiliary-variable M-H acceptance probability: %8.4f",Min(e(log_ratio),1.0));  
    Rprintf("\n- decision: %i",accept);
    }
  free(present);
  for (i = 0; i < ls->d; i++)
    {
    free(cf[i]);
    }
  free(cf);
  for (i = 0; i < ls->d; i++)
    {
    free(ls_theta[i]);
    }
  free(ls_theta);
  free(input_proposal);
  free(proposal);
  free(statistic);
  free(theta_present);
  free(theta_proposal);
  return accept;
}

int Sample_Ls_Theta_Between(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges,
                        int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames, char **sonames, 
                        double *input_present, int print,
                        double *scale_factor)
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/
{
  int accept, i, k, number_input;
  double **cf, *input_proposal, log_present, log_proposal, log_ratio, present_energy, present_a, proposal_energy, proposal_a, *present, *proposal, *theta, *statistic;
  number_input = Number_Input(ergm->terms,input_present);
  input_proposal = (double*) calloc(number_input,sizeof(double));
  if (input_proposal == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Between, input_proposal\n\n"); error("Error: out of memory"); }
  for (i = 0; i < number_input; i++) 
    { 
    input_proposal[i] = input_present[i];
    }
  present = (double*) calloc(ls->d,sizeof(double));  
  if (present == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Between, present\n\n"); error("Error: out of memory"); }
  statistic = (double*) calloc(ergm->d,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Between, statistic\n\n"); error("Error: out of memory"); }
  /* Generate candidate: */
  for (i = 0; i < ls->d; i++)
    {
    present[i] = ls->theta[i][ls->number];
    } 
  cf = (double**) calloc(ls->d,sizeof(double*));
  if (cf == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Between, cf\n\n"); error("Error: out of memory"); }
  for (i = 0; i < ls->d; i++)
    {
    cf[i] = (double*) calloc(ls->d,sizeof(double));
    }
  for (i = 0; i < ls->number_between; i++)
    {
    k = ls->between[i];
    cf[k][k] = scale_factor[0] * prior->cf2[k][k];
    }
  proposal = Sample_MVN(ls->d,present,cf); /* Random walk Metropolis-Hastings algorithm */
  /* Log acceptance ratio: */
  log_ratio = 0.0;
  log_proposal = MVN_PDF(ls->d,proposal,prior->mean2,prior->precision2); 
  log_present = MVN_PDF(ls->d,present,prior->mean2,prior->precision2); 
  log_ratio = log_ratio + (log_proposal - log_present); /* Log prior ratio */
  /*
  Rprintf("\n- log ratio: %8.4f",log_ratio);  
  */
  /* Log likelihood ratio: */
  theta = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input given present */
  present_energy = Minus_Energy(ergm->d,input_present,theta,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
  /*
  Rprintf("\n- present_energy = %8.4f",present_energy);  
  */
  present_a = Between_Block_Partition_Function(ls,ergm,input_present,theta,dn,directed,bipartite,nterms,funnames,sonames);
  /*
  Rprintf("\n- present_a = %8.4f",present_a);  
  */
  log_present = present_energy - present_a;
  /*
  Rprintf("\n- log present = %8.4f",log_present);  
  */
  for (i = 0; i < ls->d; i++)
    {
    ls->theta[i][ls->number] = proposal[i];
    } 
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_proposal); /* Set input given proposal */
  proposal_energy = Minus_Energy(ergm->d,input_proposal,theta,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,statistic); 
  /*
  Rprintf("\n- proposal_energy = %8.4f",proposal_energy);  
  */
  proposal_a = Between_Block_Partition_Function(ls,ergm,input_proposal,theta,dn,directed,bipartite,nterms,funnames,sonames);
  /*
  Rprintf("\n- proposal_a = %8.4f",proposal_a);  
  */
  for (i = 0; i < ls->d; i++)
    {
    ls->theta[i][ls->number] = present[i];
    } 
  log_proposal = proposal_energy - proposal_a;
  /*
  Rprintf("\n- log proposal = %8.4f",log_proposal);  
  */
  log_ratio = log_ratio + (log_proposal - log_present);
  /*
  Rprintf("\n- log ratio = %8.4f",log_ratio);  
  */
  accept = MH_Decision(log_ratio);
  if (accept == 1) /* Proposal accepted */
    {
    for (i = 0; i < ls->d; i++)
      {
      ls->theta[i][ls->number] = proposal[i];
      } 
    }
  if (print >= 1)
    {
    Rprintf("\nSample between-block parameters:");
    Rprintf("\n- M-H acceptance probability: %8.4f",Min(e(log_ratio),1.0));  
    Rprintf("\n- decision: %i",accept);
    }
  for (i = 0; i < ls->d; i++)
    {
    free(cf[i]);
    }
  free(cf);
  free(input_proposal);
  free(present);
  free(proposal);
  free(statistic);
  free(theta);
  return accept;
}

void Gibbs_Parameters(ergmstructure *ergm, latentstructure *ls, priorstructure *prior)
/*
input: ergm structure, latent structure, prior
output: non-structural parameters not showing up in the ergm pmf
*/
{
  int i;
  double *theta;
  for (i = 0; i < ls->number; i++)
    {
    if (ls->size[i] < ls->minimum_size) /* Structural parameter not showing up in ergm pmf */
      {
      theta = Sample_MVN(ls->d,prior->mean2,prior->cf2); /* Sample structural parameter from full conditional (conditional Gaussian prior given non-structural parameters) */
      Set_Column(ls->d,ls->theta,i,theta); /* Set ls_theta[][i] to theta */
      free(theta); 
      }
    }
}

double* Gibbs_Parameters_Means(priorstructure *prior, latentstructure *ls)
/*
input: prior structure, latent structure
output: means of parameters
*/
{
  int i, k;
  double sum, mean, numerator, denominator, precision, *sample, std;  
  sample = (double*) calloc(ls->d,sizeof(double));
  if (sample == NULL) { Rprintf("\n\ncalloc failed: Gibbs_Parameters_Means, sample\n\n"); error("Error: out of memory"); }
  for (i = 0; i < ls->d; i++)
    {
    sum = 0.0;
    for (k = 0; k < ls->number; k++) 
      {
      sum = sum + ls->theta[i][k];
      }
    numerator = (prior->mean2_precision[i] * prior->mean2_mean[i]) + (prior->precision2[i][i] * sum);
    denominator = prior->mean2_precision[i] + (ls->number * prior->precision2[i][i]);
    mean = numerator / denominator;
    precision = denominator;
    std = sqrt(1.0 / precision);
    sample[i] = mean + (norm_rand() * std); 
    /*
    Rprintf("\nprior->mean2_precision[%i] * prior->mean2_mean[%i] = %-8.4f",i,i,prior->mean2_precision[i] * prior->mean2_mean[i]);
    Rprintf("\nprior->precision2[%i][%i] * sum = %-8.4f",i,i,prior->precision2[i][i] * sum);
    Rprintf("\nnumerator = %-8.4f",numerator);
    Rprintf("\ndenominator = %-8.4f",denominator);
    Rprintf("\nmean of full conditional Gaussian of prior->mean2[%i] = %-8.4f",i,mean);
    Rprintf("\nprecision of full conditional Gaussian of prior->mean2[%i] = %-8.4f",i,precision);
    Rprintf("\ndraw from full conditional Gaussian of prior->mean2[%i] = %-8.4f",i,sample[i]);
    */
    }
  return sample;
}

double* Gibbs_Parameters_Means_Conditional(priorstructure *prior, latentstructure *ls)
/*
input: prior structure, latent structure
output: means of parameters
*/
{
  int i, k;
  /* PLEASE NOTE: FUNCTION NOT THOROUGHLY TESTED AND CURRENTLY NOT USED */
  double sum, mean, numerator, denominator, precision, *sample, std, tau;
  sample = (double*) calloc(ls->d,sizeof(double));
  if (sample == NULL) { Rprintf("\n\ncalloc failed: Gibbs_Parameters_Means_Conditional, sample\n\n"); error("Error: out of memory"); }
  tau = 1.0; /* Must match Gibbs_Parameters_Precisions_Marginal */
  for (i = 0; i < ls->d; i++)
    {
    sum = 0.0;
    for (k = 0; k < ls->number; k++) 
      {
      sum = sum + ls->theta[i][k];
      }
    numerator = (prior->mean2_precision[i] * prior->mean2_mean[i]) + sum;
    denominator = prior->mean2_precision[i] + ls->number;
    mean = numerator / denominator;
    precision = (prior->mean2_precision[i] + ls->number) * prior->precision2[i][i];
    std = sqrt(1.0 / precision);
    sample[i] = mean + (norm_rand() * std); 
    /*
    Rprintf("\nprior->mean2_precision[%i] * prior->mean2_mean[%i] = %-8.4f",i,i,prior->mean2_precision[i] * prior->mean2_mean[i]);
    Rprintf("\nprior->precision2[%i][%i] * sum = %-8.4f",i,i,prior->precision2[i][i] * sum);
    Rprintf("\nnumerator = %-8.4f",numerator);
    Rprintf("\ndenominator = %-8.4f",denominator);
    Rprintf("\nmean of full conditional Gaussian of prior->mean2[%i] = %-8.4f",i,mean);
    Rprintf("\nprecision of full conditional Gaussian of prior->mean2[%i] = %-8.4f",i,precision);
    Rprintf("\ndraw from full conditional Gaussian of prior->mean2[%i] = %-8.4f",i,sample[i]);
    */
    }
  return sample;
}

double* Gibbs_Parameters_Precisions(priorstructure *prior, latentstructure *ls)
/*
input: prior structure, latent structure
output: precisions of parameters
*/
{
  int i, k;
  double rate, shape, *sample, sum, d;
  sample = (double*) calloc(ls->d,sizeof(double)); 
  if (sample == NULL) { Rprintf("\n\ncalloc failed: Gibbs_Parameters_Precisions, sample\n\n"); error("Error: out of memory"); }
  for (i = 0; i < ls->d; i++)
    {
    shape = prior->precision2_shape + (ls->number / 2.0);
    sum = 0.0;
    for (k = 0; k < ls->number; k++)
      {
      d = ls->theta[i][k] - prior->mean2[i];
      sum = sum + (d * d);
      }
    rate = prior->precision2_rate + (sum / 2.0);
    sample[i] = rgamma(shape,1.0/rate); 
    /*
    Rprintf("\nshape of full conditional Gamma of prior->precision[%i] = %-8.4f",i,shape);
    Rprintf("\nrate of full conditional Gamma of prior->precision[%i] = %-8.4f",i,rate);
    Rprintf("\nmean of full conditional Gamma of prior->precision[%i] = %-8.4f",i,shape / rate);
    Rprintf("\nvariance of full conditional Gamma of prior->precision[%i] = %-8.4f",i,shape / (rate * rate));  
    Rprintf("\ndraw from full conditional Gaussian of prior->precision[%i] = %-8.4f",i,sample[i]);
    */
    } 
  return sample;
}

double* Gibbs_Parameters_Precisions_Marginal(priorstructure *prior, latentstructure *ls)
/*
input: prior structure, latent structure
output: precisions of parameters
*/
{
  int i, k;
  double d, m, rate, s, shape, *sample, t, tau, numerator, denominator;
  /* PLEASE NOTE: FUNCTION NOT THOROUGHLY TESTED AND CURRENTLY NOT USED */
  sample = (double*) calloc(ls->d,sizeof(double)); 
  if (sample == NULL) { Rprintf("\n\ncalloc failed: Gibbs_Parameters_Precisions_Marginal, sample\n\n"); error("Error: out of memory"); }
  tau = 1.0; /* Must match Gibbs_Parameters_Means_Conditional */
  tau = 1.0;
  for (i = 0; i < ls->d; i++)
    {
    shape = prior->precision2_shape + (ls->number / 2.0);
    m = 0.0;
    for (k = 0; k < ls->number; k++)
      {
      m = m + ls->theta[i][k];
      }
    m = m / ls->number;
    s = 0.0;
    for (k = 0; k < ls->number; k++)
      {
      d = ls->theta[i][k] - m;
      s = s + (d * d);
      }
    s = s / 2.0;
    d = m - prior->mean2[i];
    numerator = tau * ls->number * (d * d);
    denominator = 2.0 * (tau + ls->number);
    t = numerator / denominator;
    rate = prior->precision2_rate + s + t;
    sample[i] = rgamma(shape,1.0/rate); 
    /*
    Rprintf("\nshape of full conditional Gamma of prior->precision[%i] = %-8.4f",i,shape);
    Rprintf("\nrate of full conditional Gamma of prior->precision[%i] = %-8.4f",i,rate);
    Rprintf("\nmean of full conditional Gamma of prior->precision[%i] = %-8.4f",i,shape / rate);
    Rprintf("\nvariance of full conditional Gamma of prior->precision[%i] = %-8.4f",i,shape / (rate * rate));  
    Rprintf("\ndraw from full conditional Gaussian of prior->precision[%i] = %-8.4f",i,sample[i]);
    */
    } 
  return sample;
}

void Initial_State(int *parallel, double *alpha, int *indicator, priorstructure_ls *prior_ls, priorstructure *prior, latentstructure *ls, ergmstructure *ergm, double *theta, double *scale_factor)
/* 
input: clustering parameter, priors, latent structure, ergm structure, user-specified initial value of non-structural parameters
*/
{  
  int i, k;
  double *sample, *shape1, *shape2;
  if (*parallel == 1) ls->alpha = *alpha; /* Clustering parameter */
  else ls->alpha = rgamma(prior_ls->alpha_shape,1.0/prior_ls->alpha_rate); 
  shape1 = (double*) calloc(ls->number-1,sizeof(double)); /* Components 0..ls->number-2 suffice */
  if (shape1 == NULL) { Rprintf("\n\ncalloc failed: Initial_State, shape1\n\n"); error("Error: out of memory"); }
  shape2 = (double*) calloc(ls->number-1,sizeof(double)); /* Components 0..ls->number-2 suffice */
  if (shape2 == NULL) { Rprintf("\n\ncalloc failed: Initial_State, shape2\n\n"); error("Error: out of memory"); }
  for (i = 0; i < (ls->number - 1); i++)
    {
    shape1[i] = 1.0; /* First shape of Beta distribution */
    shape2[i] = ls->alpha; /* Second shape of Beta distribution */
    }
  sample = Stick_Breaking(shape1,shape2,ls); /* Construct category probability vector by stick-breaking */
  Set_D_D(ls->number,ls->p,sample);
  free(sample);
  free(shape1);
  free(shape2);
  for (i = 0; i < ls->n; i++) /* For each node i, sample category k */
    {
    k = ls->indicator[i]; 
    ls->size[k] = ls->size[k] + 1; /* ls-size was set to 0 by calloc */
    }
  if (ergm->d1 > 0)
    {
    sample = Sample_MVN(ergm->d1,prior->mean1,prior->cf1);
    Set_D_D(ergm->d1,ergm->theta,sample);
    free(sample);
    }
  for (i = 0; i < ls->number; i++) 
    {
    sample = Sample_MVN(ls->d,prior->mean2,prior->cf2); /* Random walk Metropolis-Hastings algorithm */
    Set_Column(ls->d,ls->theta,i,sample); /* Set ls_theta[][i] to proposal */
    free(sample);
    }
}

int Sample_CRP(latentstructure *ls)
/*
input: latent structure ls
output: partition of set of nodes drawn from Chinese restaurant process with scaling parameter ls->alpha
*/
{
  int i, k, number;
  double *p;
  for (i = 0; i < ls->number; i++)
    {
    ls->size[i] = 0;
    }
  p = (double*) calloc(ls->n,sizeof(double));
  if (p == NULL) { Rprintf("\n\ncalloc failed: Sample_CRP, p\n\n"); error("Error: out of memory"); }
  k = 0;
  ls->indicator[0] = k; /* First guest sits down at first table */
  ls->size[k] = ls->size[k] + 1; /* Increment number of guests at first table */
  number = 1; /* Number of occupied tables */
  for (i = 1; i < ls->n; i++) /* Guests 2, ..., n sit down at table... */
    {
    for (k = 0; k < number; k++)
      {
      p[k] = ls->size[k] / (i + ls->alpha); /* ...occupied table k with probability p[k] proportional ls->size[k]; note: if nodes are labeled 0..n-1, the denominator is i + ls->alpha rather than i - 1 + ls->alpha */
      }
    p[number] = ls->alpha / (i + ls->alpha); /* ...unoccupied table number+1 with probability p[number+1] proportional ls->alpha; note: if nodes are labeled 0..n-1, the denominator is i + ls->alpha rather than i - 1 + ls->alpha */
    for (k = (number + 1); k < ls->n; k++)
      {
      p[k] = 0.0;
      }
    k = Sample_Discrete(p); /* Sample one value of discrete random variable with PMF p by inverse transform method */
    ls->indicator[i] = k;
    ls->size[k] = ls->size[k] + 1;
    if (k == number) number = number + 1;
    }
  free(p);
  return number;
}

int Sample_Graph_Edge_Independence(latentstructure *ls, double *ln_p, int *heads, int *tails)
/*
input: latent structure; probability of edge between nodes i and j on log scale
output: graph sampled from PMF p and number of edges
*/
{
  int index, i, j, number_edges;
  double u;
  /* 
  Note 1: if i < j, edge (i, j) should be stored as (i, j)
  Note 2: i, j should be integers between 1 and n
  */  
  number_edges = 0;
  index = 0;
  for (i = 1; i < ls->n; i++)
    {
    for (j = i + 1; j < ls->n + 1; j++)
      {
      u = unif_rand(); /* Sample uniform[0,1] */
      if (ln(u) < ln_p[index]) 
        {
        number_edges = number_edges + 1;
        heads[number_edges] = i; /* number_edges is between 1..degrees of freedom; heads[0] is supposed to contain the number of edges */
        tails[number_edges] = j; /* number_edges is between 1..degrees of freedom; tails[0] is supposed to contain the number of edges */
        }
      index = index + 1;
      }
    }    
  heads[0] = number_edges; /* heads[0] is supposed to contain the number of edges */
  tails[0] = number_edges; /* tails[0] is supposed to contain the number of edges */
  return number_edges;
}

void Dirichlet(int *n,
               int *number,
               double *alpha,
               double *eta_mean,
               double *eta_sd,
               int *indicator,
               double *eta)
/*
input: R input
output: draw from truncated Dirichlet process prior:
- indicator of block membership
- parameters
*/
{ 
  int i;
  double *p, *shape1, *shape2;
  /**************/
  /* Initialize */
  /**************/
  GetRNGstate();
  epsilon = DBL_EPSILON;
  maximum = DBL_MAX;
  shape1 = (double*) calloc(*number-1,sizeof(double)); /* Components 0..ls->number-2 suffice */
  if (shape1 == NULL) { Rprintf("\n\ncalloc failed: Dirichlet, shape1\n\n"); error("Error: out of memory"); }
  shape2 = (double*) calloc(*number-1,sizeof(double)); /* Components 0..ls->number-2 suffice */
  if (shape2 == NULL) { Rprintf("\n\ncalloc failed: Dirichlet, shape2\n\n"); error("Error: out of memory"); }
  /********************************************/
  /* Sample truncated Dirichlet process prior */
  /********************************************/
  for (i = 0; i < (*number - 1); i++)
    {
    shape1[i] = 1.0; /* First shape of Beta distribution */
    shape2[i] = *alpha; /* Second shape of Beta distribution */
    }
  p = Stick_Breaking_External(shape1,shape2,*number,*n); /* Construct category probability vector by stick-breaking */ 
  /*
  Print_D(*number,p);
  */
  for (i = 0; i < *n; i++) indicator[i] = Sample_Discrete(p);
  for (i = 0; i < *number; i++) eta[i] = *eta_mean + (*eta_sd*norm_rand());
  /*
  Print_I(*n,indicator);
  Print_D(*number,eta);
  */
  /************/
  /* Finalize */
  /************/
  free(shape1);
  free(shape2);
  PutRNGstate();
}

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
             int *max_iterations, int *between, double *mean_between, int *output, double *mcmc, int *sample_heads, int *sample_tails, int *prior_assumptions, int *status)
/*
input: R input
output: simulated graph
*/
{
  int null = 0;
  int coordinate, degenerate_draws, dim, dim1, dim2, edges, element, h, i, j, dyad_dependence, *n_edges, *pseudo_indicator, iteration, k, max_iteration, *mdnedges, *mheads, *mtails, minimum_size, n, *newnetworkheads, *newnetworktails, number, print, threshold, terms, *verbose, parametric, hyper_prior, underflow;
  double between_edge_parameter, *draw, *p, **parameter, *pp, progress, *shape1, *shape2, prob_between, sum;	
  priorstructure_ls *prior_ls;
  latentstructure *ls;
  priorstructure *prior;
  ergmstructure *ergm;
  /**************/
  /* Initialize */
  /**************/
  GetRNGstate();
  print = *v; /* Console: no print; 0: short print; 1: long print */
  verbose = &null;
  epsilon = DBL_EPSILON;
  maximum = DBL_MAX;
  if (print >= 2)
    {
    Rprintf("\nMachine precision:");
    Rprintf("\n- epsilon = %e",epsilon);
    Rprintf("\n- maximum = %e",maximum);
    Rprintf("\n");
    }
  terms = (int)*nterms; /* Number of ergm terms */
  dim = (int)*d;
  dim1 = (int)*d1;
  dim2 = (int)*d2;
  n = (int)*dn; /* Number of nodes */
  number = (int)*max_number; /* Number of categories */
  max_iteration = (int)*max_iterations; /* Number of draws from posterior */
  dyad_dependence = (int)*dyaddependence; /* Conditional PMF of graph given latent structure: dyad-dependent or not */
  /*
  if (*directed == 1) 666666  
    {
    dyad_dependence = 1;
    }
  */
  minimum_size = (int)*min_size; /* Minimum size of category so that structural parameters show up in ergm pmf */
  ergm = Initialize_Ergm(terms,hierarchical,dim,dim1,dim2,structural); /* Ergm structure and non-structural parameters */
  prior = Initialize_Prior(ergm->d1,ergm->d2,m2_mean,m2_precision,*p2_shape,*p2_rate,m1,m2,cf1,cf2,p1,p2); /* Prior: non-structural, structural parameters */
  if (dyad_dependence == 0) threshold = n + 1; /* Minimum size of category so that structural parameters show up in ergm pmf */
  else 
    {
    if (*directed == 0) threshold = 8;
    else threshold = 6;
    }
  ls = Initialize_Latentstructure(number,n,indicator,minimum_size,threshold,ergm->d2,between,scaling); /* Latent structure and structural parameters */
  prior_ls = Initialize_Prior_ls(*alpha_shape,*alpha_rate); /* Prior: clustering parameter */
  mdnedges = &null;
  mheads = NULL;
  mtails = NULL;
  shape1 = (double*) calloc(ls->number-1,sizeof(double)); /* Components 0..ls->number-2 suffice */
  if (shape1 == NULL) { Rprintf("\n\ncalloc failed: Simulation, shape1\n\n"); error("Error: out of memory"); }
  shape2 = (double*) calloc(ls->number-1,sizeof(double)); /* Components 0..ls->number-2 suffice */
  if (shape2 == NULL) { Rprintf("\n\ncalloc failed: Simulation, shape2\n\n"); error("Error: out of memory"); }
  p = (double*) calloc(*maxpossibleedges,sizeof(double));
  if (p == NULL) { Rprintf("\n\ncalloc failed: Simulation, p\n\n"); error("Error: out of memory"); }
  parameter = (double**) calloc(ls->d,sizeof(double*));
  if (parameter == NULL) { Rprintf("\n\ncalloc failed: Simulation, parameter\n\n"); error("Error: out of memory"); }
  for (i = 0; i < ls->d; i++)
    {
    parameter[i] = (double*) calloc(ls->number+1,sizeof(double));
    if (parameter[i] == NULL) { Rprintf("\n\ncalloc failed: Simulation, parameter[%i]\n\n"); error("Error: out of memory"); }
    }
  pseudo_indicator = (int*) calloc(ls->n,sizeof(int));
  if (pseudo_indicator == NULL) { Rprintf("\n\ncalloc failed: Simulation, pseudo_indicator\n\n"); error("Error: out of memory"); }
  pp = (double*) calloc(ergm->d,sizeof(double));
  if (pp == NULL) { Rprintf("\n\ncalloc failed: Simulation, pp\n\n"); error("Error: out of memory"); }
  /* 
  Rprintf("\nh_ergm_c: d = %i and pp = %8.4f %8.4f \n",ergm->d, pp[0], pp[1]);
  */ 
  newnetworkheads = (int*) calloc(*maxpossibleedges+1,sizeof(int));
  if (newnetworkheads == NULL) { Rprintf("\n\ncalloc failed: Simulation, newnetworkheads\n\n"); error("Error: out of memory"); }
  newnetworktails = (int*) calloc(*maxpossibleedges+1,sizeof(int));
  if (newnetworktails == NULL) { Rprintf("\n\ncalloc failed: Simulation, newnetworktails\n\n"); error("Error: out of memory"); }
  hyper_prior = prior_assumptions[0]; /* If 1, hierarchical prior, otherwise non-hierarchical prior */
  parametric = prior_assumptions[1]; /* If 1, parametric prior */
  /****************/
  /* Sample graph */
  /****************/
  degenerate_draws = 0;
  k = -1;
  for (i = 0; i < ergm->d1; i++)
    {
    k = k + 1;
    ergm->theta[i] = eta[k];
    /*
    Rprintf("\nergm->theta[%i] = %8.4f",i,ergm->theta[i]);
    */
    }
  for (i = 0; i < ls->d; i++)
    {
    for (j = 0; j < (ls->number + 1); j++)
      {
      k = k + 1;
      ls->theta[i][j] = eta[k];
      /*
      Rprintf("\neta[%i] = %8.4f",k,eta[k]);
      Rprintf("\nls->theta[%i][%i] = %8.4f",i,j,ls->theta[i][j]);
      */
      }
    }
  coordinate = -1;
  element = -1;
  for (iteration = 0; iteration < max_iteration; iteration++)
    {
    progress = (iteration * 100.0) / max_iteration;
    if (print >= 1) Rprintf("\nProgress: %5.2f%%",progress);
    if (ls->number_fixed < ls->n)
      {	    
      if (*null_alpha == 1) ls->alpha = rgamma(prior_ls->alpha_shape,1.0/prior_ls->alpha_rate); 
      if (parametric == 0) ls->p = Stick_Breaking(shape1,shape2,ls); /* Construct category probability vector by stick-breaking */ 
      else Sample_Dirichlet(ls->number,ls->alpha,ls->p); /* Category probability vector */
      sum = 0.0;
      underflow = 0;
      for (i = 0; i < ls->number; i++)
         {
         if (ls->p[i] < epsilon) 
           {
           underflow = 1;
           ls->p[i] = epsilon;
           }
         sum = sum + ls->p[i];
         }
      if (underflow == 1) /* Renormalize */
        {
        for (i = 0; i < ls->number; i++) 
          {
          ls->p[i] = ls->p[i] / sum;
          }
        }
      for (i = 0; i < ls->n; i++)
        {
        if (ls->fixed[i] == 0) /* If block membership is not specified by user, sample block membership */
	  {
	  k = Sample_Discrete(ls->p);
          ls->indicator[i] = k;
	  } 
        }
      }
    if (*null_eta == 1)
      {
      if (ergm->d1 > 0) 
        {
        draw = Sample_MVN(ergm->d1,prior->mean1,prior->cf1);
        Set_D_D(ergm->d1,ergm->theta,draw);
        free(draw);
        }
      for (i = 0; i < ls->number; i++) 
        {
        draw = Sample_MVN(ls->d,prior->mean2,prior->cf2);
        Set_Column(ls->d,ls->theta,i,draw); /* Set ls_theta[][i] to proposal */ 
        free(draw);
        }
      for (i = 0; i < ls->d; i++)
        {
        if (between[i] == 1) 
          { 
          prob_between = *mean_between; /* Between-block probability of an edge as specified in hergm.preprocess() either by using user input or by using the expected number of edges <= 3 rule outline above*/ 
	  between_edge_parameter = log(prob_between / (1.0 - prob_between)); 
          ls->theta[i][ls->number] = between_edge_parameter; 
          } 
        }
      }
    for (k = 0; k < ls->number; k++)
      {
      ls->size[k] = 0;
      }
    for (i = 0; i < ls->n; i++)
      {
      k = ls->indicator[i];
      ls->size[k] = ls->size[k] + 1;
      }
    Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,inputs);
    Set_Parameter(ergm->d,ergm->structural,ergm->theta,theta); 
    for (i = 0; i < *maxpossibleedges; i++)
      {
      newnetworkheads[i] = 0;
      newnetworktails[i] = 0;
      }
    for (i = 0; i < ergm->d; i++)
      {
      pp[i] = 0.0;
      }
    if ((*dyaddependence == 0) && (*directed == 0)) /* Undirected dyad-independence ERGM */ 
      { 
      P_Edge_Independence(nterms,d,inputs,theta,dn,directed,bipartite,funnames,sonames,p);
      Sample_Graph_Edge_Independence(ls,p,newnetworkheads,newnetworktails);
      for (i = 0; i < ls->n; i++) /* Identical indicators */
        {
        pseudo_indicator[i] = 1;
        }
      for (i = 0; i < ls->d; i++) /* Identical structural parameters, so that structural function of graph, structural parameters reduces to corresponding non-structural function of graph */
        {
        for (k = 0; k < (ls->number + 1); k++)
          {
          parameter[i][k] = 1.0; 
          }
        }
      Set_Input(terms,hierarchical,number,ls->n,pseudo_indicator,parameter,inputs);
      n_edges = newnetworkheads; /* First element of newnetworkheads = newnetworkheads[0] is number of edges */
      mheads = (int*) calloc(*n_edges,sizeof(int));
      if (mheads == NULL) { Rprintf("\n\ncalloc failed: Simulation, mheads\n\n"); error("Error: out of memory"); }
      mtails = (int*) calloc(*n_edges,sizeof(int));
      if (mtails == NULL) { Rprintf("\n\ncalloc failed: Simulation, mtails\n\n"); error("Error: out of memory"); }
      for (i = 0; i < *n_edges; i++) /* Since first element of newnetworkheads and newnetworktails is number of edges, heads and tails must be extracted */
        {
        mheads[i] = newnetworkheads[i+1];
        mtails[i] = newnetworktails[i+1];
        }
      int timings = 0, time = 0, lasttoggle = 0;
      /*
      Rprintf("\n\nh_ergm.c:");
      Rprintf("\ntimings=%p",&timings);
      Rprintf("\ntimings=%i",timings);
      */
      network_stats_wrapper(mtails,mheads,&timings,&time,&lasttoggle,n_edges,dn,directed,bipartite,nterms,funnames,sonames,inputs,pp); /* Compute non-structural function of graph */
      free(mheads);
      free(mtails);
      }
    else Sample_Graph(ls->number,ls->n,ls->d,ergm->terms,ergm->hierarchical,ergm->d,pp, /* Sample graph */
                         heads,tails,dnedges,maxpossibleedges,dn,directed,bipartite,nterms,funnames,
                         sonames,MHproposaltype,MHproposalpackage,inputs,theta,samplesize,
                         sample,burnin,interval,newnetworkheads,newnetworktails,verbose,
                         attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,
                         maxedges,mheads,mtails,mdnedges,status);
    degenerate_draws = degenerate_draws + *status;
    if (print >= 0)
      {
      if (*status == 1) Rprintf("\nWARNING: number of edges %i is outside of (1, %i).",newnetworkheads[0],*maxedges-1);
      else if (*status == 2) Rprintf("\nWARNING: Simulation: M-H proposal failed.");
      }
    if (ergm->d1 > 0)
      {
      if (print >= 1) Rprintf("\nparameters:");
      for (i = 0; i < ergm->d1; i++) /* Non-structural parameters */
        {
        if (print >= 1) Rprintf(" %8.4f",ergm->theta[i]);
        coordinate = coordinate + 1;
        mcmc[coordinate] = ergm->theta[i];
        }  
      }
    if (ergm->d2 > 0)
      {
      if ((print >= 1) && (*null_eta == 1)) Rprintf("\nmeans of block parameters:");
      for (i = 0; i < ls->d; i++) /* Structural parameters */
        {
        if ((print >= 1) && (*null_eta == 1))  Rprintf(" %8.4f",prior->mean2[i]);
        coordinate = coordinate + 1;	
        mcmc[coordinate] = prior->mean2[i];
        }
      if ((print >= 1) && (*null_eta == 1))  Rprintf("\nprecisions of block parameters:");
      for (i = 0; i < ls->d; i++) /* Structural parameters */
        {
        if ((print >= 1) && (*null_eta == 1))  Rprintf(" %8.4f",prior->precision2[i][i]);
        coordinate = coordinate + 1;	
        mcmc[coordinate] = prior->precision2[i][i];
        }
      if (print >= 1) 
        {
        Rprintf("\nblock parameters:");
        if (ergm->d2 > 1) Rprintf("\n");
        }
      for (h = 0; h < ls->d; h++) /* Structural parameters */
        {
        for (i = 0; i < ls->number; i++) 
          {
          if (print >= 1) Rprintf(" %8.4f",ls->theta[h][i]);
          coordinate = coordinate + 1;	
          mcmc[coordinate] = ls->theta[h][i];
          }
        if ((print >= 1) && (ls->number_between > 0)) Rprintf(" %8.4f",ls->theta[h][ls->number]); /* Second condition ensures that between-category parameters are not written to screen when model without between-category parameters */
        coordinate = coordinate + 1;	
        mcmc[coordinate] = ls->theta[h][ls->number];
        if (print >= 1) Rprintf("\n");
        }
      if (print >= 1) Rprintf("block indicators:");
      for (i = 0; i < ls->n; i++) /* Category indicators */
        {
        if (print >= 1) 
          {
          Rprintf(" %i",ls->indicator[i]+1);
          }
        coordinate = coordinate + 1;
        mcmc[coordinate] = ls->indicator[i];
        }
      if (print >= 1) Rprintf("\nblock sizes:");  
      for (i = 0; i < ls->number; i++)
        {
        if (print >= 1) Rprintf(" %3i",ls->size[i]);
        coordinate = coordinate + 1;
        mcmc[coordinate] = ls->size[i];
        } 
      if ((print >= 1) && (ls->number_fixed < ls->n)) Rprintf("\nblock probabilities:");
      for (i = 0; i < ls->number; i++) /* Category probability vector */
        {
        if ((print >= 1) && (ls->number_fixed < ls->n))Rprintf(" %6.4f",ls->p[i]);
        coordinate = coordinate + 1;
        mcmc[coordinate] = ls->p[i];
        }
      if ((print >= 1) && (ls->number_fixed < ls->n)) Rprintf("\nblock probabilities prior parameter: %6.4f",ls->alpha); /* Clustering parameter */
      coordinate = coordinate + 1;
      mcmc[coordinate] = ls->alpha;
      } 
    if (print >= 1) Rprintf("\nstatistics:");
    for (i = 0; i < ergm->d; i++) /* Statistics */
      {
      if (print >= 1) Rprintf(" %6.0f",pp[i]);
      coordinate = coordinate + 1;
      mcmc[coordinate] = pp[i];
      }
    if (*output == 1)
      { 
      element = element + 1;
      edges = newnetworkheads[0];
      sample_heads[element] = edges;
      sample_tails[element] = edges;
      for (i = 0; i < edges; i++)
        {
        element = element + 1;
        sample_heads[element] = newnetworkheads[i+1];
        sample_tails[element] = newnetworktails[i+1];
        }
      }
    if (print >= 1) Rprintf("\n");
    }
  /************/
  /* Finalize */
  /************/
  if (degenerate_draws > 0) Rprintf("\nWARNING: %i generated networks were extreme in terms of the number of edges. The corresponding draws should be discarded.\n\n", degenerate_draws);
  free(p);
  free(pseudo_indicator);
  for (i = 0; i < ls->d; i++)
    {
    free(parameter[i]);
    }
  free(parameter);
  free(newnetworkheads);
  free(newnetworktails);
  free(pp);
  free(shape1);
  free(shape2);
  Finalize_Ergm(ergm);
  Finalize_Latentstructure(ls,dim2);
  Finalize_Prior_ls(prior_ls);
  Finalize_Priorstructure(prior,dim1,dim2);
  PutRNGstate();
}

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
             int *max_iterations, int *between, int *output, double *mcmc, double *scalefactor, double *mh_accept, double *q_i, int *parallel, double *temperature, int *prior_assumptions, int *status)
/*
input: R input
output: MCMC sample of unknowns from posterior
*/
{
  int null = 0;
  int *available, degenerate_draws, model, batch, n_batches, batch_size, coordinate, console, dyad_dependence, dim, dim1, dim2, h, i, parametric, hyper_prior, *mdnedges, *mheads, *mtails, number_local_mh, iteration, max_iteration, minimum_size, n, number, print, store, threshold, terms, *verbose, underflow, update_block, update_node, update_size;
  double ls_alpha, accept, *local_mh, *local_mh_accept, *ls_p, *pp, *prior_mean2, *prior_precision2, progress, *scale_factor, sample_size_indicators, sum;	
  priorstructure_ls *prior_ls;
  latentstructure *ls;
  priorstructure *prior;
  ergmstructure *ergm;
  /**************/
  /* Initialize */
  /**************/
  GetRNGstate();
  console = *v; /* Console: -1: no print; 0: short print; 1: long print */
  verbose = &null;
  epsilon = DBL_EPSILON;
  maximum = DBL_MAX;
  if (console >= 2)
    {
    Rprintf("\nMachine precision:");
    Rprintf("\n- epsilon = %e",epsilon);
    Rprintf("\n- maximum = %e",maximum);
    Rprintf("\n");
    }
  model = (int)*model_type;
  terms = (int)*nterms; 
  dim = (int)*d;
  dim1 = (int)*d1;
  dim2 = (int)*d2;
  number = (int)*max_number; /* Number of categories */
  n = (int)*dn; /* Number of nodes */
  max_iteration = *max_iterations; /* Number of draws from posterior */
  if (max_iteration <= 10000) n_batches = max_iteration;
  else n_batches = 10000;
  if (max_iteration == n_batches) batch_size = 1; 
  else batch_size = trunc(max_iteration / n_batches);
  mdnedges = &null;
  ergm = Initialize_Ergm(terms,hierarchical,dim,dim1,dim2,structural); /* Ergm structure and non-structural parameters */
  prior = Initialize_Prior(ergm->d1,ergm->d2,m2_mean,m2_precision,*p2_shape,*p2_rate,m1,m2,cf1,cf2,p1,p2); /* Prior: non-structural, structural parameters */
  dyad_dependence = (int)*dyaddependence; /* Conditional PMF of graph given latent structure: dyad-dependent or not */
  /*
  if (*directed == 1) 666666 
    {
    dyad_dependence = 1;
    }
  */
  minimum_size = (int)*min_size; /* Minimum size of category so that structural parameters show up in ergm pmf */
  if (dyad_dependence == 0) threshold = n + 1; /* Minimum size of category so that structural parameters show up in ergm pmf */
  else 
    {
    if (*directed == 0) threshold = 6;
    else threshold = 5;
    }
  /*
  Print_I(n,indicator);
  */
  ls = Initialize_Latentstructure(number,n,indicator,minimum_size,threshold,ergm->d2,between,scaling); /* Latent structure and structural parameters */
  /*
  Print_I(n,ls->indicator);
  */
  prior_ls = Initialize_Prior_ls(*alpha_shape,*alpha_rate); /* Prior: clustering parameter */
  hyper_prior = prior_assumptions[0]; /* If 1, hierarchical prior, otherwise non-hierarchical prior */
  parametric = prior_assumptions[1]; /* If 1, parametric prior */
  mheads = NULL;
  mtails = NULL;
  pp = (double*) calloc(ergm->d,sizeof(double));
  if (pp == NULL) { Rprintf("\n\ncalloc failed: Inference, pp\n\n"); error("Error: out of memory"); }
  scale_factor = (double*) calloc(2+ls->n-ls->minimum_size,sizeof(double));   
  if (scale_factor == NULL) { Rprintf("\n\ncalloc failed: Inference, scale_factor\n\n"); error("Error: out of memory"); }
  /*************************/
  /* MCMC sample posterior */
  /*************************/
  if (console >= 1)
    {
    Rprintf("\nNumber of draws from posterior: %i",n_batches * batch_size);
    Rprintf("\nNumber of batches: %i",n_batches);
    Rprintf("\nSize of batches: %i",batch_size);
    Rprintf("\n");
    }
  if (ls->number == 1) 
    {
    for (i = 0; i < ls->n; i++)
      {
      ls->indicator[i] = 0;
      }
    ls->size[0] = ls->n;
    }
  else Initial_State(parallel,alpha,ls->indicator,prior_ls,prior,ls,ergm,theta,scale_factor);
  if (dyad_dependence == 0) number_local_mh = 2;
  else number_local_mh = 2 + (ls->n - ls->minimum_size);
  local_mh = (double*) calloc(number_local_mh,sizeof(double));
  if (local_mh == NULL) { Rprintf("\n\ncalloc failed: Inference, local_mh\n\n"); error("Error: out of memory"); }
  local_mh_accept = (double*) calloc(number_local_mh,sizeof(double));
  if (local_mh_accept == NULL) { Rprintf("\n\ncalloc failed: Inference, local_mh_accept\n\n"); error("Error: out of memory"); }
  /* 
  Rprintf("\n\nls->n = %i ls->minimum_size = %i number_local_mh = %i",ls->n,ls->minimum_size,number_local_mh);
  Print_D(number_local_mh,scale_factor);
  Print_D(number_local_mh,scalefactor);
  */
  Set_D_D(number_local_mh,scale_factor,scalefactor); /* Metropolis-Hasting algorithm: scale factor */
  if (ls->n < 20) sample_size_indicators = round(ls->n / 3);
  else if (ls->n < 50) sample_size_indicators = round(ls->n / 6);
  else sample_size_indicators = round(ls->n / 12);
  if (sample_size_indicators > (ls->n - ls->number_fixed)) sample_size_indicators = ls->n - ls->number_fixed;
  degenerate_draws = 0;
  coordinate = -1;
  if (console == 0) 
    {
    Rprintf("\nSizes of blocks");
    for (i = 0; i < ls->number; i++) 
      {
      Rprintf(" %i", i+1);
      }
    Rprintf("\n");
    }
  for (batch = 0; batch < n_batches; batch++) /* Batch */
    {
    progress = (batch * 100.0) / n_batches;
    if (console >= 1) Rprintf("\nProgress: %5.2f%% of %i",progress,max_iteration);
    else if (console == 0) Rprintf("%4.1f%%",progress);
    for (iteration = 0; iteration < batch_size; iteration++) /* Iteration within batch */
      {
      if (iteration == 0)
        { 
        store = 1;
        print = console;
        }
      else 
        {
        store = 0;
        print = -1;
        }
      /* Generate MCMC sample: */
      if (dyad_dependence == 0) /* MCMC exploiting dyad-independence conditional on latent structure */
        {
        Gibbs_Indicators_Independence(ls,ergm,heads,tails,inputs_h,dnedges,dn,directed,bipartite,nterms,funnames,sonames,q_i); 
        if (ergm->d1 > 0)
          {
          local_mh[0] = local_mh[0] + 1;
          local_mh_accept[0] = local_mh_accept[0] + Sample_Ergm_Theta_Independence(ergm,ls,prior,heads,tails,dnedges,dn,directed,bipartite, 
                                   nterms,funnames,sonames,inputs,print,scale_factor);
          }
        accept = Sample_Ls_Theta_Independence(ergm,ls,prior,heads,tails,dnedges,dn,directed,bipartite,nterms,
                                   funnames,sonames,inputs,inputs_h,print,scale_factor); /* M-H exploiting dyad-independence conditional on latent structure */
        for (i = 1; i < number_local_mh; i++)
          {
          local_mh[i] = local_mh[i] + 1; 
          local_mh_accept[i] = local_mh_accept[i] + accept;
          }
        }
      else /* Dyad-dependence conditional on latent structure */
        {
        if ((ls->number > 1) && (ls->n - ls->number_fixed > 0)) /* More than one block and not all block memberships are fixed */
          {
	  available = (int*) calloc(ls->n,sizeof(int));
          if (available == NULL) { Rprintf("\n\ncalloc failed: Inference, available\n\n"); error("Error: out of memory"); }
          Set_I_I(ls->n,available,ls->fixed); /* The fixed elements of ls->indicator are not available for updating, so that ls->fixed[i] = 1 implies available[i] = 1 */
	  for (i = 0; i < sample_size_indicators; i++)
            {
            do update_node = trunc(unif_rand() * ls->n); /* Sample indicator of node to be updated */
            while (available[update_node] == 1);
            available[update_node] = 1;
            Sample_Indicators_Dependence(model,ergm,ls,prior, /* Auxiliary-variable M-H */
                         heads,tails,dnedges,maxpossibleedges,dn,directed,bipartite,nterms,funnames,
                         sonames,MHproposaltype,MHproposalpackage,sample,burnin,interval, 
                         verbose,attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,
                         maxedges,mheads,mtails,mdnedges,inputs,print,newnetworkheads,newnetworktails,scale_factor,update_node,temperature,status);
            degenerate_draws = degenerate_draws + *status;
            }
          free(available);
	  }
        if (ergm->d1 > 0)
          {
          local_mh[0] = local_mh[0] + 1; 
          local_mh_accept[0] = local_mh_accept[0] + Sample_Ergm_Theta_Dependence(model,ergm,ls,prior, /* Auxiliary-variable M-H */
                         heads,tails,dnedges,maxpossibleedges,dn,directed,bipartite,nterms,funnames,
                         sonames,MHproposaltype,MHproposalpackage,sample,burnin,interval, 
                         verbose,attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,
                         maxedges,mheads,mtails,mdnedges,inputs,print,newnetworkheads,newnetworktails,scale_factor,status);
          degenerate_draws = degenerate_draws + *status;
          }
        if (ergm->d2 > 0)
	  {
	  for (i = 0; i < ls->number; i++)
            {
            if (ls->size[i] >= ls->minimum_size) 
              { 
              update_block = i; /* Parameter vector of block to be updated */
              update_size = ls->size[update_block];
              local_mh[update_size-ls->minimum_size+1] = local_mh[update_size-ls->minimum_size+1] + 1; 
              local_mh_accept[update_size-ls->minimum_size+1] = local_mh_accept[update_size-ls->minimum_size+1] + Sample_Ls_Theta_Dependence(model,ergm,ls,prior, /* Auxiliary-variable M-H */
                         heads,tails,dnedges,maxpossibleedges,dn,directed,bipartite,nterms,funnames,
                         sonames,MHproposaltype,MHproposalpackage,sample,burnin,interval, 
                         verbose,attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,
                         maxedges,mheads,mtails,mdnedges,inputs,print,newnetworkheads,newnetworktails,scale_factor,update_block,status);
	      degenerate_draws = degenerate_draws + *status;
              }
            }
	  }
        }
      if (ls->number_between > 0)
        {
        local_mh[0] = local_mh[0] + 1;
        local_mh_accept[0] = local_mh_accept[0] + Sample_Ls_Theta_Between(ergm,ls,prior,
                              heads,tails,dnedges,dn,directed,bipartite,
                              nterms,funnames,sonames,inputs,print,scale_factor);
        }
      Gibbs_Parameters(ergm,ls,prior); /* Structural parameters not showing up in ergm pmf */
      if (parametric == 0)
         {
         ls_p = Sample_P(ls); /* Category probability vector */
	 Set_D_D(ls->number,ls->p,ls_p);
         free(ls_p);
         }
      else Sample_Dirichlet(ls->number,ls->alpha,ls->p); /* Category probability vector */ 
      sum = 0.0;
      underflow = 0;
      for (i = 0; i < ls->number; i++)
         {
         if (ls->p[i] < epsilon) 
           {
           underflow = 1;
           ls->p[i] = epsilon;
           }
         sum = sum + ls->p[i];
         }
      if (underflow == 1) /* Renormalize */
        {
        for (i = 0; i < ls->number; i++) 
          {
          ls->p[i] = ls->p[i] / sum;
          }
        }
      if (hyper_prior == 1) /* Hyper prior: mean and precisions of Gaussian baseline distribution have non-degenerate prior */
        {
        ls_alpha = Sample_Alpha(prior_ls,ls); /* Clustering parameter */
        ls->alpha = ls_alpha;
        /*
        prior_mean2 = Gibbs_Parameters_Means_Conditional(prior,ls); 
        */
        prior_mean2 = Gibbs_Parameters_Means(prior,ls);  
        Set_D_D(ls->d,prior->mean2,prior_mean2);
        /*
        prior_precision2 = Gibbs_Parameters_Precisions_Marginal(prior,ls); 
        */
        prior_precision2 = Gibbs_Parameters_Precisions(prior,ls); 
        for (i = 0; i < ls->d; i++)
          {
          prior->precision2[i][i] = prior_precision2[i];
          }
        free(prior_mean2);
        free(prior_precision2); 
        }
      /* Internal storage and console output: */
      if (store == 1) 
        {
        /* Posterior prediction when full output is desired */
        if (*output == 1)
          {
          Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,inputs);
          Set_Parameter(ergm->d,ergm->structural,ergm->theta,theta);
          for (i = 0; i < *maxedges; i++)
            {
            newnetworkheads[i] = 0;
            newnetworktails[i] = 0;
            }
          for (i = 0; i < ergm->d; i++)
            {
            pp[i] = 0.0;
            }
          Sample_Graph(ls->number,ls->n,ls->d,ergm->terms,ergm->hierarchical,ergm->d,pp, /* Posterior prediction of graph */
                           heads,tails,dnedges,maxpossibleedges,dn,directed,bipartite,nterms,funnames,
                           sonames,MHproposaltype,MHproposalpackage,inputs,theta,samplesize,
                           sample,burnin,interval,newnetworkheads,newnetworktails,verbose,
                           attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,
                           maxedges,mheads,mtails,mdnedges,status);
	  degenerate_draws = degenerate_draws + *status;
          }
        /* Store and output MCMC sample: */ 
        if (ergm->d1 > 0)
          {
          if (print >= 1) Rprintf("\nparameters:");
          for (i = 0; i < ergm->d1; i++) /* Non-structural parameters */
            {
            if (print >= 1) Rprintf(" %8.4f",ergm->theta[i]);
            coordinate = coordinate + 1;
            mcmc[coordinate] = ergm->theta[i];
            }  
          }
        if (ergm->d2 > 0)
          {
          if ((hyper_prior == 1) && (print >= 1)) Rprintf("\nmeans of block parameters:");
          for (i = 0; i < ls->d; i++) /* Structural parameters */
            {
            if ((hyper_prior == 1) && (print >= 1)) Rprintf(" %8.4f",prior->mean2[i]);
            coordinate = coordinate + 1;	
            mcmc[coordinate] = prior->mean2[i];
            }
          if ((hyper_prior == 1) && (print >= 1)) Rprintf("\nprecisions of block parameters:");
          for (i = 0; i < ls->d; i++) /* Structural parameters */
            {
            if ((hyper_prior == 1) && (print >= 1)) Rprintf(" %8.4f",prior->precision2[i][i]);
            coordinate = coordinate + 1;	
            mcmc[coordinate] = prior->precision2[i][i];
            }
	  if (print >= 1) 
            {
            Rprintf("\nblock parameters:");
            if (ergm->d2 > 1) Rprintf("\n");
            }
          for (h = 0; h < ls->d; h++) /* Structural parameters */
            {
            for (i = 0; i < ls->number; i++) 
              {
              if (print >= 1) Rprintf(" %8.4f",ls->theta[h][i]);
              coordinate = coordinate + 1;	
              mcmc[coordinate] = ls->theta[h][i];
              }
            if ((print >= 1) && (ls->number_between > 0)) Rprintf(" %8.4f",ls->theta[h][ls->number]); /* Second condition ensures that between-category parameters are not written to screen when model without between-category parameters */
            coordinate = coordinate + 1;	
            mcmc[coordinate] = ls->theta[h][ls->number];
            if (print >= 1) Rprintf("\n");
            }
          if (print >= 1) Rprintf("block indicators:");
          for (i = 0; i < ls->n; i++) /* Category indicators */
            {
            if (print >= 1) 
              {
     	      Rprintf(" %i",ls->indicator[i]+1);
              }
            coordinate = coordinate + 1; 
	    mcmc[coordinate] = ls->indicator[i];
            }
          if (print >= 1) Rprintf("\nblock sizes:");  
          for (i = 0; i < ls->number; i++)
            {
            if (print >= 1) Rprintf(" %3i",ls->size[i]);
            coordinate = coordinate + 1;
            mcmc[coordinate] = ls->size[i];
            }
          if (print >= 1) Rprintf("\nblock probabilities:");
          for (i = 0; i < ls->number; i++) /* Category probability vector */
            {
            if (print >= 1) Rprintf(" %6.4f",ls->p[i]);
            coordinate = coordinate + 1;
             mcmc[coordinate] = ls->p[i];
            }
          if ((hyper_prior == 1) && (print >= 1)) Rprintf("\nblock probabilities prior parameter: %6.4f",ls->alpha); /* Clustering parameter */
          coordinate = coordinate + 1;
          mcmc[coordinate] = ls->alpha;
          }
        if (*output == 1)
          {
          if ((*output == 1) && (print >= 1)) Rprintf("\nposterior prediction of statistics:");
          for (i = 0; i < ergm->d; i++) /* Posterior prediction of statistics */
            {
            if ((*output == 1) && (print >= 1)) Rprintf(" %6.0f",pp[i]);
            coordinate = coordinate + 1;
            mcmc[coordinate] = pp[i];
            }
          }
        if (print >= 1) Rprintf("\n");
        else if (print == 0) 
          {
          for (i = 0; i < ls->number; i++) Rprintf("%4i",ls->size[i]);
          if (*output == 1)
            {
            for (i = 0; i < ergm->d; i++) Rprintf("%6.0f",pp[i]); /* Posterior prediction of statistics */
            } 
          Rprintf("\n");
          }
        }
      }
    }
  /************/
  /* Finalize */
  /************/
  for (i = 0; i < number_local_mh; i++)
    {
    if (local_mh[i] == 0) local_mh_accept[i] = 1.0;
    else local_mh_accept[i] = local_mh_accept[i] / local_mh[i];
    }
  Set_D_D(number_local_mh,mh_accept,local_mh_accept);
  if (console >= 1)
    {
    Rprintf("\nNumber of draws from posterior: %i",n_batches * batch_size);
    Rprintf("\nThinning: every %i-th draw recorded\n\n",batch_size);
    if (degenerate_draws > 0) Rprintf("\nWARNING: %i generated networks were extreme in terms of the number of edges. The corresponding draws should be discarded.\n\n", degenerate_draws);
    }
  free(scale_factor);
  free(local_mh_accept);
  free(local_mh);
  free(pp);
  free(mheads);
  free(mtails);
  Finalize_Ergm(ergm);
  Finalize_Latentstructure(ls,dim2);
  Finalize_Prior_ls(prior_ls);
  Finalize_Priorstructure(prior,dim1,dim2);
  PutRNGstate();
}

