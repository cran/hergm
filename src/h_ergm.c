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
  if (b == NULL) { Rprintf("\n\ncalloc failed: Stick_Breaking, b\n\n"); exit(1); }
  p = (double*) calloc(ls->number,sizeof(double));
  if (p == NULL) { Rprintf("\n\ncalloc failed: Stick_Breaking, p\n\n"); exit(1); }
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

double* Sample_P(latentstructure *ls)
/*
input: latent structure
output: category probability vector
*/
{
  int i, rest;
  double *p, *shape1, *shape2;
  shape1 = (double*) calloc(ls->number-1,sizeof(double)); /* Element 0..ls->number-2 required */
  if (shape1 == NULL) { Rprintf("\n\ncalloc failed: Sample_P, sample1\n\n"); exit(1); }
  shape2 = (double*) calloc(ls->number-1,sizeof(double)); /* Element 0..ls->number-2 required */
  if (shape2 == NULL) { Rprintf("\n\ncalloc failed: Sample_P, sample2\n\n"); exit(1); }
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

void P_Edge_Independence(int *number_terms, int *number_parameters, double *input, double *theta,  int *n, int *directed, int *bipartite, char **funnames, char **sonames, double *p)
/*
input: undirected graph; number of terms; number of parameters;  input vector; parameter vector; number of nodes; other variables
output: probabilities of edges between nodes i and j on log scale, computed under the assumption of conditional edge-independence given latent structure,
and ordered in accordance with i < j
*/
{
  int one = 1;
  int h, i, j, *number_edges, *heads, *tails;
  double log_odds, *statistic;
  number_edges = &one;
  statistic = (double*) calloc(*number_parameters,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: P_Edge_Independence, statistic\n\n"); exit(1); }
  /* 
  Note 1: if undirected graph and i < j, undirected edge (i, j) should be stored as (i, j)
  Note 2: i, j should be integers between 1 and n
  */  
  h = -1;
  for (i = 1; i < *n + 1; i++) /* Row i */
    {
    heads = &i; 
    for (j = i + 1; j < *n + 1; j++) /* Row i, column j > i (undirected, directed graph) */
      {
      h = h + 1;
      tails = &j;
      log_odds = Minus_Energy(*number_parameters,input,theta,heads,tails,number_edges,n,directed,bipartite,number_terms,funnames,sonames,statistic); /* Compute log-odds of probability of edge statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
      p[h] = -ln(1.0 + e(-log_odds));
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
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: Partition_Function_Edge_Independence, statistic\n\n"); exit(1); }
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
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: Partition_Function_Dyad_Independence, statistic\n\n"); exit(1); }
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
      if (heads == NULL) { Rprintf("\n\ncalloc failed: Partition_Function_Dyad_Independence, heads 1\n\n"); exit(1); }
      tails = (int*) calloc(*number_edges,sizeof(int)); 
      if (tails == NULL) { Rprintf("\n\ncalloc failed: Partition_Function_Dyad_Independence, tails 1\n\n"); exit(1); }
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
      if (heads == NULL) { Rprintf("\n\ncalloc failed: Partition_Function_Dyad_Independence, heads 2\n\n"); exit(1); }
      tails = (int*) calloc(*number_edges,sizeof(int)); 
      if (tails == NULL) { Rprintf("\n\ncalloc failed: Partition_Function_Dyad_Independence, tails 2\n\n"); exit(1); }
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
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: PMF_Independence, statistic\n\n"); exit(1); }
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
  if (p_i == NULL) { Rprintf("\n\ncalloc failed: Gibbs_Indicators_Independence, p_i\n\n"); exit(1); }
  sample = (int*) calloc(ls->n,sizeof(int));
  if (sample == NULL) { Rprintf("\n\ncalloc failed: Gibbs_Indicators_Independence, sample\n\n"); exit(1); }
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

int Sample_Ergm_Theta_Independence(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
                        int *heads, int *tails, int *dnedges, int *dn, int *directed, int *bipartite, 
                        int *nterms, char **funnames, char **sonames, 
                        double *input, int print, int n_between, double scale_factor)
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/
{
  int accept, i;
  double **cf, *ergm_theta, log_present, log_proposal, log_ratio, *theta_present, *theta_proposal;
  log_ratio = 0.0;
  /* Propose ergm->theta: random walk (ratio of proposal pdfs cancels): */
  cf = Scale(ergm->d1,ergm->d1,prior->cf1,scale_factor); /* Rescale Cholesky factor of Gaussian prior */        
  ergm_theta = Sample_MVN(ergm->d1,ergm->theta,cf); /* Random walk Metropolis-Hastings */
  log_proposal = MVN_PDF(ergm->d1,ergm_theta,prior->mean1,prior->precision1); /* Prior pdf: proposal */
  log_present = MVN_PDF(ergm->d1,ergm->theta,prior->mean1,prior->precision1); /* Prior pdf: present */
  log_ratio = log_ratio + (log_proposal - log_present);
  /*
  Rprintf("\n- log_ratio (parameters) = %8.4f",log_ratio);  
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
  if (print == 1)
    {
    Rprintf("\nSample parameters:");
    Rprintf("\n- log_ratio = %8.4f",log_ratio);  
    Rprintf("\n- decision = %i",accept);
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
                        int *nterms, char **funnames, char **sonames, double *input_proposal, double *input_present, int print, int n_between, double scale_factor)
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/
{
  int accept, i;
  double **cf, *present, log_present, log_proposal, log_ratio, **ls_theta, *proposal, *theta;
  /* Proposal:
  note 1: all ls->theta such that ls->size >= ls->threshold and all ergm->theta are updated by random walk Metropolis-Hastings algorithm
  note 2: ratio of proposal pdfs cancels under random walk Metropolis-Hastings algorithm */
  log_ratio = 0.0;
  /* Propose ls->theta: */
  ls_theta = (double**) calloc(ls->d,sizeof(double*)); /* Remark: since memory is allocated to ls_theta by calloc, between-block parameters are 0 */
  if (ls_theta == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Independence, ls_theta\n\n"); exit(1); }
  for (i = 0; i < ls->d; i++)
    {
    ls_theta[i] = (double*) calloc(ls->number+1,sizeof(double));
    if (ls_theta[i] == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Independence, ls_theta[%i]\n\n",i); exit(1); }
    }
  present = (double*) calloc(ls->d,sizeof(double));  
  if (present == NULL) { Rprintf("\n\ncalloc failed: Sample_Ls_Theta_Independence, present\n\n"); exit(1); }
  cf = Scale(ls->d,ls->d,prior->cf2,scale_factor); /* Rescale Cholesky factor of Gaussian prior */ 
  for (i = 0; i < ls->number; i++) 
    {
    Get_Column(ls->d,present,ls->theta,i); /* Set mean to ls->theta[][i] */
    if (ls->size[i] < ls->threshold) Set_Column(ls->d,ls_theta,i,present); /* Set proposal = present */ 
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
  if (print == 1)
    {
    Rprintf("\nSample block parameters:");
    Rprintf("\n- log_ratio = %8.4f",log_ratio);  
    Rprintf("\n- decision = %i",accept);
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
                        double *input_proposal, double *input_present, int print, int n_between, double scale_factor)
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/
{
  int accept, i;
  double **cf, *present, *ergm_theta, log_present, log_proposal, log_ratio, **ls_theta, *mean, *proposal, *theta_present, *theta_proposal;
  /* Proposal:
  note 1: all ls->theta such that ls->size >= ls->threshold and all ergm->theta are updated by random walk Metropolis-Hastings algorithm
  note 2: ratio of proposal pdfs cancels under random walk Metropolis-Hastings algorithm */
  log_ratio = 0.0;
  /* Propose ergm->theta: */
  if (ergm->d1 > 0)
    {
    cf = Scale(ergm->d1,ergm->d1,prior->cf1,scale_factor); /* Rescale Cholesky factor of Gaussian prior */        
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
  if (ls_theta == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Independence, ls_theta\n\n"); exit(1); }
  for (i = 0; i < ls->d; i++)
    {
    ls_theta[i] = (double*) calloc(ls->number+1,sizeof(double));
    if (ls_theta[i] == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Independence, ls_theta[%i]\n\n",i); exit(1); }
    }
  present = (double*) calloc(ls->d,sizeof(double));  
  if (present == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Independence, present\n\n"); exit(1); }
  cf = Scale(ls->d,ls->d,prior->cf2,scale_factor); /* Rescale Cholesky factor of Gaussian prior */ 
  for (i = 0; i < ls->number; i++) 
    {
    Get_Column(ls->d,present,ls->theta,i); /* Set mean to ls->theta[][i] */
    if (ls->size[i] < ls->threshold) Set_Column(ls->d,ls_theta,i,present); /* Set proposal = present */ 
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
  if (print == 1)
    {
    Rprintf("\nSample parameters:");
    Rprintf("\n- log_ratio = %8.4f",log_ratio);  
    Rprintf("\n- decision = %i",accept);
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
                        int *newnetworkheads, int *newnetworktails, int n_between, double scale_factor, double *q_i)
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/
{
  int accept, i, k, n_input, proposal_n_edges, *proposal_heads, *proposal_tails, *ls_indicator, *ls_size, *sample_i, sample_size;
  double **cf, *present, *ergm_theta, *input_proposal, log_denominator, log_numerator, log_present, log_proposal, log_ratio, **ls_theta, *mean, *proposal, *theta_present, *theta_proposal, *statistic, sum, u;
  n_input = Number_Input(ergm->terms,input_present);
  input_proposal = (double*) calloc(n_input,sizeof(double));
  if (input_proposal == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Dependence, input_proposal\n\n"); exit(1); }
  ls_indicator = (int*) calloc(ls->n,sizeof(int));
  if (ls_indicator == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Dependence, ls_indicator\n\n"); exit(1); }
  ls_size = (int*) calloc(ls->number,sizeof(int));
  if (ls_size == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Dependence, ls_size\n\n"); exit(1); }
  ls_theta = (double**) calloc(ls->d,sizeof(double*));
  if (ls_theta == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Dependence, ls_theta\n\n"); exit(1); }
  for (i = 0; i < ls->d; i++)
    {
    ls_theta[i] = (double*) calloc(ls->number+1,sizeof(double));
    if (ls_theta[i] == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Dependence, ls_theta[%i]\n\n",i); exit(1); }
    }
  log_ratio = 0.0;
  u = unif_rand(); /* Decide whether indicators or parameters are to be updated */
  if (u < 0.5) /* Update indicators */
    {
    /* Propose indicators:
    - node i sampled with probability q_i[i]
    - for given node i with ls->indicator[i] = l, category k sampled with probability ls->p[k]
    - acceptance probability: proposal ratio x prior ratio x likelihood ratio, where
      * proposal ratio = (q_i[i] x ls->p[l]) / (q_i[i] * ls->p[k]) = ls->p[l] / ls->p[k]
      * prior ratio = ls->p[k] / ls->p[l]
      * proposal ratio x prior ratio = (ls->p[l] / ls->p[k]) x (ls->p[k] / ls->p[l]) = 1 and its logarithm 0
      * likelihood ratio to be computed below
    */
    sample_i = (int*) calloc(ls->n,sizeof(int));
    if (sample_i == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Dependence, sample_i\n\n"); exit(1); }
    sample_size = trunc(ls->n / 10.0);
    if (sample_size < 1) sample_size = ls->n;
    for (k = 0; k < sample_size; k++)
      {
      i = Sample_Discrete(q_i);
      sample_i[i] = 1;
      }
    for (i = 0; i < ls->n; i++) /* Given node, sample category */
      {
      if (sample_i[i] == 1) k = Sample_Discrete(ls->p); /* Given node, sample category */ 
      else k = ls->indicator[i];
      ls_indicator[i] = k;
      ls_size[k] = ls_size[k] + 1;
      }
    free(sample_i);
    }
  else /* Update parameters */
    {
    /* Propose ergm->theta: */
    if (ergm->d1 > 0)
      {
      cf = Scale(ergm->d1,ergm->d1,prior->cf1,scale_factor); /* Rescale Cholesky factor of Gaussian prior */        
      ergm_theta = Sample_MVN(ergm->d1,ergm->theta,cf); /* Random walk Metropolis-Hastings */
      log_proposal = MVN_PDF(ergm->d1,ergm_theta,prior->mean1,prior->precision1); /* Prior pdf: proposal */
      log_present = MVN_PDF(ergm->d1,ergm->theta,prior->mean1,prior->precision1); /* Prior pdf: present */
      log_ratio = log_ratio + (log_proposal - log_present);
      /*
      Rprintf("\n- log_ratio (parameters) = %8.4f",log_ratio);  
      */
      for (i = 0; i < ergm->d1; i++)
        { 
        free(cf[i]);
        }
      free(cf);
      }
    /* Propose ls->theta: */
    present = (double*) calloc(ls->d,sizeof(double));  
    if (present == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Dependence, present\n\n"); exit(1); }
    cf = Scale(ls->d,ls->d,prior->cf2,scale_factor); /* Rescale Cholesky factor of Gaussian prior */ 
    for (i = 0; i < ls->number; i++) 
      {
      Get_Column(ls->d,present,ls->theta,i); /* Set mean to ls->theta[][i] */
      if (ls->size[i] < ls->threshold) Set_Column(ls->d,ls_theta,i,present); /* Set proposal = present */ 
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
      /*
      Rprintf("\n- log_ratio (block parameters) = %8.4f",log_ratio);  
      */
      }
    free(present);
    for (i = 0; i < ls->d; i++)
      {
      free(cf[i]);
      }
    free(cf);
    }
  if (u < 0.5) /* Updated indicators but not parameters */
    {
    ergm_theta = (double*) calloc(ergm->d1,sizeof(double));
    if (ergm_theta == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Dependence, ergm_theta\n\n"); exit(1); }
    for (i = 0; i < ergm->d1; i++)
      {
      ergm_theta[i] = ergm->theta[i];
      }
    for (i = 0; i < ls->d; i++)
      { 
      for (k = 0; k < ls->number; k++)
        {
        ls_theta[i][k] = ls->theta[i][k];
        }
      }
    }
  else /* Updated parameters but not indicators */
    {
    for (i = 0; i < ls->n; i++) 
      {
      ls_indicator[i] = ls->indicator[i];
      }
    for (k = 0; k < ls->number; k++)
      {
      ls_size[k] = ls->size[k];
      } 
    }
  /* Propose auxiliary variable: */
  for (i = 0; i < n_input; i++) 
    {
    input_proposal[i] = input_present[i];
    }
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls_indicator,ls_theta,input_proposal); /* Set input given ls_theta */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input given ls->theta */
  theta_proposal = Get_Parameter(ergm->d,ergm->structural,ergm_theta); /* Set parameter; note: if ergm_d1 == 0, ergm_theta is not used */
  theta_present = Get_Parameter(ergm->d,ergm->structural,ergm->theta); /* Set parameter; note: if ergm->d1 == 0, ergm_theta is not used */
  sample_size = 1; /* One sample point is all that is required */
  MCMC_wrapper(heads,tails,dnedges,  /* Sample one graph from posterior predictive distribution given input and theta */
                  maxpossibleedges,
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
                  mheads,mtails,mdnedges);
  proposal_n_edges = newnetworkheads[0]; /* Number of simulated edges */
  proposal_heads = (int*) calloc(proposal_n_edges,sizeof(int)); /* Proposed heads for auxiliary variable */
  if (proposal_heads == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Dependence, proposal_heads\n\n"); exit(1); }
  proposal_tails = (int*) calloc(proposal_n_edges,sizeof(int)); /* Proposed tails for auxiliary variable */
  if (proposal_tails == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Dependence, proposal_tails\n\n"); exit(1); }
  for (i = 0; i < proposal_n_edges; i++)  
    {
    proposal_heads[i] = newnetworkheads[i+1]; /* Note: while heads corresponds to the list of observed heads, newnetworkheads contains the number of simulated edges as well as the list of simulated heads: to use auxiliary->heads here, one must not store the number of simulated edges */
    proposal_tails[i] = newnetworktails[i+1]; /* Note: while tails corresponds to the list of observed tails, newnetworktails contains the number of simulated edges as well as the list of simulated tails: to use auxiliary->tails here, one must not store the number of simulated edges */
    }
  /* Ratio of proposal pmfs of auxiliary graph under proposal / present */
  statistic = (double*) calloc(ergm->d,sizeof(double));
  if (statistic == NULL) { Rprintf("\n\ncalloc failed: Sample_Parameters_Dependence, statistic\n\n"); exit(1); }
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
  Rprintf("\n- log_ratio (auxiliary graph) = %8.4f",log_ratio);  
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
  Rprintf("\n- log_ratio (observed graph) = %8.4f",log_ratio);  
  */
  accept = MH_Decision(log_ratio);
  if (accept == 1) /* Proposal accepted */
    {
    Set_I_I(ls->n,ls->indicator,ls_indicator);
    Set_I_I(ls->number,ls->size,ls_size);
    if (ergm->d1 > 0) Set_D_D(ergm->d1,ergm->theta,ergm_theta);
    Set_DD_DD(ls->d,ls->number+1,ls->theta,ls_theta);
    }
  if (print == 1)
    {
    Rprintf("\nSample_Parameters_Dependence:");
    Rprintf("\n- log_ratio = %8.4f",log_ratio);  
    Rprintf("\n- decision = %i",accept);
    }
  if (ergm->d1 > 0) free(ergm_theta);
  free(ls_indicator);
  free(ls_size);
  free(proposal_heads);
  free(proposal_tails);
  free(statistic);
  free(theta_present);
  free(theta_proposal);
  for (i = 0; i < ls->d; i++)
    {
    free(ls_theta[i]);
    }
  free(ls_theta);
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
    if (ls->size[i] < ls->threshold) /* Structural parameter not showing up in ergm pmf */
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
  if (sample == NULL) { Rprintf("\n\ncalloc failed: Gibbs_Parameters_Means, sample\n\n"); exit(1); }
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

double* Gibbs_Parameters_Precisions(priorstructure *prior, latentstructure *ls)
/*
input: prior structure, latent structure
output: precisions of parameters
*/
{
  int i, k;
  double rate, shape, *sample, sum, d;
  sample = (double*) calloc(ls->d,sizeof(double)); 
  if (sample == NULL) { Rprintf("\n\ncalloc failed: Gibbs_Parameters_Precisions, sample\n\n"); exit(1); }
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

void Initial_State(int *parallel, double *alpha, int *indicator, priorstructure_ls *prior_ls, priorstructure *prior, latentstructure *ls, ergmstructure *ergm, double *theta, double scale_factor)
/* 
input: clustering parameter, priors, latent structure, ergm structure, user-specified initial value of non-structural parameters
*/
{  
  int i, k;
  double **cf, *sample, *shape1, *shape2, sum;
  if (*parallel == 1) ls->alpha = *alpha; /* Clustering parameter */
  else ls->alpha = rgamma(prior_ls->alpha_shape,1.0/prior_ls->alpha_rate); 
  shape1 = (double*) calloc(ls->number-1,sizeof(double)); /* Components 0..ls->number-2 suffice */
  if (shape1 == NULL) { Rprintf("\n\ncalloc failed: Initial_State, shape1\n\n"); exit(1); }
  shape2 = (double*) calloc(ls->number-1,sizeof(double)); /* Components 0..ls->number-2 suffice */
  if (shape2 == NULL) { Rprintf("\n\ncalloc failed: Initial_State, shape2\n\n"); exit(1); }
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
    k = Sample_Discrete(ls->p);
    ls->indicator[i] = k; 
    ls->size[k] = ls->size[k] + 1; /* ls-size was set to 0 by calloc */
    }
  if (ergm->d1 > 0)
    {
    cf = Scale(ergm->d1,ergm->d1,prior->cf1,scale_factor); /* Rescale Cholesky factor of Gaussian prior */        
    sample = Sample_MVN(ergm->d1,prior->mean1,cf);
    Set_D_D(ergm->d1,ergm->theta,sample);
    free(sample);
    for (i = 0; i < ergm->d1; i++)
      { 
      free(cf[i]);
      }
    free(cf);
    }
  cf = Scale(ls->d,ls->d,prior->cf2,scale_factor); /* Rescale Cholesky factor of Gaussian prior */ 
  for (i = 0; i < ls->number; i++) 
    {
    sample = Sample_MVN(ls->d,prior->mean2,cf); /* Random walk Metropolis-Hastings algorithm */
    Set_Column(ls->d,ls->theta,i,sample); /* Set ls_theta[][i] to proposal */
    free(sample);
    }
  for (i = 0; i < ls->d; i++)
    {
    free(cf[i]);
    }
  free(cf);
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
  if (p == NULL) { Rprintf("\n\ncalloc failed: Sample_CRP, p\n\n"); exit(1); }
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

int Sample_Graph_Edge_Independence(int *directed, latentstructure *ls, double *ln_p, int *heads, int *tails)
/*
input: latent structure; probability of edge between nodes i and j on log scale
output: graph sampled from PMF p and number of edges
*/
{
  int h, i, j, number_edges;
  double u;
  /* 
  Note 1: if i < j, edge (i, j) should be stored as (i, j)
  Note 2: i, j should be integers between 1 and n
  */  
  number_edges = 0;
  h = -1;
  for (i = 1; i < ls->n + 1; i++)
    {
    for (j = i + 1; j < ls->n + 1; j++)
      {
      u = unif_rand(); /* Sample uniform[0,1] */
      h = h + 1;
      if (ln(u) < ln_p[h]) 
        {
        number_edges = number_edges + 1;
        heads[number_edges] = i; /* number_edges is between 1..degrees of freedom; heads[0] is supposed to contain the number of edges */
        tails[number_edges] = j; /* number_edges is between 1..degrees of freedom; tails[0] is supposed to contain the number of edges */
        }
      }
    }    
  heads[0] = number_edges; /* heads[0] is supposed to contain the number of edges */
  tails[0] = number_edges; /* tails[0] is supposed to contain the number of edges */
  return number_edges;
}

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
             int *max_iterations, int *n_between_block_parameters, int *output, double *mcmc, int *sample_heads, int *sample_tails, int *call_RNGstate, int *hyperprior)
/*
input: R input
output: simulated graph
*/
{
  int null = 0;
  int coordinate, *degree, *degree_freq, dim, dim1, dim2, edges, element, h, i, j, hyper_prior, dyad_dependence, *n_edges, *pseudo_indicator, iteration, k, max_iteration, *mdnedges, *mheads, *mtails, n, *newnetworkheads, *newnetworktails, number, print, threshold, terms, *verbose;
  double *draw, *p, **parameter, *pp, progress, *shape1, *shape2, sum;	
  priorstructure_ls *prior_ls;
  latentstructure *ls;
  priorstructure *prior;
  ergmstructure *ergm;
  /**************/
  /* Initialize */
  /**************/
  print = *v; /* Console: no print; 0: short print; 1: long print */
  verbose = &null;
  epsilon = DBL_EPSILON;
  maximum = DBL_MAX;
  if (print >= 2)
    {
    Rprintf("\nMachine precision:");
    Rprintf("\n- epsilon = %e",epsilon);
    Rprintf("\n- maximum = %e",maximum);
    /*
    Rprintf("\n- ln(epsilon) = %e",ln(epsilon));
    Rprintf("\n- ln(maximum) = %e",ln(maximum));
    Rprintf("\n- exp(-maximum) = %e",e(-maximum));
    Rprintf("\n- exp(-epsilon)= %e",e(-epsilon));
    Rprintf("\n- exp(+epsilon) = %e",e(+epsilon));
    Rprintf("\n- exp(+maximum) = %e",e(+maximum));
    */
    Rprintf("\n");
    }
  terms = (int)*nterms; /* Number of ergm terms */
  dim = (int)*d;
  dim1 = (int)*d1;
  dim2 = (int)*d2;
  threshold = (int)*min_size; /* Minimum size of category so that structural parameters show up in ergm pmf */
  n = (int)*dn; /* Number of nodes */
  number = (int)*max_number; /* Number of categories */
  max_iteration = *max_iterations; /* Number of draws from posterior */
  ergm = Initialize_Ergm(terms,hierarchical,dim,dim1,dim2,structural); /* Ergm structure and non-structural parameters */
  prior = Initialize_Prior(ergm->d1,ergm->d2,m2_mean,m2_precision,*p2_shape,*p2_rate,m1,m2,b,cf1,cf2,p1,p2); /* Prior: non-structural, structural parameters */
  ls = Initialize_Latentstructure(number,n,threshold,ergm->d2); /* Latent structure and structural parameters */
  prior_ls = Initialize_Prior_ls(*alpha_shape,*alpha_rate); /* Prior: clustering parameter */
  mdnedges = &null;
  mheads = NULL;
  mtails = NULL;
  shape1 = (double*) calloc(ls->number-1,sizeof(double)); /* Components 0..ls->number-2 suffice */
  if (shape1 == NULL) { Rprintf("\n\ncalloc failed: Simulation, shape1\n\n"); exit(1); }
  shape2 = (double*) calloc(ls->number-1,sizeof(double)); /* Components 0..ls->number-2 suffice */
  if (shape2 == NULL) { Rprintf("\n\ncalloc failed: Simulation, shape2\n\n"); exit(1); }
  p = (double*) calloc(*maxpossibleedges,sizeof(double));
  if (p == NULL) { Rprintf("\n\ncalloc failed: Simulation, p\n\n"); exit(1); }
  parameter = (double**) calloc(ls->d,sizeof(double*));
  if (parameter == NULL) { Rprintf("\n\ncalloc failed: Simulation, parameter\n\n"); exit(1); }
  for (i = 0; i < ls->d; i++)
    {
    parameter[i] = (double*) calloc(ls->number+1,sizeof(double));
    if (parameter[i] == NULL) { Rprintf("\n\ncalloc failed: Simulation, parameter[%i]\n\n"); exit(1); }
    }
  pseudo_indicator = (int*) calloc(ls->n,sizeof(int));
  if (pseudo_indicator == NULL) { Rprintf("\n\ncalloc failed: Simulation, pseudo_indicator\n\n"); exit(1); }
  pp = (double*) calloc(ergm->d,sizeof(double));
  if (pp == NULL) { Rprintf("\n\ncalloc failed: Simulation, pp\n\n"); exit(1); }
  newnetworkheads = (int*) calloc(*maxpossibleedges,sizeof(int));
  if (newnetworkheads == NULL) { Rprintf("\n\ncalloc failed: Simulation, newnetworkheads\n\n"); exit(1); }
  newnetworktails = (int*) calloc(*maxpossibleedges,sizeof(int));
  if (newnetworktails == NULL) { Rprintf("\n\ncalloc failed: Simulation, newnetworktails\n\n"); exit(1); }
  hyper_prior = (int)*hyperprior; /* Means and precisions of Gaussian baseline distribution have non-degenerate prior */
  dyad_dependence = (int)*dyaddependence; /* Conditional PMF of graph given latent structure: dyad-dependent or not */
  if (print >= 1)
    {
    if (hyper_prior == 0) Rprintf("\nPosterior prediction.");
    else Rprintf("\nSimulation.");
    if (dyad_dependence == 0) Rprintf("\nConditional dyad-independence model.\n");
    else Rprintf("\nDyad-dependence model.\n");
    }
  /****************/
  /* Sample graph */
  /****************/
  if (*call_RNGstate == 1) GetRNGstate();
  ls->alpha = *alpha; /* Clustering parameter */
  for (i = 0; i < (ls->number - 1); i++)
    {
    shape1[i] = 1.0; /* First shape of Beta distribution */
    shape2[i] = ls->alpha; /* Second shape of Beta distribution */
    }
  if (hyper_prior == 0) /* Posterior prediction */
    {
    k = -1;
    for (i = 0; i < ergm->d1; i++)
      {
      k = k + 1;
      ergm->theta[i] = eta[k];
      }
    for (i = 0; i < ls->d; i++)
      {
      for (j = 0; j < ls->number; j++)
        {
        k = k + 1;
        ls->theta[i][j] = eta[k];
        }
      }
    Set_I_I(ls->n,ls->indicator,indicator);
    }
  coordinate = -1;
  element = -1;
  for (iteration = 0; iteration < max_iteration; iteration++)
    {
    progress = (iteration * 100.0) / max_iteration;
    if (print == 1) Rprintf("\nProgress: %5.2f%%",progress);
    if (hyper_prior == 0) /* Posterior prediction */
      {
      for (k = 0; k < ls->number; k++)
        {
        ls->size[k] = 0;
        }
      for (i = 0; i < ls->n; i++)
        {
        k = ls->indicator[i];
        ls->size[k] = ls->size[k] + 1;
        }
      }
    else /* Simulation */
      {
      ls->p = Stick_Breaking(shape1,shape2,ls); /* Construct category probability vector by stick-breaking */
      if (ls->number < ls->n) /* Sample partition of set of nodes (truncated case) */
        {
        for (k = 0; k < ls->number; k++)
          {
          ls->size[k] = 0;
          }
        for (i = 0; i < ls->n; i++)
          {
          k = Sample_Discrete(ls->p);
          ls->indicator[i] = k;
          ls->size[k] = ls->size[k] + 1;
          }
        }
      else Sample_CRP(ls); /* Sample partition of set of nodes (untruncated case) */
      if (ergm->d1 > 0) 
        {
        draw = Sample_MVN(ergm->d1,prior->mean1,prior->cf1);
        Set_D_D(ergm->d1,ergm->theta,draw);
        free(draw);
        }
      for (i = 0; i < ls->number; i++) 
        {
        draw = Sample_MVN(ls->d,prior->mean2,prior->cf2); /* Random walk Metropolis-Hastings algorithm */
        Set_Column(ls->d,ls->theta,i,draw); /* Set ls_theta[][i] to proposal */
        free(draw);
        }
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
      Sample_Graph_Edge_Independence(directed,ls,p,newnetworkheads,newnetworktails);
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
      if (mheads == NULL) { Rprintf("\n\ncalloc failed: Simulation, mheads\n\n"); exit(1); }
      mtails = (int*) calloc(*n_edges,sizeof(int));
      if (mtails == NULL) { Rprintf("\n\ncalloc failed: Simulation, mtails\n\n"); exit(1); }
      for (i = 0; i < *n_edges; i++) /* Since first element of newnetworkheads and newnetworktails is number of edges, heads and tails must be extracted */
        {
        mheads[i] = newnetworkheads[i+1];
        mtails[i] = newnetworktails[i+1];
        }
      network_stats_wrapper(mheads,mtails,n_edges,dn,directed,bipartite,nterms,funnames,sonames,inputs,pp); /* Compute non-structural function of graph */
      if (print == 1)
        {
        degree = Degree_Sequence(ls->n,*directed,*n_edges,mheads,mtails);
        degree_freq = Degree_Freq(ls->n,degree);
        }
      free(mheads);
      free(mtails);
      }
    else Sample_Graph(ls->number,ls->n,ls->d,ergm->terms,ergm->hierarchical,ergm->d,pp, /* Sample graph */
                         heads,tails,dnedges,maxpossibleedges,dn,directed,bipartite,nterms,funnames,
                         sonames,MHproposaltype,MHproposalpackage,inputs,theta,samplesize,
                         sample,burnin,interval,newnetworkheads,newnetworktails,verbose,
                         attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,
                         maxedges,mheads,mtails,mdnedges);
    if (ergm->d1 > 0)
      {
      if (print == 1) Rprintf("\nparameters:");
      for (i = 0; i < ergm->d1; i++) /* Non-structural parameters */
        {
        if (print == 1) Rprintf(" %8.4f",ergm->theta[i]);
        coordinate = coordinate + 1;
        mcmc[coordinate] = ergm->theta[i];
        }  
      }
    if ((hyper_prior == 1) && (print == 1)) Rprintf("\nmeans of block parameters:");
    for (i = 0; i < ls->d; i++) /* Structural parameters */
      {
      if ((hyper_prior == 1) && (print == 1)) Rprintf(" %8.4f",prior->mean2[i]);
      coordinate = coordinate + 1;	
      mcmc[coordinate] = prior->mean2[i];
      }
    if ((hyper_prior == 1) && (print == 1)) Rprintf("\nprecisions of block parameters:");
    for (i = 0; i < ls->d; i++) /* Structural parameters */
      {
      if ((hyper_prior == 1) && (print == 1)) Rprintf(" %8.4f",prior->precision2[i][i]);
      coordinate = coordinate + 1;	
      mcmc[coordinate] = prior->precision2[i][i];
      }
    if (print == 1) Rprintf("\nblock parameters:\n");
    for (h = 0; h < ls->d; h++) /* Structural parameters */
      {
      for (i = 0; i < ls->number; i++) 
        {
        if (print == 1) Rprintf(" %8.4f",ls->theta[h][i]);
        coordinate = coordinate + 1;	
        mcmc[coordinate] = ls->theta[h][i];
        }
     if ((print == 1) && (ls->theta[h][ls->number] != 0)) Rprintf(" %8.4f",ls->theta[h][ls->number]); /* Second condition ensures that between-category parameters are not written to screen when model without between-category parameters */
      coordinate = coordinate + 1;	
      mcmc[coordinate] = ls->theta[h][ls->number];
      if (print == 1) Rprintf("\n");
      }
    if (print == 1) Rprintf("block indicators:");
    for (i = 0; i < ls->n; i++) /* Category indicators */
      {
      if (print == 1) 
        {
        if (*dyaddependence == 0) Rprintf(" (%i:%i)%i",i+1,degree[i],ls->indicator[i]+1);
        else Rprintf(" (%i)%i",i+1,ls->indicator[i]+1);
        }
      coordinate = coordinate + 1;
      mcmc[coordinate] = ls->indicator[i];
      }
    if (print == 1) Rprintf("\nblock sizes:");  
    for (i = 0; i < ls->number; i++)
      {
      if (print == 1) Rprintf(" %3i",ls->size[i]);
      coordinate = coordinate + 1;
      mcmc[coordinate] = ls->size[i];
      } 
    if ((hyper_prior == 1) && (print == 1)) Rprintf("\nblock probabilities:");
    for (i = 0; i < ls->number; i++) /* Category probability vector */
      {
      if ((hyper_prior == 1) && (print == 1)) Rprintf(" %6.4f",ls->p[i]);
      coordinate = coordinate + 1;
      mcmc[coordinate] = ls->p[i];
      }
    if (print == 1) Rprintf("\nblock probabilities prior parameter: %6.4f",ls->alpha); /* Clustering parameter */
    coordinate = coordinate + 1;
    mcmc[coordinate] = ls->alpha;
    if ((*dyaddependence == 0) && (*directed == 0) && (print == 1)) 
      {
      Rprintf("\ndegree distribution:\n");
      i = 0;
      for (k = 0; k < (ls->n - 1); k++)
        {
        if (degree_freq[k] > 0) 
          {
          i = i + 1;
          Rprintf("%4i",k);
          }
        } 
      Rprintf("\n");
      for (k = 0; k < (ls->n - 1); k++)
        {
        if (degree_freq[k] > 0) Rprintf("%4i",degree_freq[k]);
        } 
      Rprintf("\n%i of %i possible values",i,ls->n);
      free(degree);
      free(degree_freq);
      }
    if (print == 1) Rprintf("\nstatistics:");
    for (i = 0; i < ergm->d; i++) /* Statistics */
      {
      if (print == 1) Rprintf(" %6.0f",pp[i]);
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
    if (print == 1) Rprintf("\n");
    }
  if (*call_RNGstate == 1) PutRNGstate();
  /************/
  /* Finalize */
  /************/
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
}

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
             int *max_iterations, int *n_between_block_parameters, int *output, double *mcmc, double *scalefactor, double *mh_accept, double *q_i, int *call_RNGstate, int *parallel, int *hyperprior)
/*
input: R input
output: MCMC sample of unknowns from posterior
*/
{
  int null = 0;
  int batch, n_batches, batch_size, coordinate, console, *degree, *degree_freq, dyad_dependence, dim, dim1, dim2, h, i, j, k, hyper_prior, *mdnedges, *mheads, *mtails, n_input, iteration, max_iteration, n, n_between, number, print, store, threshold, terms, *verbose;
  double ls_alpha, accept, *block_degree_freq, local_mh_accept, *ls_p, *pp, *prior_mean2, *prior_precision2, progress, rate, shape, scale_factor, u;	
  priorstructure_ls *prior_ls;
  latentstructure *ls;
  priorstructure *prior;
  ergmstructure *ergm;
  /**************/
  /* Initialize */
  /**************/
  /*
  SEXP NEWRANDOMSEED;
  SEXP rs;
  PROTECT(rs = install(".Random.seed"));
  */
  console = *v; /* Console: -1: no print; 0: short print; 1: long print */
  verbose = &null;
  epsilon = DBL_EPSILON;
  maximum = DBL_MAX;
  if (console >= 2)
    {
    Rprintf("\nMachine precision:");
    Rprintf("\n- epsilon = %e",epsilon);
    Rprintf("\n- maximum = %e",maximum);
    /*
    Rprintf("\n- ln(epsilon) = %e",ln(epsilon));
    Rprintf("\n- ln(maximum) = %e",ln(maximum));
    Rprintf("\n- exp(-maximum) = %e",e(-maximum));
    Rprintf("\n- exp(-epsilon)= %e",e(-epsilon));
    Rprintf("\n- exp(+epsilon) = %e",e(+epsilon));
    Rprintf("\n- exp(+maximum) = %e",e(+maximum));
    */
    Rprintf("\n");
    }
  terms = (int)*nterms; /* Number of ergm terms */
  dim = (int)*d;
  dim1 = (int)*d1;
  dim2 = (int)*d2;
  threshold = (int)*min_size; /* Minimum size of category so that structural parameters show up in ergm pmf */
  number = (int)*max_number; /* Number of categories */
  n = (int)*dn; /* Number of nodes */
  max_iteration = *max_iterations; /* Number of draws from posterior */
  if (max_iteration <= 12000) n_batches = max_iteration;
  else n_batches = 12000;
  if (max_iteration == n_batches) batch_size = 1; 
  else batch_size = trunc(max_iteration / n_batches);
  n_between = (int)*n_between_block_parameters; /* Number of between-category parameters */
  mdnedges = &null;
  scale_factor = *scalefactor; /* Metropolis-Hasting algorithm: scale factor */
  ergm = Initialize_Ergm(terms,hierarchical,dim,dim1,dim2,structural); /* Ergm structure and non-structural parameters */
  prior = Initialize_Prior(ergm->d1,ergm->d2,m2_mean,m2_precision,*p2_shape,*p2_rate,m1,m2,b,cf1,cf2,p1,p2); /* Prior: non-structural, structural parameters */
  ls = Initialize_Latentstructure(number,n,threshold,ergm->d2); /* Latent structure and structural parameters */
  prior_ls = Initialize_Prior_ls(*alpha_shape,*alpha_rate); /* Prior: clustering parameter */
  hyper_prior = (int)*hyperprior; /* Means and precisions of Gaussian baseline distribution have non-degenerate prior */
  dyad_dependence = (int)*dyaddependence; /* Conditional PMF of graph given latent structure: dyad-dependent or not */
  if (console >= 1)
    {
    if (hyper_prior == 0) Rprintf("\nDirichlet process prior without hyper prior.");
    else Rprintf("\nDirichlet process prior: Gaussian baseline distribution with hyper prior.");
    if (dyad_dependence == 0) Rprintf("\nConditional dyad-independence model.\n");
    else Rprintf("\nDyad-dependence model.\n");
    }
  mheads = NULL;
  mtails = NULL;
  pp = (double*) calloc(ergm->d,sizeof(double));
  if (pp == NULL) { Rprintf("\n\ncalloc failed: Inference, pp\n\n"); exit(1); }
  /*************************/
  /* MCMC sample posterior */
  /*************************/
  if (*call_RNGstate == 1) GetRNGstate();
  if (console >= 1)
    {
    Rprintf("\nNumber of draws from posterior: %i",n_batches * batch_size);
    Rprintf("\nNumber of batches: %i",n_batches);
    Rprintf("\nSize of batches: %i",batch_size);
    Rprintf("\n");
    }
  if ((dyad_dependence == 0) && (console == 1)) 
    {
    degree = Degree_Sequence(ls->n,*directed,*dnedges,heads,tails);
    degree_freq = Degree_Freq(ls->n,degree);
    Rprintf("\nObserved degree distribution:\n");
    i = 0;
    for (k = 0; k < ls->n; k++)
      {
      if (degree_freq[k] > 0) 
        {
        i = i + 1;
        Rprintf("%4i",k);
        }
      } 
    Rprintf("\n");
    for (k = 0; k < ls->n; k++)
      {
      if (degree_freq[k] > 0) Rprintf("%4i",degree_freq[k]);
      } 
    Rprintf("\n%i of %i possible values observed\n",i,ls->n);
    }
  Initial_State(parallel,alpha,indicator,prior_ls,prior,ls,ergm,theta,scale_factor);
  local_mh_accept = 0.0;
  coordinate = -1;
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
      /*
      NEWRANDOMSEED = findVar(rs, R_GlobalEnv); 
      Rprintf("%d %d %d %d\n",INTEGER(NEWRANDOMSEED)[0],
                              INTEGER(NEWRANDOMSEED)[1],
 			      INTEGER(NEWRANDOMSEED)[2],
 			      INTEGER(NEWRANDOMSEED)[3]);
      */
      if (hyper_prior == 1) /* Hyper prior: mean and precisions of Gaussian baseline distribution have non-degenerate prior */
        {
        prior_mean2 = Gibbs_Parameters_Means(prior,ls); /* Sample means of parameters */
        Set_D_D(ls->d,prior->mean2,prior_mean2);
        prior_precision2 = Gibbs_Parameters_Precisions(prior,ls); /* Sample precisions of parameters */
        for (i = 0; i < ls->d; i++)
          {
          prior->precision2[i][i] = prior_precision2[i];
          }
        free(prior_mean2);
        free(prior_precision2); 
        } 
      if (dyad_dependence == 0) /* MCMC exploiting dyad-independence conditional on latent structure */
        {
        if (ergm->d1 > 0)
          {
          u = unif_rand();
          if (u < 0.00) accept = Sample_Ergm_Theta_Independence(ergm,ls,prior,heads,tails,dnedges,dn,directed,bipartite, 
                                   nterms,funnames,sonames,inputs,print,n_between,scale_factor);
          else if (u < 0.00) Sample_Ls_Theta_Independence(ergm,ls,prior,heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,
                                   inputs,inputs_h,print,n_between,scale_factor);
          else accept = Sample_Parameters_Independence(ergm,ls,prior, /* M-H exploiting dyad-independence conditional on latent structure */
                                   heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,inputs,inputs_h,print,n_between,scale_factor);
          }
        else accept = Sample_Parameters_Independence(ergm,ls,prior, /* M-H exploiting dyad-independence conditional on latent structure */
                                   heads,tails,dnedges,dn,directed,bipartite,nterms,funnames,sonames,inputs,inputs_h,print,n_between,scale_factor);
        Gibbs_Indicators_Independence(ls,ergm,heads,tails,inputs_h,dnedges,dn,directed,bipartite,nterms,funnames,sonames,q_i); 
        }
      else accept = Sample_Parameters_Dependence(ergm,ls,prior, /* Auxiliary-variable M-H */
                           heads,tails,dnedges,maxpossibleedges,dn,directed,bipartite,nterms,funnames,
                           sonames,MHproposaltype,MHproposalpackage,sample,burnin,interval, 
                           verbose,attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,
                           maxedges,mheads,mtails,mdnedges,inputs,print,newnetworkheads,newnetworktails,n_between,scale_factor,q_i);
      local_mh_accept = local_mh_accept + accept;
      Gibbs_Parameters(ergm,ls,prior); /* Structural parameters not showing up in ergm pmf */
      ls_p = Sample_P(ls); /* Category probability vector */ 
      Set_D_D(ls->number,ls->p,ls_p);
      free(ls_p);
      ls_alpha = Sample_Alpha(prior_ls,ls); /* Clustering parameter */
      ls->alpha = ls_alpha;
      if (store == 1) 
        {
        /* Posterior prediction when full output is desired */
        if (*output == 1)
          {
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
          Sample_Graph(ls->number,ls->n,ls->d,ergm->terms,ergm->hierarchical,ergm->d,pp, /* Posterior prediction of graph */
                           heads,tails,dnedges,maxpossibleedges,dn,directed,bipartite,nterms,funnames,
                           sonames,MHproposaltype,MHproposalpackage,inputs,theta,samplesize,
                           sample,burnin,interval,newnetworkheads,newnetworktails,verbose,
                           attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,
                           maxedges,mheads,mtails,mdnedges);
          }
        /* Store and output MCMC sample: */ 
        if (ergm->d1 > 0)
          {
          if (print == 1) Rprintf("\nparameters:");
          for (i = 0; i < ergm->d1; i++) /* Non-structural parameters */
            {
            if (print == 1) Rprintf(" %8.4f",ergm->theta[i]);
            coordinate = coordinate + 1;
            mcmc[coordinate] = ergm->theta[i];
            }  
          }
        if (print == 1) Rprintf("\nmeans of block parameters:");
        for (i = 0; i < ls->d; i++) /* Structural parameters */
          {
          if (print == 1) Rprintf(" %8.4f",prior->mean2[i]);
          coordinate = coordinate + 1;	
          mcmc[coordinate] = prior->mean2[i];
          }
        if (print == 1) Rprintf("\nprecisions of block parameters:");
        for (i = 0; i < ls->d; i++) /* Structural parameters */
          {
          if (print == 1) Rprintf(" %8.4f",prior->precision2[i][i]);
          coordinate = coordinate + 1;	
          mcmc[coordinate] = prior->precision2[i][i];
          }
        if (print == 1) Rprintf("\nblock parameters:\n");
        for (h = 0; h < ls->d; h++) /* Structural parameters */
          {
          for (i = 0; i < ls->number; i++) 
            {
            if (print == 1) Rprintf(" %8.4f",ls->theta[h][i]);
            coordinate = coordinate + 1;	
            mcmc[coordinate] = ls->theta[h][i];
            }
          if ((print == 1) && (ls->theta[h][ls->number] != 0)) Rprintf(" %8.4f",ls->theta[h][ls->number]); /* Second condition ensures that between-category parameters are not written to screen when model without between-category parameters */
          coordinate = coordinate + 1;	
          mcmc[coordinate] = 0.0;
          if (print == 1) Rprintf("\n");
          }
        if (print == 1) Rprintf("block indicators:");
        for (i = 0; i < ls->n; i++) /* Category indicators */
          {
          if (print == 1) 
            {
            if (dyad_dependence == 0) Rprintf(" (%i:%i)%i",i+1,degree[i],ls->indicator[i]+1);
            else Rprintf(" (%i)%i",i+1,ls->indicator[i]+1);
            }
          coordinate = coordinate + 1;
          mcmc[coordinate] = ls->indicator[i];
          }
        if (print == 1) Rprintf("\nblock sizes:");  
        for (i = 0; i < ls->number; i++)
          {
          if (print == 1) Rprintf(" %3i",ls->size[i]);
          coordinate = coordinate + 1;
          mcmc[coordinate] = ls->size[i];
          } 
        if (print == 1) Rprintf("\nblock probabilities:");
        for (i = 0; i < ls->number; i++) /* Category probability vector */
          {
          if (print == 1) Rprintf(" %6.4f",ls->p[i]);
          coordinate = coordinate + 1;
          mcmc[coordinate] = ls->p[i];
          }
        if (print == 1) Rprintf("\nblock probabilities prior parameter: %6.4f",ls->alpha); /* Clustering parameter */
        coordinate = coordinate + 1;
        mcmc[coordinate] = ls->alpha;
        if ((dyad_dependence == 0) && (print == 1))
          {
          block_degree_freq = Block_Degree_Freq(ls->n,degree,ls->number,ls->size,ls->indicator);
          Rprintf("\nobserved degree by block:");
          for (i = 0; i < ls->number; i++)
            {
            Rprintf(" %6.4f",block_degree_freq[i]);
            }
          free(block_degree_freq); 
          }
        if ((*output == 1) && (print == 1)) Rprintf("\nposterior prediction of statistics:");
        for (i = 0; i < ergm->d; i++) /* Posterior prediction of statistics */
          {
          if ((*output == 1) && (print == 1)) Rprintf(" %6.0f",pp[i]);
          coordinate = coordinate + 1;
          mcmc[coordinate] = pp[i];
          }
        if (print == 1) Rprintf("\n");
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
  if (*call_RNGstate == 1) PutRNGstate();
  /************/
  /* Finalize */
  /************/
  local_mh_accept = local_mh_accept / max_iteration;
  *mh_accept = local_mh_accept;
  if (console >= 1)
    {
    Rprintf("\n");
    Rprintf("\nNumber of draws from posterior: %i",n_batches * batch_size);
    Rprintf("\nThinning: every %i-th draw recorded",batch_size);
    Rprintf("\nAcceptance rate of Metropolis-Hastings algorithm: %6.4f",local_mh_accept);
    }
  free(pp);
  free(mheads);
  free(mtails);
  if ((dyad_dependence == 0) && (print == 1))     
    {
    free(degree);
    free(degree_freq);
    }
  Finalize_Ergm(ergm);
  Finalize_Latentstructure(ls,dim2);
  Finalize_Prior_ls(prior_ls);
  Finalize_Priorstructure(prior,dim1,dim2);
  /*
  UNPROTECT(1);
  */
}

