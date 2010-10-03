#include "h_ergm.h"

void Sample_Alpha(priorstructure_ls *prior_ls, latentstructure *ls)
/*
input: prior, latent structure
output: clustering parameter
*/
{
  double shape, rate;
  /* GetRNGstate(); */
  shape = prior_ls->alpha_shape + (ls->number - 1.0); /* Shape of full conditional Gamma */
  rate = prior_ls->alpha_rate - ln(ls->p[ls->number-1]); /* Rate (inverse scale) of full conditional Gamma */
  ls->alpha = rgamma(shape,1.0/rate); 
  /* PutRNGstate(); */
}

void Stick_Breaking(double *shape1, double *shape2, latentstructure *ls)
/*
input: shape parameters of Beta distribution, latent structure
output: category probability vector
*/
{
  int i;
  double *b, c;
  /* GetRNGstate(); */
  b = D(ls->number);
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
  ls->p[0] = b[0];
  /*
  Rprintf("\nb[%i] = %f, p[%i] = %f",0,b[0],0,ls->p[0]);
  */
  c = 1.0;
  for (i = 1; i < ls->number; i++)
    {
    c = c * (1.0 - b[i-1]);
    ls->p[i] = b[i] * c; 
    /*
    Rprintf("\nb[%i] = %f c = %f p[%i] = %f",i,b[i],c,i,ls->p[i]);
    */
    }
  /*
  c = 0;
  for (i = 0; i < ls->number; i++)
    {
    c = c + ls->p[i];
    }
  Rprintf("\nlength of ls->p = %f",c);
  */
  /* PutRNGstate(); */
}

void Sample_P(latentstructure *ls)
/*
input: latent structure
output: category probability vector
*/
{
  int i, rest;
  double *shape1, *shape2;
  shape1 = D(ls->number-1); /* Element 0..ls->number-2 required */
  shape2 = D(ls->number-1); /* Element 0..ls->number-2 required */
  rest = ls->n;
  for (i = 0; i < (ls->number - 1); i++)
    {
    rest = rest - ls->size[i]; /* Number of nodes in category i + 1, ..., ls->number */
    shape1[i] = 1.0 + ls->size[i]; /* First shape parameter of Beta distribution */
    shape2[i] = ls->alpha + rest; /* Second shape parameter of Beta distribution */
    }
  Stick_Breaking(shape1,shape2,ls); /* Construct category probability vector by stick-breaking */
}

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
                       double *input_present, double *theta, int print)
/*
input: latent structure
output: category indicators
*/
{
  int i, j, k, present_i, present_j, present_k, sample;
  double u, *p;
  /* GetRNGstate(); */
  u = unif_rand();
  /* PutRNGstate(); */
  if (u < 1.0) sample = 1; /* Sample blocks of one node indicator */
  else if (u < 1.0) sample = 2; /* Sample block of two node indicators */
  else sample = 3; /* Sample block of three node indicators */
  p = D(ls->n);
  for (i = 0; i < ls->n; i++)
    {
    p[i] = 1.0 / ls->n;
    }
  i = Sample_Discrete(p); /* Sample node i */
  do /* Sample node j */
    {
    j = Sample_Discrete(p);
    }
  while (i == j);
  do /* Sample node k */
    {
    k = Sample_Discrete(p);
    }
  while ((i == k) || (j == k));
  present_i = ls->indicator[i];
  present_j = ls->indicator[j];
  present_k = ls->indicator[k];
  if (sample == 1) /* Sample blocks of one node indicator */
    {
    Gibbs_Indicators_1(i,j,k,ls,ergm,heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite, 
                       nterms,funnames,sonames,MHproposaltype,MHproposalpackage,samplesize, 
                       burnin,interval,newnetworkheads,newnetworktails,verbose,attribs,
                       maxout,maxin,minout,minin,condAllDegExact,attriblength,maxedges,
                       mheads,mtails,mdnedges,input_present,theta); /* Sample indicator of node i */
    }
  else if (sample == 2) /* Sample block of two node indicators */
    {
    Gibbs_Indicators_2(i,j,ls,ergm,heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite, 
                       nterms,funnames,sonames,MHproposaltype,MHproposalpackage,samplesize, 
                       burnin,interval,newnetworkheads,newnetworktails,verbose,attribs,
                       maxout,maxin,minout,minin,condAllDegExact,attriblength,maxedges,
                       mheads,mtails,mdnedges,input_present,theta); /* Sample block of indicators of i and j */
    }
  else /* Sample block of three node indicators */
    {
    Gibbs_Indicators_3(i,j,k,ls,ergm,heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite, 
                       nterms,funnames,sonames,MHproposaltype,MHproposalpackage,samplesize, 
                       burnin,interval,newnetworkheads,newnetworktails,verbose,attribs,
                       maxout,maxin,minout,minin,condAllDegExact,attriblength,maxedges,
                       mheads,mtails,mdnedges,input_present,theta); /* Sample block of indicators of i, j, and k */
    }  
  if (print == 1)
    {
    Rprintf("\nSample_Indicators: Gibbs sampling block of %i indicators:",sample);
    Rprintf("\nindicator[%i] %i > %i",i,present_i,ls->indicator[i]);
    Rprintf("\nindicator[%i] %i > %i",j,present_j,ls->indicator[j]);
    if (sample != 2) Rprintf("\nindicator[%i] %i > %i",k,present_k,ls->indicator[k]);
    }
}

void Set_Input_i(int i, double *input_present, int *index, int n_input, double **input, ergmstructure *ergm, latentstructure *ls)
/*
input: node, input (vector), starting index of row of input (matrix), number of columns of input (matrix), input (matrix), ergm structure, latent structure
output: input where indicator of node is set to all possible values 
*/
{
  int k, l, m, present;
  present = ls->indicator[i]; /* Store indicator */
  m = *index; /* Get row index of input */
  for (k = 0; k < ls->number; k++) /* Set input to input where indicator of node i is set to category k */
    {
    m = m + 1;
    ls->indicator[i] = k; /* Set indicator */
    Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input_present given latent structure  */
    for (l = 0; l < n_input; l++) 
      {
      input[m][l] = input_present[l]; /* Store input_present in row k of input */
      }
    }
  *index = m; /* Set row index of input */
  ls->indicator[i] = present; /* Reset indicator */
}

void Gibbs_Indicator_i(int i, double *p, latentstructure *ls)
/* 
input: node, full conditional, latent structure
output: category indicator
*/
{
  int k, l, present;
  double sum;  
  sum = 0;
  for (k = 0; k < ls->number; k++) /* Full conditional may be unnormalized: compute sum... */
    {
    sum = sum + p[k];
    }
  for (k = 0; k < ls->number; k++) /* Full conditional may be unnormalized: ...normalize */
    {
    p[k] = p[k] / sum;
    }
  present = ls->indicator[i]; /* Old category */
  k = Sample_Discrete(p); /* Sample full conditional */
  ls->indicator[i] = k; /* New category */
  ls->size[present] = ls->size[present] - 1; /* Decrement size of old category */
  ls->size[k] = ls->size[k] + 1; /* Increment size of new category */
}

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
                       double *input_present, double *theta)
/*
input: node, latent structure
output: category indicators
*/
{
  int generating, k, l, m, *m_p, n_input, number;
  double **input, *p, *p_i;
  number = 1 + (3 * ls->number);  /* Dimension: 1 generating indicator (present indicators) + 3 x ls->number proposed indicators (ls->number proposed indicators for 3 nodes) */
  n_input = Number_Input(ergm->terms,input_present); /* Number of input parameters */
  input = DD(number,n_input);
  p = D(number);
  p_i = D(ls->number);
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input_present given latent structure  */
  m = 0;
  for (l = 0; l < n_input; l++) 
    {
    input[m][l] = input_present[l]; /* Store input_present in row &m = 0 of input */
    }
  m_p = &m; /* Pointer m_p points to index of row m; to be passed to function Set_Input_i; function Set_Input_i is to manipulate m */ 
  Set_Input_i(h,input_present,m_p,n_input,input,ergm,ls); /* Set rows of input corresponding to changing indicators of node h */
  Set_Input_i(i,input_present,m_p,n_input,input,ergm,ls); /* Set rows of input corresponding to changing indicators of node i */
  Set_Input_i(j,input_present,m_p,n_input,input,ergm,ls); /* Set rows of input corresponding to changing indicators of node j */
  Set_Parameter(ergm->d,ergm->structural,ergm->theta,theta); 
  generating = 0;
  Full_Conditional_Indicator(heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite,nterms,funnames,
                            sonames,MHproposaltype,MHproposalpackage,samplesize,burnin,interval, 
                            newnetworkheads,newnetworktails,verbose,attribs,maxout,maxin,minout,
                            minin,condAllDegExact,attriblength,maxedges,mheads,mtails,mdnedges,
                            number,n_input,input,ergm->d,theta,generating,p); /* Approximate full conditional of category indicator of node i */
  l = 0; /* 0 correspond to present indicators */
  for (k = 0; k < ls->number; k++)
    {
    l = l + 1;
    p_i[k] = p[l]; /* Get unnormalized full conditional of node h */
    }
  Gibbs_Indicator_i(h,p_i,ls); /* Sample full conditional of indicator of h */
  for (k = 0; k < ls->number; k++)
    {
    l = l + 1;
    p_i[k] = p[l]; /* Get unnormalized full conditional of node i */
    }
  Gibbs_Indicator_i(i,p_i,ls); /* Sample full conditional of indicator of i */
  for (k = 0; k < ls->number; k++)
    {
    l = l + 1;
    p_i[k] = p[l]; /* Get unnormalized full conditional of node j */
    }
  Gibbs_Indicator_i(j,p_i,ls); /* Sample full conditional of indicator of j */
}

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
                       double *input_present, double *theta)
/*
input: two nodes, latent structure
output: category indicators
*/
{
  int generating, k, k_i, k_j, l, n_input, number, present, present_i, present_j, **set;
  double **input, *p;
  n_input = Number_Input(ergm->terms,input_present); /* Number of input parameters */
  number = ls->number * ls->number; /* Number of possible values of category indicators */
  input = DD(number,n_input);
  p = D(number);
  set = II(number,2); 
  present_i = ls->indicator[i];
  present_j = ls->indicator[j];
  k = -1;
  for (k_i = 0; k_i < ls->number; k_i++) /* Set input for each possible category indicator of node i */
    {
    for (k_j = 0; k_j < ls->number; k_j++) /* Set input for each possible category indicator of node j */
      {
      k = k + 1;
      if ((k_i == present_i) && (k_j == present_j)) present = k;
      set[k][0] = k_i;
      set[k][1] = k_j;
      ls->indicator[i] = k_i; /* Set indicators */
      ls->indicator[j] = k_j;
      Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input_present given latent structure  */
      for (l = 0; l < n_input; l++) 
        {
        input[k][l] = input_present[l]; /* Store input_present in row k of input */
        }
      }
    }
  ls->indicator[i] = present_i; /* Reset indicators */
  ls->indicator[j] = present_j;
  Set_Parameter(ergm->d,ergm->structural,ergm->theta,theta); 
  generating = present;
  Full_Conditional_Indicator(heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite,nterms,funnames,
                            sonames,MHproposaltype,MHproposalpackage,samplesize,burnin,interval, 
                            newnetworkheads,newnetworktails,verbose,attribs,maxout,maxin,minout,
                            minin,condAllDegExact,attriblength,maxedges,mheads,mtails,mdnedges,
                            number,n_input,input,ergm->d,theta,generating,p); /* Approximate full conditional of category indicator of node i */
  k = Sample_Discrete(p); /* Sample full conditional of category indicators */
  k_i = set[k][0]; /* New category indicators */
  k_j = set[k][1];
  ls->indicator[i] = k_i; /* Store new category indicators */
  ls->indicator[j] = k_j;
  ls->size[present_i] = ls->size[present_i] - 1; /* Decrement sizes */
  ls->size[present_j] = ls->size[present_j] - 1;
  ls->size[k_i] = ls->size[k_i] + 1; /* Increment sizes */
  ls->size[k_j] = ls->size[k_j] + 1; 
}

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
                       double *input_present, double *theta)
/*
input: three nodes, latent structure
output: category indicators
*/
{
  int generating, k, k_h, k_i, k_j, l, n_input, number, present, present_h, present_i, present_j, **set;
  double **input, *p;
  n_input = Number_Input(ergm->terms,input_present); /* Number of input parameters */
  number = ls->number * ls->number * ls->number; /* Number of possible values of category indicators */
  input = DD(number,n_input);
  p = D(number);
  set = II(number,3);
  present_h = ls->indicator[h];
  present_i = ls->indicator[i];
  present_j = ls->indicator[j];
  k = -1;
  for (k_h = 0; k_h < ls->number; k_h++) /* Set input for each possible category indicator of node h */
    {
    for (k_i = 0; k_i < ls->number; k_i++) /* Set input for each possible category indicator of node i */
      {
      for (k_j = 0; k_j < ls->number; k_j++) /* Set input for each possible category indicator of node j */
        {
        k = k + 1;
        if ((present_h == k_h) && (present_i == k_i) && (present_h == k_j)) present = k; 
        set[k][0] = k_h;
        set[k][1] = k_i;
        set[k][2] = k_j;
        ls->indicator[h] = k_h;
        ls->indicator[i] = k_i;
        ls->indicator[j] = k_j;
        Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input_present given latent structure */
        for (l = 0; l < n_input; l++) 
          {
          input[k][l] = input_present[l]; /* Store input_present in row k of input */
          }
        }
      }
    }
  ls->indicator[h] = present_h; /* Reset indicators */
  ls->indicator[i] = present_i;
  ls->indicator[j] = present_j;
  Set_Parameter(ergm->d,ergm->structural,ergm->theta,theta); 
  generating = present;
  Full_Conditional_Indicator(heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite,nterms,funnames,
                            sonames,MHproposaltype,MHproposalpackage,samplesize,burnin,interval, 
                            newnetworkheads,newnetworktails,verbose,attribs,maxout,maxin,minout,
                            minin,condAllDegExact,attriblength,maxedges,mheads,mtails,mdnedges,
                            number,n_input,input,ergm->d,theta,generating,p); /* Approximate full conditional of category indicator of node i */
  k = Sample_Discrete(p); /* Sample full conditional of category indicators */
  k_h = set[k][0]; /* New category indicators */  
  k_i = set[k][1]; 
  k_j = set[k][2];
  ls->indicator[h] = k_h; /* Store new category indicators */
  ls->indicator[i] = k_i; 
  ls->indicator[j] = k_j;
  ls->size[present_h] = ls->size[present_h] - 1; /* Decrement sizes */
  ls->size[present_i] = ls->size[present_i] - 1;
  ls->size[present_j] = ls->size[present_j] - 1;
  ls->size[k_h] = ls->size[k_h] + 1; /* Increment sizes */
  ls->size[k_i] = ls->size[k_i] + 1; 
  ls->size[k_j] = ls->size[k_j] + 1; 
}

void P_Independence(int *number_terms, int *number_parameters, double *input, double *theta,  int *n, int *flag, int *bipartite, char **funnames, char **sonames, double *p)
/*
input: number of terms; number of parameters;  input vector; parameter vector; number of nodes; other variables
output: probabilities of edges between nodes i and j on log scale, computed under the assumption of conditional dyad-independence given latent structure,
and ordered in accordance with i < j
*/
{
  int one = 1;
  int h, i, j, *number_edges, *heads, *tails;
  double log_odds, *statistic;
  number_edges = &one;
  statistic = D(*number_parameters);
  /* 
  Note 1: if i < j, edge (i, j) should be stored as (i, j)
  Note 2: i, j should be integers between 1 and n
  */  
  h = -1;
  for (i = 1; i < *n + 1; i++)
    {
    for (j = i + 1; j < *n + 1; j++)
      {
      h = h + 1;
      heads = &i; 
      tails = &j;
      log_odds = Minus_Energy(*number_parameters,input,theta,heads,tails,number_edges,n,flag,bipartite,number_terms,funnames,sonames,statistic); /* Compute log-odds of probability of edge statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
      p[h] = -ln(1.0 + e(-log_odds));
      }
    }
}

double Partition_Function_Independence(latentstructure *ls, ergmstructure *ergm, double *input, double *theta, 
                                         int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function on log scale, computed under the assumption of conditional dyad-independence given latent structure
*/
{
  int i, j, edges, *n_edges, *heads, *tails;
  double a, b, *statistic;
  edges = 1;
  n_edges = &edges;
  /*
  Rprintf("\n\nNumber of edges: %i",*n_edges);
  */
  heads = I(*n_edges);
  tails = I(*n_edges);
  statistic = D(ergm->d);
  a = 0; /* Log partition function */
  for (i = 1; i < (ls->n + 1); i++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
    {
    heads[0] = i; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
    for (j = (i + 1); j < (ls->n + 1); j++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
      {
      tails[0] = j;
      b = Minus_Energy(ergm->d,input,theta,heads,tails,n_edges,n,dflag,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
      a = a + ln(1.0 + e(b));
      }
    }
  return a;
}

double Partition_Function_Independence_Node(int node, latentstructure *ls, ergmstructure *ergm, double *input, double *theta, 
                                         int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: partition function of node i on log scale, computed under the assumption of conditional dyad-independence given latent structure
*/
{
  int i, j, edges, *n_edges, *heads, *tails;
  double a, b, *statistic;
  edges = 1;
  n_edges = &edges;
  /*
  Rprintf("\n\nNumber of edges: %i",*n_edges);
  */
  heads = I(*n_edges);
  tails = I(*n_edges);
  statistic = D(ergm->d);
  a = 0; /* Log partition function */
  i = node + 1; /* The passed argument node of PMF_Independence_Node is in 0..n-1, whereas Partition_Function_Independence_Node assumes that it is in 1..n */
  tails[0] = i; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
  for (j = 1; j < i; j++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
    {
    heads[0] = j;
    b = Minus_Energy(ergm->d,input,theta,heads,tails,n_edges,n,dflag,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
    a = a + ln(1.0 + e(b));
    }
  heads[0] = i; /* If i < j, edge (i, j) should be stored as (i, j) rather than (j, i): see MCMC.c */
  for (j = (i + 1); j < (ls->n + 1); j++) /* Important note: the C/C++ source files of the ergm package label nodes by integers 1..n rather than 0..n-1 */
    {
    tails[0] = j;
    b = Minus_Energy(ergm->d,input,theta,heads,tails,n_edges,n,dflag,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
    a = a + ln(1.0 + e(b));
    }
  return a;
}

double PMF_Independence(latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input, double *theta, 
                        int *n_edges, int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: probability mass on log scale, computed under the assumption of dyad-dependence
*/
{
  double a, log_p, *statistic, u;
  statistic = D(ergm->d);
  u = Minus_Energy(ergm->d,input,theta,heads,tails,n_edges,n,dflag,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> on log scale */
  /*
  Rprintf("\nPMF_Independence: minus potential energy function = %f",- u);
  */
  a = Partition_Function_Independence(ls,ergm,input,theta,n,dflag,bipartite,nterms,funnames,sonames); /* Log partition function */
  /*
  Rprintf("\nPMF_Independence: log partition function = %f",a);
  */
  log_p = u - a; /* Probability mass */
  return log_p;
}

double PMF_Independence_Node(int i, latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input, double *theta, 
                        int *n_edges, int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: input
output: mass of node i on log scale, computed under the assumption of dyad-dependence
*/
{
  double a, log_p, *statistic, u;
  statistic = D(ergm->d);
  u = Minus_Energy(ergm->d,input,theta,heads,tails,n_edges,n,dflag,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> on log scale */
  /*
  Rprintf("\nPMF_Independence: potential energy function = %f",- u);
  */
  a = Partition_Function_Independence_Node(i,ls,ergm,input,theta,n,dflag,bipartite,nterms,funnames,sonames); /* Log partition function */
  /*
  Rprintf("\nPMF_Independence: log partition function = %f",a);
  */
  log_p = u - a; /* Mass of node i */
  return log_p;
}

double PMF_i_k(int i, int l, latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input_proposal, 
               int *n_edges, int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: node i, catogory l, latent structure, ergm structure
output: conditional PMF of graph given latent structure 
*/
{
  int k, *proposal;
  double log_p_i_k, *theta;
  theta = D(ergm->d);
  k = ls->indicator[i]; /* Store indicator */
  ls->indicator[i] = l; /* Set indicator */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_proposal); /* Set input given indicator  */
  Set_Parameter(ergm->d,ergm->structural,ergm->theta,theta); /* Set parameter */
  log_p_i_k = PMF_Independence(ls,ergm,heads,tails,input_proposal,theta,n_edges,n,dflag,bipartite,nterms,funnames,sonames); /* Probability mass under given indicator */
  ls->indicator[i] = k; /* Reset indicator */
  /*
  Rprintf("\nPMF_i_k: %f",log_p_i_k);
  */
  return log_p_i_k;
}

double PMF_i_k_Node(int i, int l, latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input_proposal, 
                    int *n_edges, int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: node i, catogory l, latent structure, ergm structure
output: conditional PMF of graph given latent structure 
*/
{
  int k, *proposal;
  double log_p_i_k, *theta;
  theta = D(ergm->d);
  k = ls->indicator[i]; /* Store indicator */
  ls->indicator[i] = l; /* Set indicator */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_proposal); /* Set input given indicator  */
  Set_Parameter(ergm->d,ergm->structural,ergm->theta,theta); /* Set parameter */
  log_p_i_k = PMF_Independence_Node(i,ls,ergm,heads,tails,input_proposal,theta,n_edges,n,dflag,bipartite,nterms,funnames,sonames); /* Probability mass under given indicator */
  ls->indicator[i] = k; /* Reset indicator */
  /*
  Rprintf("\nPMF_i_k: %f",log_p_i_k);
  */
  return log_p_i_k;
}

void Gibbs_Indicators_Independence(latentstructure *ls, ergmstructure *ergm, int *heads, int *tails, double *input_proposal, 
                       int *n_edges, int *n, int *dflag, int *bipartite, int *nterms, char **funnames, char **sonames)
/*
input: latent structure, ergm structure
output: indicators
note: function more efficient than sister function Gibbs_Indicators_Independence
*/
{
  int i, k;
  double center, log_p_i_k, p_i_k, *p_i, sum, *selected;
  selected = D(ls->n);
  /* GetRNGstate(); */
  for (i = 0; i < ls->n; i++)
    {
    selected[i] = unif_rand();
    }
  /* PutRNGstate(); */
  p_i = D(ls->number);
  for (k = 0; k < ls->number; k++) /* Reset size */
    {
    ls->size[k] = 0;
    }
  for (i = 0; i < ls->n; i++) /* Node i */
    {
    if (selected[i] < 0.1) /* Indicator of node i updated: y/n */ /* 333 */
      {
      /*
      Rprintf("\nnode %-3i",i);
      */
      sum = 0;
      for (k = 0; k < ls->number; k++) /* Category k */
        {
        log_p_i_k = PMF_i_k_Node(i,k,ls,ergm,heads,tails,input_proposal,n_edges,n,dflag,bipartite,nterms,funnames,sonames);
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
        Rprintf(" %f",p_i[k]);
        */
        }
      k = Sample_Discrete(p_i); /* Sample full conditional of category indicators */
      ls->indicator[i] = k; /* Update indicator */ 
      }
    else k = ls->indicator[i];
    ls->size[k] = ls->size[k] + 1; /* Update size */
    }   
}

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
                       double *input_proposal, double *input_present, double *theta, int print)
/*
input: latent structure
output: category indicators
*/
{
  int accept, i, j, k, l, present_i, proposal_i, *proposal, *size;
  double log_ergm_ratio, log_ratio, *q_i, *statistic, sum, t, u;
  /* Proposal: */
  q_i = D(ls->number);
  t = 1.25; 
  sum = 0; 
  for (k = 0; k < ls->number; k++)
    {
    q_i[k] = e(t * ln(ls->p[k])); /* Melt multinomial pmf with probability vector ls->p */
    sum = sum + q_i[k];
    }
  for (k = 0; k < ls->number; k++) 
    {
    q_i[k] = q_i[k] / sum; /* Normalize */
    }
  /* Proposal: */ 
  proposal = I(ls->n);
  size = I(ls->number);
  log_ratio = 0;
  for (i = 0; i < ls->n; i++)
    {
    k = ls->indicator[i];
    /* GetRNGstate(); */
    u = unif_rand();
    /* PutRNGstate(); */
    if (u < 0.125) l = Sample_Discrete(q_i);
    else l = k;
    proposal[i] = l;
    size[l] = size[l] + 1;
    /* Add log-ratio of ls_theta pmfs: */
    log_ratio = log_ratio + (ln(q_i[k]) - ln(q_i[l]));
    }
  /* Add log-ratio of multinomial pmfs: */
  for (i = 0; i < ls->n; i++) 
    {
    k = ls->indicator[i]; /* Present category */
    l = proposal[i]; /* Proposed category */
    log_ratio = log_ratio + (ln(ls->p[l]) - ln(ls->p[k]));
    }
  /* Add log-ratio of ergm pmfs: */
  theta = D(ergm->d);
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,proposal,ls->theta,input_proposal); /* Set input given indicator  */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input given ls->indicator */
  Set_Parameter(ergm->d,ergm->structural,ergm->theta,theta); 
  log_ergm_ratio = Ratio_Ergm_Pmfs(heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite,
                   nterms,funnames,sonames,MHproposaltype,MHproposalpackage,samplesize,
                   burnin,interval,newnetworkheads,newnetworktails,verbose,attribs,maxout,maxin,minout,
                   minin,condAllDegExact,attriblength,maxedges,mheads,mtails,mdnedges,
                   input_proposal,input_present,ergm->d,theta,theta); /* Log ratio of ergm pmfs */
  log_ratio = log_ratio + log_ergm_ratio;
  accept = MH_Decision(log_ratio);
  if (accept == 1) /* Proposal accepted */
    {
    Set_I_I(ls->n,ls->indicator,proposal);
    Set_I_I(ls->number,ls->size,size);
    }
  if (print == 1)
    {
    Rprintf("\nSample_Indicators_2");
    Rprintf("\n- log_ratio: %8.4f",log_ratio);
    Rprintf("\n- M-H decision (indicators) = %i",accept);
    }
  return accept;
}

int Sample_Parameters_1(ergmstructure *ergm, latentstructure *ls, priorstructure *prior,
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
                        int *newnetworkheads, int *newnetworktails, int n_between, double scale_factor)
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/
{
  int accept, candidate, h, i, j, k, l, max_edges, n_input, proposal_n_edges, present_i, proposal_i, *proposal_heads, *proposal_tails, *ls_indicator, *ls_size, s;
  double **cf, *present, *ergm_theta, *input_proposal, log_denominator, log_numerator, log_present, log_proposal, log_ratio, **ls_theta, *mean, **precision2, *proposal, *theta_present, *theta_proposal, *q_i, *q_k, *prior_mean, **prior_cf, **prior_precision, *statistic, sum, t, u;
  max_edges = *maxpossibleedges; /* Number of possible edges */
  log_ratio = 0;
  n_input = Number_Input(ergm->terms,input_present);
  input_proposal = D(n_input);
  for (i = 0; i < n_input; i++) input_proposal[i] = input_present[i];
  /* Propose indicators: */
  q_i = D(ls->n);
  q_k = D(ls->number);
  ls_indicator = I(ls->n);
  ls_size = I(ls->number);
  for (i = 0; i < ls->n; i++) /* Discrete uniform */
    {
    q_i[i] = 1.0 / ls->n;
    }
  /* GetRNGstate(); */
  u = unif_rand();
  /* PutRNGstate(); */
  if (u < 0.5) candidate = Sample_Discrete(q_i); /* Sample node */
  else candidate = -1; /* No node sampled */
  for (i = 0; i < ls->number; i++) /* Discrete uniform */
    {
    q_k[i] = 1.0 / ls->number;
    }
  for (i = 0; i < ls->n; i++) /* Given node, sample category */
    {
    if (i == candidate) k = Sample_Discrete(q_k); /* Given node, sample category */ 
    else k = ls->indicator[i];
    ls_indicator[i] = k;
    ls_size[k] = ls_size[k] + 1;
    }
  for (i = 0; i < ls->n; i++) /* Ratio of prior pmfs: ratio of multinomial pmfs */
    {
    k = ls->indicator[i]; /* Present category */
    l = ls_indicator[i]; /* Proposed category */
    log_ratio = log_ratio + (ln(ls->p[l]) - ln(ls->p[k]));
    }
  /*
  Rprintf("\n- log_ratio (block indicators) = %8.4f",log_ratio);  
  */
  /* Propose ergm->theta: */
  if (ergm->d1 > 0) /* Michael: must be changed */
    {
    ergm_theta = D(ergm->d1);
    cf = DD(ergm->d1,ergm->d1);
    Scale(ergm->d1,ergm->d1,prior->cf1,scale_factor,cf); /* Rescale Cholesky factor of Gaussian prior */        
    Sample_MVN(ergm->d1,ergm->theta,cf,ergm_theta); /* Random walk Metropolis-Hastings */
    log_proposal = MVN_PDF(ergm->d1,ergm_theta,prior->mean1,prior->precision1); /* Prior pdf: proposal */
    log_present = MVN_PDF(ergm->d1,ergm->theta,prior->mean1,prior->precision1); /* Prior pdf: present */
    log_ratio = log_ratio + (log_proposal - log_present);
    /*
    Rprintf("\n- log_ratio (parameters) = %8.4f",log_ratio);  
    */
    }
  /* Propose ls->theta: */
  ls_theta = DD(ls->d,ls->number+1);
  proposal = D(ls->d);
  present = D(ls->d);  
  cf = DD(ls->d,ls->d);
  Scale(ls->d,ls->d,prior->cf2,scale_factor,cf); /* Rescale Cholesky factor of Gaussian prior */ 
  for (i = 0; i < ls->number; i++) 
    {
    Get_Column(ls->d,present,ls->theta,i); /* Set mean to ls->theta[][i] */
    if (ls->size[i] < ls->threshold) Set_Column(ls->d,ls_theta,i,present); /* Set proposal = present */ 
    else 
      {
      /* Generate candidate: */
      Sample_MVN(ls->d,present,cf,proposal); /* Random walk Metropolis-Hastings algorithm */
      Set_Column(ls->d,ls_theta,i,proposal); /* Set ls_theta[][i] to proposal */
      /* Add ratio of prior pdf: */
      log_proposal = MVN_PDF(ls->d,proposal,prior->mean2,prior->precision2); /* Prior pdf of proposal */
      log_present = MVN_PDF(ls->d,present,prior->mean2,prior->precision2); /* Prior pdf of present */
      log_ratio = log_ratio + (log_proposal - log_present);
      }
    /*
    Rprintf("\n- log_ratio (block parameters) = %8.4f",log_ratio);  
    */
    }
  /* Assumption:
  - n_between = number of hierarchical ergm terms with between-category parameters
  - n_between hierarchical ergm terms with between-category parameters come first, and all arrays are ordered accordingly
  - note: must not computed under parametric prior, because there is no betweeness parameter and thus adding prior pdf of betweeness parameter would bias conclusions */  
  if (n_between > 0) /* Number of between-category parameters */ 
    {
    prior_mean = D(n_between);
    prior_cf = DD(n_between,n_between);
    prior_precision = DD(n_between,n_between);
    proposal = D(n_between);
    present = D(n_between);  
    for (i = 0; i < n_between; i++)
      {
      prior_mean[i] = prior->mean2[i];
      present[i] = ls->theta[i][ls->number];
      for (j = 0; j < n_between; j++) 
        {
        cf[i][j] = prior->cf2[i][j] * scale_factor;
        if (i == j) prior_precision[i][j] = 1.0 / (prior->cf2[i][j] * prior->cf2[i][j]); /* prior->cf2 is the Cholesky factor: its square is the variance matrix, and the inverse of its square is the precision */
        else prior_precision[i][j] = 0;
        }
      }
    for (i = n_between; i < ls->d; i++) /* To be on the safe side, set unadmissible between-category parameters to 0 */
      {
      ls->theta[i][ls->number] = 0;
      }
    Sample_MVN(n_between,present,cf,proposal); /* Random walk Metropolis-Hastings algorithm */
    for (i = 0; i < n_between; i++) 
      {
      ls_theta[i][ls->number] = proposal[i];
      }
    for (i = n_between; i < ls->d; i++)
      {
      ls_theta[i][ls->number] = 0; 
      }
    }
  log_proposal = MVN_PDF(n_between,proposal,prior_mean,prior_precision); /* Proposal */
  log_present = MVN_PDF(n_between,present,prior_mean,prior_precision); /* Present */
  log_ratio = log_ratio + (log_proposal - log_present);
  /*
  Rprintf("\n- log_ratio = %8.4f",log_ratio);  
  */
  /* Propose auxiliary variable: */
  theta_proposal = D(ergm->d);
  theta_present = D(ergm->d);
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls_indicator,ls_theta,input_proposal); /* Set input given ls_theta */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input given ls->theta */
  Set_Parameter(ergm->d,ergm->structural,ergm_theta,theta_proposal); /* Set parameter; note: if ergm_d1 == 0, ergm_theta is not used */
  Set_Parameter(ergm->d,ergm->structural,ergm->theta,theta_present); /* Set parameter; note: if ergm->d1 == 0, ergm_theta is not used */
  s = 1; /* One sample point is all that is required */
  MCMC_wrapper(heads,tails,dnedges,  /* Sample one graph from posterior predictive distribution given input and theta */
                  maxpossibleedges,
                  dn,dflag,bipartite, 
                  nterms,funnames,
                  sonames, 
                  MHproposaltype,MHproposalpackage,
                  input_proposal,theta_proposal,&s,
                  sample,burnin,interval,  
                  newnetworkheads, 
                  newnetworktails, 
                  verbose, 
                  attribs,maxout,maxin,minout,
                  minin,condAllDegExact,attriblength, 
                  maxedges,
                  mheads,mtails,mdnedges);
  proposal_n_edges = newnetworkheads[0]; /* Number of simulated edges */
  proposal_heads = I(proposal_n_edges); /* Proposed heads for auxiliary variable */
  proposal_tails = I(proposal_n_edges); /* Proposed tails for auxiliary variable */
  for (i = 0; i < proposal_n_edges; i++)  
    {
    proposal_heads[i] = newnetworkheads[i+1]; /* Note: while heads corresponds to the list of observed heads, newnetworkheads contains the number of simulated edges as well as the list of simulated heads: to use auxiliary->heads here, one must not store the number of simulated edges */
    proposal_tails[i] = newnetworktails[i+1]; /* Note: while tails corresponds to the list of observed tails, newnetworktails contains the number of simulated edges as well as the list of simulated tails: to use auxiliary->tails here, one must not store the number of simulated edges */
    }
  /* Ratio of proposal pmfs of auxiliary graph under proposal / present */
  statistic = D(ergm->d);
  log_numerator = Minus_Energy(ergm->d,input_present,theta_present,
  proposal_heads,proposal_tails,&proposal_n_edges,dn,dflag,bipartite,nterms,funnames,sonames,statistic); 
  /*
  Rprintf("\n- function 1: log numerator = %8.4f",log_numerator);  
  */
  log_denominator = Minus_Energy(ergm->d,input_proposal,theta_proposal,
  proposal_heads,proposal_tails,&proposal_n_edges,dn,dflag,bipartite,nterms,funnames,sonames,statistic); 
  /*
  Rprintf("\n- function 2: log denominator = %8.4f",log_denominator);  
  */
  log_ratio = log_ratio + (log_numerator - log_denominator);
  /*
  Rprintf("\n- log_ratio (auxiliary graph) = %8.4f",log_ratio);  
  */
  /* Ratio of mass of observed graph under proposal / present */
  log_present = Minus_Energy(ergm->d,input_present,theta_present,heads,tails,dnedges,dn,dflag,bipartite,nterms,funnames,sonames,statistic); 
  /*
  Rprintf("\n- function 3: log present = %8.4f",log_present);  
  */
  log_proposal = Minus_Energy(ergm->d,input_proposal,theta_proposal,heads,tails,dnedges,dn,dflag,bipartite,nterms,funnames,sonames,statistic); 
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
    if (ergm->d1 > 0) Set_D_D(ergm->d1,ergm->theta,ergm_theta); /* Michael: must be changed */
    Set_DD_DD(ls->d,ls->number+1,ls->theta,ls_theta);
    }
  if (print == 1)
    {
    Rprintf("\nSample parameters and indicators:");
    Rprintf("\n- log_ratio = %8.4f",log_ratio);  
    Rprintf("\n- decision = %i",accept);
    }
  return accept;
}

void Sample_Parameters_2(ergmstructure *ergm, latentstructure *ls, priorstructure *prior)
/*
input: ergm structure, latent structure, prior
output: non-structural parameters not showing up in the ergm pmf
*/
{
  int i;
  double *theta;
  theta = D(ls->d);
  for (i = 0; i < ls->number; i++)
    {
    if (ls->size[i] < ls->threshold) /* Structural parameter not showing up in ergm pmf */
      {
      Sample_MVN(ls->d,prior->mean2,prior->cf2,theta); /* Sample structural parameter from full conditional (conditional Gaussian prior given non-structural parameters) */
      Set_Column(ls->d,ls->theta,i,theta); /* Set ls_theta[][i] to theta */
      }
    } 
}

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
                        double *input_proposal, double *input_present, int print, int n_between, double scale_factor)
/*
input: ergm structure, latent structure, prior
output: structural, non-structural parameters showing up in ergm pmf
*/
{
  int accept, h, i, j;
  double **cf, *present, *ergm_theta, log_present, log_proposal, log_p_present, log_p_proposal, log_denominator, log_numerator, log_ergm_ratio, log_ratio, **ls_theta, *mean, **precision2, *proposal, *theta_present, *theta_proposal;
  /* Note:
  parametric Bayesian prior can be specified by
  - specifying ls->number = 1
  - including edges_ij */
  /* Proposal:
  note 1: all ls->theta such that ls->size >= ls->threshold and all ergm->theta are updated by random walk Metropolis-Hastings algorithm
  note 2: ratio of proposal pdfs cancels under random walk Metropolis-Hastings algorithm */
  log_ratio = 0;
  /* Propose ergm->theta: */
  if (ergm->d1 > 0)
    {
    ergm_theta = D(ergm->d1);
    cf = DD(ergm->d1,ergm->d1);
    Scale(ergm->d1,ergm->d1,prior->cf1,scale_factor,cf); /* Rescale Cholesky factor of Gaussian prior */        
    Sample_MVN(ergm->d1,ergm->theta,cf,ergm_theta); /* Random walk Metropolis-Hastings */
    /* Add ratio of prior pdf: */
    log_proposal = MVN_PDF(ergm->d1,ergm_theta,prior->mean1,prior->precision1); /* Prior pdf: proposal */
    log_present = MVN_PDF(ergm->d1,ergm->theta,prior->mean1,prior->precision1); /* Prior pdf: present */
    log_ratio = log_ratio + (log_proposal - log_present);
    /*
    Rprintf("\n- log_ratio (parameters) = %8.4f",log_ratio);  
    */
    }
  /* Propose ls->theta: */
  ls_theta = DD(ls->d,ls->number+1);
  proposal = D(ls->d);
  present = D(ls->d);  
  cf = DD(ls->d,ls->d);
  Scale(ls->d,ls->d,prior->cf2,scale_factor,cf); /* Rescale Cholesky factor of Gaussian prior */ 
  for (i = 0; i < ls->number; i++) 
    {
    Get_Column(ls->d,present,ls->theta,i); /* Set mean to ls->theta[][i] */
    if (ls->size[i] < ls->threshold) Set_Column(ls->d,ls_theta,i,present); /* Set proposal = present */ 
    else 
      {
      /* Generate candidate: */
      Sample_MVN(ls->d,present,cf,proposal); /* Random walk Metropolis-Hastings algorithm */
      Set_Column(ls->d,ls_theta,i,proposal); /* Set ls_theta[][i] to proposal */
      /* Add ratio of prior pdf: */
      log_proposal = MVN_PDF(ls->d,proposal,prior->mean2,prior->precision2); /* Prior pdf of proposal */
      log_present = MVN_PDF(ls->d,present,prior->mean2,prior->precision2); /* Prior pdf of present */
      log_ratio = log_ratio + (log_proposal - log_present);
      }
    /*
    Rprintf("\n- log_ratio (block parameters) = %8.4f",log_ratio);  
    */
    }
  /* Assumption: there are no between-category parameters */
  proposal = D(ls->d);
  present = D(ls->d);  
  Set_Column(ls->d,ls_theta,ls->number,proposal); /* Set ls_theta[][i] to proposal */
  Set_Column(ls->d,ls->theta,ls->number,present); /* Set ls_theta[][i] to proposal */
  /* Decide: */
  theta_proposal = D(ergm->d);
  theta_present = D(ergm->d);
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls_theta,input_proposal); /* Set input given ls_theta */
  Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,input_present); /* Set input given ls->theta */
  Set_Parameter(ergm->d,ergm->structural,ergm_theta,theta_proposal); /* Set parameter; note: if ergm_d1 == 0, ergm_theta is not used */
  Set_Parameter(ergm->d,ergm->structural,ergm->theta,theta_present); /* Set parameter; note: if ergm->d1 == 0, ergm_theta is not used */  
  /* Compare naive M-H with importance sampler to exact M-H under dyad-independence conditional on latent structure: 
  log_ergm_ratio = Ratio_Ergm_Pmfs(heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite,
                   nterms,funnames,sonames,MHproposaltype,MHproposalpackage,samplesize,
                   burnin,interval,newnetworkheads,newnetworktails,verbose,attribs,maxout,maxin,minout,
                   minin,condAllDegExact,attriblength,maxedges,mheads,mtails,mdnedges,
                   input_proposal,input_present,ergm->d,theta_proposal,theta_present); 
  Rprintf("\nlog_ergm_ratio version A: %f",log_ergm_ratio);
  */
  log_p_proposal = PMF_Independence(ls,ergm,heads,tails,input_proposal,theta_proposal,dnedges,dn,dflag,bipartite,nterms,funnames,sonames); /* Probability mass under proposed parameter */
  log_p_present = PMF_Independence(ls,ergm,heads,tails,input_present,theta_present,dnedges,dn,dflag,bipartite,nterms,funnames,sonames); /* Probability mass under present parameters */
  log_ergm_ratio = log_p_proposal - log_p_present;
  /*
  Rprintf("\nlog_ergm_ratio version B: %f",log_ergm_ratio);
  */
  log_ratio = log_ratio + log_ergm_ratio; 
  accept = MH_Decision(log_ratio);
  if (accept == 1) /* Proposal accepted: set ergm->theta and ls->theta to proposal */
    {
    if (ergm->d1 > 0) Set_D_D(ergm->d1,ergm->theta,ergm_theta); /* Michael: must be changed */
    Set_DD_DD(ls->d,ls->number+1,ls->theta,ls_theta);
    }
  if (print == 1)
    {
    Rprintf("\nSample parameters:");
    Rprintf("\n- log_ratio = %8.4f",log_ratio);  
    Rprintf("\n- decision = %i",accept);
    }
  return accept;
}

void Initial_State(int *parallel, double *alpha, int *indicator, priorstructure_ls *prior_ls, priorstructure *prior, latentstructure *ls, ergmstructure *ergm, double *theta)
/* 
input: clustering parameter, priors, latent structure, ergm structure, user-specified initial value of non-structural parameters
*/
{  
  int i, k;
  double *sample, *shape1, *shape2, sum;
  if (*parallel == 1) ls->alpha = *alpha; /* Clustering parameter */
  else ls->alpha = rgamma(prior_ls->alpha_shape,1.0/prior_ls->alpha_rate); 
  shape1 = D(ls->number-1); /* Components 0..ls->number-2 suffice */
  shape2 = D(ls->number-1); /* Components 0..ls->number-2 suffice */
  for (i = 0; i < (ls->number - 1); i++)
    {
    shape1[i] = 1.0; /* First shape of Beta distribution */
    shape2[i] = ls->alpha; /* Second shape of Beta distribution */
    }
  Stick_Breaking(shape1,shape2,ls); /* Construct category probability vector by stick-breaking */
  for (i = 0; i < ls->n; i++) /* For each node i, sample category k */
    {
    if (*parallel == 1) k = indicator[i]; 
    else k = Sample_Discrete(ls->p);
    ls->indicator[i] = k; 
    ls->size[k] = ls->size[k] + 1; /* ls-size was set to 0 by S_alloc */
    }
  if (*parallel > 1)
    {
    if (ergm->d1 > 0) Sample_MVN(ergm->d1,prior->mean1,prior->cf1,ergm->theta);
    for (i = 0; i < ls->number; i++) 
      {
      sample = D(ls->d);
      Sample_MVN(ls->d,prior->mean2,prior->cf2,sample); /* Random walk Metropolis-Hastings algorithm */
      Set_Column(ls->d,ls->theta,i,sample); /* Set ls_theta[][i] to proposal */
      }
    }
  Get_Parameter(ergm->d,ls->number,ergm->structural,theta,ergm->theta,ls->theta); /* Set ergm->theta and ls->theta parameter values inputted by user */
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
  p = D(ls->n);
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
  return number;
}

int Sample_Graph_Independence(latentstructure *ls, double *ln_p, int *heads, int *tails)
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
  /* GetRNGstate(); */
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
  /* PutRNGstate(); */
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
             int *max_iterations, int *n_between_block_parameters, int *output, double *mcmc, int *sample_heads, int *sample_tails, int *call_RNGstate)
/*
input: R input
output: simulated graph
*/
{
  int null = 0;
  int coordinate, dim, dim1, dim2, edges, element, h, i, *n_edges, *pseudo_indicator, iteration, k, max_iteration, *mdnedges, *mheads, *mtails, n, *newnetworkheads, *newnetworktails, number, print, threshold, terms, *verbose;
  double *ln_p, **parameter, *pp, progress, *shape1, *shape2, sum;	
  char *vmax;
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
  if (print == 1)
    {
    Rprintf("\nMachine precision:");
    Rprintf("\n- epsilon = %e",epsilon);
    Rprintf("\n- maximum = %e",maximum);
    Rprintf("\n- ln(epsilon) = %e",ln(epsilon));
    Rprintf("\n- ln(maximum) = %e",ln(maximum));
    Rprintf("\n- exp(-maximum) = %e",e(-maximum));
    Rprintf("\n- exp(-epsilon)= %e",e(-epsilon));
    Rprintf("\n- exp(+epsilon) = %e",e(+epsilon));
    Rprintf("\n- exp(+maximum) = %e",e(+maximum));
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
  mdnedges = &null;
  ergm = Initialize_Ergm(terms,hierarchical,dim,dim1,dim2,structural); /* Ergm structure and non-structural parameters */
  prior = Initialize_Prior(ergm->d1,ergm->d2,m1,m2,b,cf1,cf2,p1,p2); /* Prior: non-structural, structural parameters */
  ls = Initialize_Latentstructure(number,n,threshold,ergm->d2); /* Latent structure and structural parameters */
  prior_ls = Initialize_Prior_ls(*alpha_shape,*alpha_rate); /* Prior: clustering parameter */
  mheads = NULL;
  mtails = NULL;
  /****************/
  /* Sample graph */
  /****************/
  if (*call_RNGstate == 1) GetRNGstate();
  ls->alpha = *alpha; /* Clustering parameter */
  coordinate = -1;
  element = -1;
  for (iteration = 0; iteration < max_iteration; iteration++)
    {
    /* vmax = vmaxget(); */
    progress = (iteration * 100.0) / max_iteration;
    if (print == 1) Rprintf("\nProgress: %5.2f%%",progress);
    shape1 = D(ls->number-1); /* Components 0..ls->number-2 suffice */
    shape2 = D(ls->number-1); /* Components 0..ls->number-2 suffice */
    for (i = 0; i < (ls->number - 1); i++)
      {
      shape1[i] = 1.0; /* First shape of Beta distribution */
      shape2[i] = ls->alpha; /* Second shape of Beta distribution */
      }
    Stick_Breaking(shape1,shape2,ls); /* Construct category probability vector by stick-breaking */
    Sample_CRP(ls);
    if (ergm->d1 > 0) Sample_MVN(ergm->d1,prior->mean1,prior->cf1,ergm->theta);
    for (i = 0; i < ls->number; i++) 
      {
      sample = D(ls->d);
      Sample_MVN(ls->d,prior->mean2,prior->cf2,sample); /* Random walk Metropolis-Hastings algorithm */
      Set_Column(ls->d,ls->theta,i,sample); /* Set ls_theta[][i] to proposal */
      }
    pp = D(ergm->d);
    Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,inputs);
    Set_Parameter(ergm->d,ergm->structural,ergm->theta,theta); 
    newnetworkheads = I(*maxpossibleedges);
    newnetworktails = I(*maxpossibleedges);
    if (*dyaddependence == 0) 
      { 
      ln_p = D(*maxpossibleedges);
      P_Independence(nterms,d,inputs,theta,dn,dflag,bipartite,funnames,sonames,ln_p);
      Sample_Graph_Independence(ls,ln_p,newnetworkheads,newnetworktails);
      pseudo_indicator = I(ls->n);
      for (i = 0; i < ls->n; i++) /* Identical indicators */
        {
        pseudo_indicator[i] = 1;
        }
      parameter = DD(ls->d,ls->number+1);
      for (i = 0; i < ls->d; i++) /* Identical structural parameters, so that structural function of graph, structural parameters reduces to corresponding non-structural function of graph */
        {
        for (k = 0; k < (ls->number + 1); k++)
          {
          parameter[i][k] = 1.0; 
          }
        }
      Set_Input(terms,hierarchical,number,ls->n,pseudo_indicator,parameter,inputs);
      n_edges = newnetworkheads; /* First element of newnetworkheads = newnetworkheads[0] is number of edges */
      mheads = I(*n_edges);
      mtails = I(*n_edges);
      for (i = 0; i < *n_edges; i++) /* Since first element of newnetworkheads and newnetworktails is number of edges, heads and tails must be extracted */
        {
        mheads[i] = newnetworkheads[i+1];
        mtails[i] = newnetworktails[i+1];
        }
      network_stats_wrapper(mheads,mtails,n_edges,dn,dflag,bipartite,nterms,funnames,sonames,inputs,pp); /* Compute non-structural function of graph */
      }
    else Sample_Graph(ls->number,ls->n,ls->d,ergm->terms,ergm->hierarchical,ergm->d,pp, /* Sample graph */
                         heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite,nterms,funnames,
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
      if (print == 1) Rprintf(" %4i",ls->indicator[i]+1);
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
    /* vmaxset(vmax); */
    }
  if (*call_RNGstate == 1) PutRNGstate();
  /************/
  /* Finalize */
  /************/
  Rprintf("\n");
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
             int *max_iterations, int *n_between_block_parameters, int *output, double *mcmc, double *scalefactor, double *mh_accept, int *call_RNGstate, int *parallel)
/*
input: R input
output: MCMC sample of unknowns from posterior
*/
{
  int null = 0;
  int batch, n_batches, batch_size, coordinate, console, dependence, dim, dim1, dim2, h, i, j, k, *mdnedges, *mheads, *mtails, n_input, iteration, max_iteration, n, n_between, number, print, store, threshold, terms, *verbose;
  double accept, accept_ergm_theta, *accept_ls_theta, accept_indicators, local_mh_accept, *pp, progress, rate, shape, scale_factor, u;	
  char *vmax; 
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
  if (console == 1)
    {
    Rprintf("\nMachine precision:");
    Rprintf("\n- epsilon = %e",epsilon);
    Rprintf("\n- maximum = %e",maximum);
    Rprintf("\n- ln(epsilon) = %e",ln(epsilon));
    Rprintf("\n- ln(maximum) = %e",ln(maximum));
    Rprintf("\n- exp(-maximum) = %e",e(-maximum));
    Rprintf("\n- exp(-epsilon)= %e",e(-epsilon));
    Rprintf("\n- exp(+epsilon) = %e",e(+epsilon));
    Rprintf("\n- exp(+maximum) = %e",e(+maximum));
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
  if (max_iteration <= 6000) n_batches = max_iteration;
  else n_batches = 6000;
  if (max_iteration == n_batches) batch_size = 1; 
  else batch_size = trunc(max_iteration / n_batches) + 1;
  n_between = (int)*n_between_block_parameters; /* Number of between-category parameters */
  mdnedges = &null;
  scale_factor = *scalefactor; /* Metropolis-Hasting algorithm: scale factor */
  ergm = Initialize_Ergm(terms,hierarchical,dim,dim1,dim2,structural); /* Ergm structure and non-structural parameters */
  prior = Initialize_Prior(ergm->d1,ergm->d2,m1,m2,b,cf1,cf2,p1,p2); /* Prior: non-structural, structural parameters */
  ls = Initialize_Latentstructure(number,n,threshold,ergm->d2); /* Latent structure and structural parameters */
  prior_ls = Initialize_Prior_ls(*alpha_shape,*alpha_rate); /* Prior: clustering parameter */
  dependence = *dyaddependence; /* Conditional PMF of graph given latent structure: dyad-dependent or not */
  if (console >= 0)
    {
    if (dependence == 0) Rprintf("\nConditional dyad-independence model.\n");
    else Rprintf("\nConditional dyad-dependence model.\n");
    }
  /*************************/
  /* MCMC sample posterior */
  /*************************/
  if (*call_RNGstate == 1) GetRNGstate();
  if (console >= 0)
    {
    Rprintf("\nNumber of draws from posterior: %i",n_batches * batch_size);
    Rprintf("\nNumber of batches: %i",n_batches);
    Rprintf("\nSize of batches: %i",batch_size);
    Rprintf("\n");
    }
  Initial_State(parallel,alpha,indicator,prior_ls,prior,ls,ergm,theta);
  accept_ergm_theta = 0;
  accept_ls_theta = D(ls->number+1);
  accept_indicators = 0;
  local_mh_accept = 0;
  coordinate = -1;
  for (batch = 0; batch < n_batches; batch++) /* Batch */
    {
    progress = (batch * 100.0) / n_batches;
    if (console == 1) Rprintf("\nProgress: %5.2f%%",progress);
    else if (console == 0) Rprintf("%4.1f%%",progress);
    for (iteration = 0; iteration < batch_size; iteration++) /* Iteration within batch */
      {
      /* vmax = vmaxget(); */
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
      if (dependence != 0) accept = Sample_Parameters_1(ergm,ls,prior, /* Auxiliary-variable M-H */
                          heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite,nterms,funnames,
                          sonames,MHproposaltype,MHproposalpackage,sample,burnin,interval, 
                          verbose,attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,
                          maxedges,mheads,mtails,mdnedges,inputs,print,newnetworkheads,newnetworktails,n_between,scale_factor);
      else accept = Sample_Parameters_3(ergm,ls,prior, /* M-H exploiting dyad-independence conditional on latent structure */
                          heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite,nterms,funnames,sonames, 
                          MHproposaltype,MHproposalpackage,samplesize,burnin,interval,
                          newnetworkheads,newnetworktails,verbose,attribs,maxout,maxin,minout,minin,condAllDegExact,
                          attriblength,maxedges,mheads,mtails,mdnedges,inputs,inputs_h,print,n_between,scale_factor); 
      accept_ergm_theta = accept_ergm_theta + accept;
      local_mh_accept = local_mh_accept + accept;
      if (ls->number > 1) /* Non-parametric prior, otherwise paramtric prior */
        {
        Sample_Parameters_2(ergm,ls,prior); /* Structural parameters not showing up in ergm pmf */
        for (i = 0; i < ls->number; i++) /* Acceptance rate of Metropolis-Hastings algorithm for updating structural parameters */
          {
          if (ls->size[i] >= ls->threshold) 
            {
            accept_ls_theta[i] = accept_ls_theta[i] + accept;
            }
          else
            {
            accept_ls_theta[i] = accept_ls_theta[i] + 1; /* Sampled from full conditional */
            }
          }
        accept_ls_theta[ls->number] = accept_ls_theta[ls->number] + 1;
        if (dependence == 0) 
          {
          Gibbs_Indicators_Independence(ls,ergm,heads,tails,inputs_h,dnedges,dn,dflag,bipartite,nterms,funnames,sonames); 
          accept_indicators = accept_indicators + 1.0;
          }
        Sample_P(ls); /* Category probability vector */ 
        Sample_Alpha(prior_ls,ls); /* Clustering parameter */
        if (ls->alpha < epsilon) ls->alpha = epsilon;
        else if (ls->alpha > maximum) ls->alpha = maximum;
        }
      if (store == 1) 
        {
        /* Posterior predict when output is to be created */
        pp = D(ergm->d); /* Memory must be allocated here */
        if (*output == 1)
          {
          theta = D(ergm->d);
          Set_Input(ergm->terms,ergm->hierarchical,ls->number,ls->n,ls->indicator,ls->theta,inputs);
          Set_Parameter(ergm->d,ergm->structural,ergm->theta,theta); 
          Sample_Graph(ls->number,ls->n,ls->d,ergm->terms,ergm->hierarchical,ergm->d,pp, /* Posterior prediction of graph */
                           heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite,nterms,funnames,
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
          if (print == 1) Rprintf(" %4i",ls->indicator[i]+1);
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
      /* vmaxset(vmax); */
      }
    }
  if (*call_RNGstate == 1) PutRNGstate();
  /************/
  /* Finalize */
  /************/
  accept_ergm_theta = accept_ergm_theta / (n_batches * batch_size); /* Acceptance rate of M-H algorithm for updating non-structural, structural parameters showing up in ergm pmf */
  for (i = 0; i < (ls->number + 1); i++)
    {
    accept_ls_theta[i] = accept_ls_theta[i] / (n_batches * batch_size); /* Acceptance rate of M-H algorithm for updating non-structural, structural parameters showing up in ergm pmf */
    }
  accept_indicators = accept_indicators / (n_batches * batch_size); /* Acceptance rate of M-H algorithm for updating indicators */
  local_mh_accept = local_mh_accept / max_iteration;
  mh_accept[0] = local_mh_accept;
  if (console >= 0)
    {
    Rprintf("\n");
    Rprintf("\nNumber of draws from posterior: %i",n_batches * batch_size);
    Rprintf("\nThinning: every %i-th draw recorded",batch_size);
    Rprintf("\nAcceptance rate of Metropolis-Hastings algorithm: %6.4f",local_mh_accept);
    }
  /*
  UNPROTECT(1);
  */
}

