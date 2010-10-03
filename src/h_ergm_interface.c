#include "h_ergm_interface.h"

int Number_Input(int terms, double *input)
/*
input: number of ergm terms, input parameters
output: number of input parameters
*/
{
  int i, k, number;
  k = -1; 
  for (i = 0; i < terms; i++) 
    {                        
    k = k + 3; 
    number = input[k]; /* Element 3: total number of input parameters */
    k = k + number;
    }
  k = k + 1; /* Since input starts at 0, number of elements is not k but k + 1 */
  return k;
}

void Set_Input(int terms, int *hierarchical, int max_number, int n, int *indicator, double **theta, double *input)
/*
input: number of ergm terms, indicator of hierarchical ergm terms, number of nodes, node-bound category indicator, category-bound parameter
output: input parameters
*/
{
  int h, i, j, k, number;
  h = -1; /* Hierarchical ergm term h */
  k = -1; /* Input parameter k */
  for (i = 0; i < terms; i++) /* For given ergm term... */
    {
    if (hierarchical[i] == 0) /* ...if non-hierarchical, go to following ergm term */
      {
      k = k + 3; /* Elements 1, 2, 3 unchanged */
      number = input[k]; /* Element 3: total number of input parameters */
      k = k + number; /* If number > 0, elements 4, ..., 3 + number unchanged */ 
      }
    else /* ...if hierarchical, set input parameters */ 
      {
      h = h + 1; /* Hierarchical ergm term h */
      k = k + 1; 
      input[k] = 0; /* Element 1 */
      k = k + 1;
      input[k] = 1; /* Element 2: one change statistic */
      k = k + 1;
      input[k] = 1.0 + n + (max_number + 1); /* Element 3: total number of input parameters: (maximum) number of categories, n node-bound category indicators, (max_number + 1) category-bound parameters */
      k = k + 1; 
      input[k] = max_number; /* Elements 4: (maximum) number of categories */
      for (j = 0; j < n; j++) /* Elements 4 + n: category indicators */
        {
        k = k + 1;
        input[k] = indicator[j]; 
        }
      for (j = 0; j < max_number; j++) /* Elements 4 + n + max_number: within-category parameters */
        {
        k = k + 1; 
        input[k] = theta[h][j]; 
        }
      k = k + 1;
      input[k] = theta[h][max_number]; /* Element 4 + n + max_number + 1: between-category parameter */       
      }
    }
}

void Set_Parameter(int d, int *structural, double *theta, double *parameter)
/*
input: number of ergm terms, indicator of structural parameters
output: parameter
*/
{
  int i, k;
  k = -1;
  for (i = 0; i < d; i++)
    {
    if (structural[i] == 0) /* Non-structural parameter */
      {
      k = k + 1;
      parameter[i] = theta[k];
      }
    else /* Structural parameter */
      {
      parameter[i] = 1; /* Structural parameters enter ergm pmf through input parameters of "change statistics" */
      }
    }
}

void Get_Parameter(int d, int max_number, int *structural, double *parameter, double *ergm_theta, double **ls_theta)
/*
input: number of ergm terms, number of categories, indicator of structural parameters, initial parameter specified by user
output: non-structural, structural parameters
*/
{
  int i, k, k1, k2;
  k1 = -1;
  k2 = -1;
  for (i = 0; i < d; i++)
    {
    if (structural[i] == 0) /* Non-structural parameter */
      {
      k1 = k1 + 1;
      ergm_theta[k1] = parameter[i];
      }
    else /* Structural parameter */
      {
      k2 = k2 + 1;
      for (k = 0; k < max_number; k++) /* Set structural parameter for all categories */
        {
        ls_theta[k2][k] = parameter[i];
        }
      if (k2 == 0) ls_theta[k2][max_number] = parameter[i] - 1; /* Between-category parameter of hierarchical ergm term edges_ij */
      else ls_theta[k2][max_number] = 0; /* Between-category parameter of all other hierarchical ergm terms */
      }
    }
}

void MCMCSample_h (char *MHproposaltype, char *MHproposalpackage,
  double *theta, double *networkstatistics, 
  long int samplesize, long int burnin, 
  long int interval, int hammingterm, int fVerbose,
  Network *nwp, Model *m, DegreeBound *bd,
  double *input,
  int nmax,
  int *newnetworkheads, int *newnetworktails, int *n_nodes, int *dflag, int *bipartite,
  int *n_terms, char **funnames, char **sonames, double *input_h, double *networkstatistics_h) {
  long int staken, tottaken, ptottaken, originterval;
  int i, j, components, diam;
  MHproposal MH;
  /* Michael */
  int k;
  Edge e;
  Edge *n_edges; 
  double *statistics_h;
  statistics_h = D(m->n_stats);

  originterval = interval;
  components = diam = 0;
  nwp->duration_info.MCMCtimer=0;
  
/*  if (fVerbose) {
    Rprintf("Simulating %d stats on %ld networks using %s",
             m->n_stats, burnin + samplesize*interval, MHproposaltype);
  } */
  MH_init(&MH,
	  MHproposaltype, MHproposalpackage,
	  fVerbose,
	  nwp, bd);

  /*********************
  networkstatistics are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats networkstatistics should 
  all be zero
  *********************/
/*for (j=0; j < m->n_stats; j++) */
/*  networkstatistics[j] = 0.0; */
/* Rprintf("\n"); */
/* for (j=0; j < m->n_stats; j++){ */
/*   Rprintf("j %d %f\n",j,networkstatistics[j]); */
/* } */
/* Rprintf("\n"); */

  /*********************
   Burn in step.  While we're at it, use burnin statistics to 
   prepare covariance matrix for Mahalanobis distance calculations 
   in subsequent calls to M-H
   *********************/
/*  Catch massive number of edges caused by degeneracy */
   if(nwp->nedges > (50000-1000)){burnin=1;}
   MetropolisHastings(&MH, theta, networkstatistics, burnin, &staken,
		      hammingterm, fVerbose, nwp, m, bd);  
   /* Michael */
   e = EdgeTree2EdgeList((Vertex*)newnetworkheads,(Vertex*)newnetworktails,(Network*)nwp,nmax); /* Extract number of edges and heads and tails from sampled network */
   n_edges = &e;
   for (j = 0; j < m->n_stats; j++) 
     {
     statistics_h[j] = 0;
     }
   network_stats_wrapper(newnetworkheads,newnetworktails,n_edges,n_nodes,dflag,bipartite,n_terms,funnames,sonames,input_h,statistics_h); /* Compute statistics of sampled network given input parameters */ 
   for (j = 0; j < m->n_stats; j++)
     {
     networkstatistics_h[j] = statistics_h[j];
     }  
/*   if (fVerbose){ 
       Rprintf(".");
     } */
  
  if (samplesize>1){
    staken = 0;
    tottaken = 0;
    ptottaken = 0;
    
    /* Now sample networks */
    k = m->n_stats; /* Michael */
    for (i=1; i < samplesize; i++){
      /* Set current vector of stats equal to previous vector */
      for (j=0; j<m->n_stats; j++){
        networkstatistics[j+m->n_stats] = networkstatistics[j];
      }
      networkstatistics += m->n_stats;
      /* This then adds the change statistics to these values */
      
      /* Catch massive number of edges caused by degeneracy */
      if(nwp->nedges > (50000-1000)){interval=1;}
      MetropolisHastings (&MH, theta, networkstatistics, interval, &staken,
                           hammingterm, fVerbose, nwp, m, bd);
      /* Michael */
      e = EdgeTree2EdgeList((Vertex*)newnetworkheads,(Vertex*)newnetworktails,(Network*)nwp,nmax); /* Extract number of edges and heads and tails from sampled network */
      n_edges = &e;
      for (j = 0; j < m->n_stats; j++) 
        {
        statistics_h[j] = 0;
        }
      network_stats_wrapper(newnetworkheads,newnetworktails,n_edges,n_nodes,dflag,bipartite,n_terms,funnames,sonames,input_h,statistics_h); /* Compute statistics of sampled network given input parameters */ 
     for (j = 0; j < m->n_stats; j++)
       {
       networkstatistics_h[k] = statistics_h[j];
       k = k + 1;
       }  
 
      tottaken += staken;

#ifdef Win32
      if( ((100*i) % samplesize)==0 && samplesize > 500){
	R_FlushConsole();
    	R_ProcessEvents();
      }
#endif
      /*if (fVerbose){
        if( ((3*i) % samplesize)==0 && samplesize > 500){
        Rprintf("Sampled %d from Metropolis-Hastings\n", i);}
      }
      if( ((3*i) % samplesize)==0 && tottaken == ptottaken){
        ptottaken = tottaken; 
        Rprintf("Warning:  Metropolis-Hastings algorithm has accepted only "
        "%d steps out of a possible %d\n",  ptottaken-tottaken, i); 
      } 
      if (fVerbose && (i % dotinterval)==0) { 
        Rprintf(".");  
      } */
    }
    /*********************
    Below is an extremely crude device for letting the user know
    when the chain doesn't accept many of the proposed steps.
    *********************/
    if (fVerbose){
      Rprintf("%s sampler accepted %6.3f%% of %d proposed steps.\n",
      MHproposaltype, tottaken*100.0/(1.0*originterval*samplesize), originterval*samplesize); 
    }
  }else{
    if (fVerbose){
      Rprintf("%s sampler accepted %6.3f%% of %d proposed steps.\n",
      MHproposaltype, staken*100.0/(1.0*burnin), burnin); 
    }
  }
  MH_free(&MH);
}

void MCMCSample_h_2 (char *MHproposaltype, char *MHproposalpackage,
  double *theta, double *networkstatistics, 
  long int samplesize, long int burnin, 
  long int interval, int hammingterm, int fVerbose,
  Network *nwp, Model *m, DegreeBound *bd,
  double *input,
  int nmax,
  int *newnetworkheads, int *newnetworktails, int *n_nodes, int *dflag, int *bipartite,
  int *n_terms, char **funnames, char **sonames, int n_number, int n_input_h, double **inputs, double **networkstatistics_h) {
  long int staken, tottaken, ptottaken, originterval;
  int i, j, components, diam;
  MHproposal MH;
  /* Michael */
  int h, k, k_0, number;
  Edge e, *n_edges; 
  double *input_h;
  double *statistics_h;
  input_h = D(n_input_h);
  statistics_h = D(m->n_stats);
  k = 0;

  originterval = interval;
  components = diam = 0;
  nwp->duration_info.MCMCtimer=0;
  
/*  if (fVerbose) {
    Rprintf("Simulating %d stats on %ld networks using %s",
             m->n_stats, burnin + samplesize*interval, MHproposaltype);
  } */
  MH_init(&MH,
	  MHproposaltype, MHproposalpackage,
	  fVerbose,
	  nwp, bd);

  /*********************
  networkstatistics are modified in groups of m->n_stats, and they
  reflect the CHANGE in the values of the statistics from the
  original (observed) network.  Thus, when we begin, the initial 
  values of the first group of m->n_stats networkstatistics should 
  all be zero
  *********************/
/*for (j=0; j < m->n_stats; j++) */
/*  networkstatistics[j] = 0.0; */
/* Rprintf("\n"); */
/* for (j=0; j < m->n_stats; j++){ */
/*   Rprintf("j %d %f\n",j,networkstatistics[j]); */
/* } */
/* Rprintf("\n"); */

  /*********************
   Burn in step.  While we're at it, use burnin statistics to 
   prepare covariance matrix for Mahalanobis distance calculations 
   in subsequent calls to M-H
   *********************/
/*  Catch massive number of edges caused by degeneracy */
   if(nwp->nedges > (50000-1000)){burnin=1;}
   MetropolisHastings(&MH, theta, networkstatistics, burnin, &staken,
		      hammingterm, fVerbose, nwp, m, bd);  
   /* Michael */
   e = EdgeTree2EdgeList((Vertex*)newnetworkheads,(Vertex*)newnetworktails,(Network*)nwp,nmax); /* Extract number of edges and heads and tails from sampled network */
   n_edges = &e;
   k_0 = k; /* Store k so that k can be reset */
   for (number = 0; number < n_number; number++)
     {
     k = k_0; /* Reset k */
     input_h = inputs[number]; /* Set pointer input_h to row number of inputs */
     for (h = 0; h < m->n_stats; h++) /* Reset statistic */
       {
       statistics_h[h] = 0;
       }
     network_stats_wrapper(newnetworkheads,newnetworktails,n_edges,n_nodes,dflag,bipartite,n_terms,funnames,sonames,input_h,statistics_h); /* Compute statistics of sampled network given input parameters */ 
     for (h = 0; h < m->n_stats; h++) /* Store statistic */
       {
       networkstatistics_h[number][k] = statistics_h[h];
       k = k + 1;
       }  
     }
   /* At the end of i-loop, k was incremented by m->n_stats */
/*   if (fVerbose){ 
       Rprintf(".");
     } */
  
  if (samplesize>1){
    staken = 0;
    tottaken = 0;
    ptottaken = 0;
    
    /* Now sample networks */
    for (i=1; i < samplesize; i++){
      /* Set current vector of stats equal to previous vector */
      for (j=0; j<m->n_stats; j++){
        networkstatistics[j+m->n_stats] = networkstatistics[j];
      }
      networkstatistics += m->n_stats;
      /* This then adds the change statistics to these values */
      
      /* Catch massive number of edges caused by degeneracy */
      if(nwp->nedges > (50000-1000)){interval=1;}
      MetropolisHastings (&MH, theta, networkstatistics, interval, &staken,
                           hammingterm, fVerbose, nwp, m, bd);
   /* Michael */
   e = EdgeTree2EdgeList((Vertex*)newnetworkheads,(Vertex*)newnetworktails,(Network*)nwp,nmax); /* Extract number of edges and heads and tails from sampled network */
   n_edges = &e;
   k_0 = k; /* Store k so that k can be reset */
   for (number = 0; number < n_number; number++)
     {
     k = k_0; /* Reset k */
     input_h = inputs[number]; /* Set pointer input_h to row number of inputs */
     for (h = 0; h < m->n_stats; h++) /* Reset statistic */
       {
       statistics_h[h] = 0;
       }
     network_stats_wrapper(newnetworkheads,newnetworktails,n_edges,n_nodes,dflag,bipartite,n_terms,funnames,sonames,input_h,statistics_h); /* Compute statistics of sampled network given input parameters */ 
     for (h = 0; h < m->n_stats; h++) /* Store statistic */
       {
       networkstatistics_h[number][k] = statistics_h[h];
       k = k + 1;
       }  
     }
   /* At the end of i-loop, k was incremented by m->n_stats */
 
      tottaken += staken;

#ifdef Win32
      if( ((100*i) % samplesize)==0 && samplesize > 500){
	R_FlushConsole();
    	R_ProcessEvents();
      }
#endif
      /*if (fVerbose){
        if( ((3*i) % samplesize)==0 && samplesize > 500){
        Rprintf("Sampled %d from Metropolis-Hastings\n", i);}
      }
      if( ((3*i) % samplesize)==0 && tottaken == ptottaken){
        ptottaken = tottaken; 
        Rprintf("Warning:  Metropolis-Hastings algorithm has accepted only "
        "%d steps out of a possible %d\n",  ptottaken-tottaken, i); 
      } 
      if (fVerbose && (i % dotinterval)==0) { 
        Rprintf(".");  
      } */
    }
    /*********************
    Below is an extremely crude device for letting the user know
    when the chain doesn't accept many of the proposed steps.
    *********************/
    if (fVerbose){
      Rprintf("%s sampler accepted %6.3f%% of %d proposed steps.\n",
      MHproposaltype, tottaken*100.0/(1.0*originterval*samplesize), originterval*samplesize); 
    }
  }else{
    if (fVerbose){
      Rprintf("%s sampler accepted %6.3f%% of %d proposed steps.\n",
      MHproposaltype, staken*100.0/(1.0*burnin), burnin); 
    }
  }
  MH_free(&MH);
}

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
                   double *inputs_h, double *sample_h) {
  int directed_flag, hammingterm, formationterm;
  Vertex n_nodes, nmax, bip, hhead, htail;
  Edge n_edges, n_medges, nddyads, kedge;
  Network nw[2];
  DegreeBound *bd;
  Model *m;
  ModelTerm *thisterm;

  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_edges = (Edge)*dnedges; /* coerce double *dnedges to type Edge */
  n_medges = (Edge)*mdnedges; /* coerce double *mdnedges to type Edge */
  nmax = (Edge)*maxedges; /* coerce double *maxedges to type Edge */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */
  
  /* GetRNGstate(); */  
  
  directed_flag = *dflag;

  m=ModelInitialize(*funnames, *sonames, inputs, *nterms);

  /* Form the missing network */
  nw[0]=NetworkInitialize(heads, tails, n_edges, 
                          n_nodes, directed_flag, bip, 0);
  if (n_medges>0) {
   nw[1]=NetworkInitialize(mheads, mtails, n_medges,
                           n_nodes, directed_flag, bip, 0);
  }

  hammingterm=ModelTermHamming (*funnames, *nterms);
  if(hammingterm>0){
/*	     Rprintf("start with setup\n"); */
   Network nwhamming;
   thisterm = m->termarray + hammingterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
   nwhamming=NetworkInitializeD(thisterm->inputparams+1, 
				thisterm->inputparams+1+nddyads, nddyads, 
        n_nodes, directed_flag, bip,0);
   nddyads=0;
   nw[1]=NetworkInitializeD(thisterm->inputparams+1, 
			   thisterm->inputparams+1+nddyads, nddyads,
         n_nodes, directed_flag, bip,0);
/*	     Rprintf("made hw[1]\n"); */
   for (kedge=1; kedge <= nwhamming.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwhamming);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
/*	     Rprintf(" in g0 not g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
   for (kedge=1; kedge <= nw[0].nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nw[0]);
     if(EdgetreeSearch(hhead, htail, nwhamming.outedges) == 0){
/*	     Rprintf("not g0  in g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
/*   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwhamming.nedges,nw[0].nedges); */
   NetworkDestroy(&nwhamming);
  }

/* Really this is a formation term */
  formationterm=ModelTermFormation (*funnames, *nterms);
  if(formationterm>0){
   Network nwformation;
   thisterm = m->termarray + formationterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
   nwformation=NetworkInitializeD(thisterm->inputparams+1,
				  thisterm->inputparams+1+nddyads, nddyads,
          n_nodes, directed_flag, bip,0);
   nddyads=0;
   nw[1]=NetworkInitializeD(thisterm->inputparams+1,
			    thisterm->inputparams+1+nddyads, nddyads,
          n_nodes, directed_flag, bip,0);
/*	     Rprintf("made hw[1]\n"); */
   for (kedge=1; kedge <= nwformation.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwformation);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
/*	     Rprintf(" in g0 not g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[0]);
     }
   }
   for (kedge=1; kedge <= nw[0].nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nw[0]);
     if(EdgetreeSearch(hhead, htail, nwformation.outedges) == 0){
/*	     Rprintf("not g0  in g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
/*   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwformation.nedges,nw[0].nedges); */
   hammingterm=1;
   NetworkDestroy(&nwformation);
/*   Rprintf("Initial number (discord) from reference %d Number of original %d\n",nw[1].nedges,nw[0].nedges); */
  }
  
  bd=DegreeBoundInitialize(attribs, maxout, maxin, minout, minin,
			   *condAllDegExact, *attriblength, nw);

  MCMCSample_h (*MHproposaltype, *MHproposalpackage,
	      theta0, sample, (long int)*samplesize,
	      (long int)*burnin, (long int)*interval,
	      hammingterm,
	      (int)*fVerbose, nw, m, bd, 
              inputs,
              nmax, 
              newnetworkheads, newnetworktails, dn, dflag, bipartite,
              nterms, funnames, sonames, inputs_h, sample_h); /* Michael */

/*   int ii;
   double mos=0.0;
   for(ii=0; ii < bd->attrcount; ii++) 
     mos += bd->maxout[ii];
   Rprintf("bd -> attrcount = %d, sum = %f\n", ii, mos); */
        
        
/* Rprintf("Back! %d %d\n",nw[0].nedges, nmax); */

  /* record new generated network to pass back to R */
  newnetworkheads[0]=newnetworktails[0]=EdgeTree2EdgeList(newnetworkheads+1,newnetworktails+1,nw,nmax-1);
  
  ModelDestroy(m);
  DegreeBoundDestroy(bd);
  NetworkDestroy(nw);
  if (n_medges>0 || hammingterm > 0  || formationterm > 0)
    NetworkDestroy(&nw[1]);
  /* PutRNGstate(); */ 
}

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
                   int n_number, int n_inputs_h, double **inputs_h, double **sample_h) {
  int directed_flag, hammingterm, formationterm;
  Vertex n_nodes, nmax, bip, hhead, htail;
  Edge n_edges, n_medges, nddyads, kedge;
  Network nw[2];
  DegreeBound *bd;
  Model *m;
  ModelTerm *thisterm;
  
  n_nodes = (Vertex)*dn; /* coerce double *dn to type Vertex */
  n_edges = (Edge)*dnedges; /* coerce double *dnedges to type Edge */
  n_medges = (Edge)*mdnedges; /* coerce double *mdnedges to type Edge */
  nmax = (Edge)*maxedges; /* coerce double *maxedges to type Edge */
  bip = (Vertex)*bipartite; /* coerce double *bipartite to type Vertex */
  
  /* GetRNGstate(); */  
  
  directed_flag = *dflag;

  m=ModelInitialize(*funnames, *sonames, inputs, *nterms);

  /* Form the missing network */
  nw[0]=NetworkInitialize(heads, tails, n_edges, 
                          n_nodes, directed_flag, bip, 0);
  if (n_medges>0) {
   nw[1]=NetworkInitialize(mheads, mtails, n_medges,
                           n_nodes, directed_flag, bip, 0);
  }

  hammingterm=ModelTermHamming (*funnames, *nterms);
  if(hammingterm>0){
/*	     Rprintf("start with setup\n"); */
   Network nwhamming;
   thisterm = m->termarray + hammingterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
   nwhamming=NetworkInitializeD(thisterm->inputparams+1, 
				thisterm->inputparams+1+nddyads, nddyads, 
        n_nodes, directed_flag, bip,0);
   nddyads=0;
   nw[1]=NetworkInitializeD(thisterm->inputparams+1, 
			   thisterm->inputparams+1+nddyads, nddyads,
         n_nodes, directed_flag, bip,0);
/*	     Rprintf("made hw[1]\n"); */
   for (kedge=1; kedge <= nwhamming.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwhamming);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
/*	     Rprintf(" in g0 not g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
   for (kedge=1; kedge <= nw[0].nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nw[0]);
     if(EdgetreeSearch(hhead, htail, nwhamming.outedges) == 0){
/*	     Rprintf("not g0  in g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
/*   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwhamming.nedges,nw[0].nedges); */
   NetworkDestroy(&nwhamming);
  }

/* Really this is a formation term */
  formationterm=ModelTermFormation (*funnames, *nterms);
  if(formationterm>0){
   Network nwformation;
   thisterm = m->termarray + formationterm - 1;
   nddyads = (Edge)(thisterm->inputparams[0]);
   nwformation=NetworkInitializeD(thisterm->inputparams+1,
				  thisterm->inputparams+1+nddyads, nddyads,
          n_nodes, directed_flag, bip,0);
   nddyads=0;
   nw[1]=NetworkInitializeD(thisterm->inputparams+1,
			    thisterm->inputparams+1+nddyads, nddyads,
          n_nodes, directed_flag, bip,0);
/*	     Rprintf("made hw[1]\n"); */
   for (kedge=1; kedge <= nwformation.nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nwformation);
     if(EdgetreeSearch(hhead, htail, nw[0].outedges) == 0){
/*	     Rprintf(" in g0 not g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[0]);
     }
   }
   for (kedge=1; kedge <= nw[0].nedges; kedge++) {
     FindithEdge(&hhead, &htail, kedge, &nw[0]);
     if(EdgetreeSearch(hhead, htail, nwformation.outedges) == 0){
/*	     Rprintf("not g0  in g hhead %d htail %d\n",hhead, htail); */
       ToggleEdge(hhead, htail, &nw[1]);
     }
   }
/*   Rprintf("Initial number of discordant %d Number of g0 ties %d Number of ties in g %d\n",nw[1].nedges, nwformation.nedges,nw[0].nedges); */
   hammingterm=1;
   NetworkDestroy(&nwformation);
/*   Rprintf("Initial number (discord) from reference %d Number of original %d\n",nw[1].nedges,nw[0].nedges); */
  }
  
  bd=DegreeBoundInitialize(attribs, maxout, maxin, minout, minin,
			   *condAllDegExact, *attriblength, nw);

  MCMCSample_h_2 (*MHproposaltype, *MHproposalpackage,
	      theta0, sample, (long int)*samplesize,
	      (long int)*burnin, (long int)*interval,
	      hammingterm,
	      (int)*fVerbose, nw, m, bd, 
              inputs,
              nmax, 
              newnetworkheads, newnetworktails, dn, dflag, bipartite,
              nterms, funnames, sonames, n_number, n_inputs_h, inputs_h, sample_h); /* Michael */

/*   int ii;
   double mos=0.0;
   for(ii=0; ii < bd->attrcount; ii++) 
     mos += bd->maxout[ii];
   Rprintf("bd -> attrcount = %d, sum = %f\n", ii, mos); */
        
        
/* Rprintf("Back! %d %d\n",nw[0].nedges, nmax); */

  /* record new generated network to pass back to R */
  newnetworkheads[0]=newnetworktails[0]=EdgeTree2EdgeList(newnetworkheads+1,newnetworktails+1,nw,nmax-1);
  
  ModelDestroy(m);
  DegreeBoundDestroy(bd);
  NetworkDestroy(nw);
  if (n_medges>0 || hammingterm > 0  || formationterm > 0)
    NetworkDestroy(&nw[1]);
  /* PutRNGstate(); */ 
}

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
                          int *mheads, int *mtails, int *mdnedges)
/*
input: (maximum) number of categories, number of nodes, number of structural parameters, number of parameters
output: one sample from posterior predictive distribution
*/
{
  int *h, *t, i, *indicator, k, *nedges, s;
  double **parameter;
  char *vmax_k;
  /*
  vmax_k = vmaxget();
  */
  s = 1; /* Sample one graph from posterior predictive distribution */
  for (i = 0; i < ergm_d; i++)
    {
    sample[i] = 0;
    }
  MCMC_wrapper(heads,tails,dnedges, /* Sample one graph from posterior predictive distribution given input and theta */
               maxpossibleedges,
               dn,dflag,bipartite,
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
               mheads,mtails,mdnedges);
  indicator = I(n);
  for (i = 0; i < n; i++) /* Identical indicators */
    {
    indicator[i] = 1;
    }
  parameter = DD(ls_d,number+1);
  for (i = 0; i < ls_d; i++) /* Identical structural parameters, so that structural function of graph, structural parameters reduces to corresponding non-structural function of graph */
    {
    for (k = 0; k < (number + 1); k++)
      {
      parameter[i][k] = 1; 
      }
    }
  Set_Input(terms,hierarchical,number,n,indicator,parameter,input);
  nedges = newnetworkheads; /* First element of newnetworkheads = newnetworkheads[0] is number of edges */
  h = I(*nedges);
  t = I(*nedges);
  for (i = 0; i < *nedges; i++) /* Since first element of newnetworkheads and newnetworktails is number of edges, heads and tails must be extracted */
    {
    h[i] = newnetworkheads[i+1];
    t[i] = newnetworktails[i+1];
    }
  network_stats_wrapper(h,t,nedges,dn,dflag,bipartite,nterms,funnames,sonames,input,statistic); /* Compute non-structural function of graph */
  /*
  vmaxset(vmax_k);
  */
}

double Minus_Energy(int d, double *input, double *parameter, 
                       int *heads, int *tails, int *nedges, 
		       int *n, int *dflag,  int *bipartite,
		       int *nterms, char **funnames,
		       char **sonames,
                       double *statistic)
/*
input: number of parameters, input parameters, parameters
output: statistic, inner product <parameter, statistic>
*/
{
  int i;
  double sum;
  for (i = 0; i < d; i++) /* Statistic must be null */
    {
    statistic[i] = 0;
    }
  network_stats_wrapper(heads,tails,nedges,n,dflag,bipartite,nterms,funnames,sonames,input,statistic); /* Compute statistic given input */
  sum = 0;
  for (i = 0; i < d; i++)
    {
    sum = sum + (parameter[i] * statistic[i]);
    }
  return sum;
} 

double Ratio_Partition_Functions_1(int s, int d, double sum_observed, double *statistic_generating, double *statistic, double *theta_generating, double *theta)
/*
input: sample size, dimension, difference of inner products under alternative, data-generating parameter for observed graph, value of statistic under data-generating, alternative parameter, value of data-generating, alternative parameter
output: ratio of partition functions of ergms under alternative and data-generating parameter on log scale
*/
{ 
  int i, j, k;
  double ratio_n_const, log, log_generating, log_ratio_n_const, log_ratio, sum;
  /* Ratio of partition functions is expectation of exponential functions:
  - estimated by: MCMC sample average
  - computations of exponential functions stabilized by: multiplying and dividing by e(sum_observed)
  */
  ratio_n_const = 0;
  k = 0; 
  for (i = 0; i < s; i++) /* Sample i */
    {
    sum = 0; 
    for (j = 0; j < d; j++) /* Given sample i... */
      {
      log = theta[j] * statistic[k]; /* ...compute inner product; note indices j, k */
      log_generating = theta_generating[j] * statistic_generating[k]; /* ...compute inner product; note indices j, k */
      sum = sum + (log - log_generating); /* ...compute difference in inner products */
      k = k + 1; /* Get following elements of statistic, statistic_generating */
      }
    /*
    Rprintf("\nLog of exponential function (uncentered) = % -f (centered) = % -f",sum,sum-sum_observed);
    */
    ratio_n_const = ratio_n_const + e(sum - sum_observed); /* e(sum) / e(sum_observed) = e(sum - sum_observed) */
    }
  ratio_n_const = ratio_n_const + 1; /* Add observed term: e(sum_observed - sum_observed) = e(0) = 1 */
  ratio_n_const = ratio_n_const / (s + 1); /* Sample size: s plus observed term */
  log_ratio_n_const = ln(ratio_n_const) + sum_observed; /* On log scale: multiply by e(sum_observed), thus adding ln(e(sum_observed)) = sum_observed */
  return log_ratio_n_const;
}

double Ratio_Partition_Functions_2(int s, int d, double sum_observed, double *statistic_generating, double *statistic, double *theta_generating, double *theta)
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
  moment1 = 0;
  moment2 = 0;
  k = 0; 
  for (i = 0; i < s; i++) /* Sample i */
    {
    sum = 0; 
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
  log_ratio_n_const = moment1 + (variance / 2); /* log expectation = log ratio of partition constants */
  return log_ratio_n_const;
}

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
                   double *input_proposal, double *input_present, int d, double *theta_proposal, double *theta_present)
/*
input: two values of input parameters, two values of parameter
output: ratio of ergm pmfs on log scale
*/
{ 
  int i, j, k, s;
  double ratio_n_const, log_present, log_proposal, log_ratio, log_ratio_n_const, numerator, *statistic_present, *statistic_proposal, *statistic, sum_observed;
  char *vmax_k;
  /*
  vmax_k = vmaxget();
  */
  s = *samplesize;  
  statistic = D(d);
  statistic_present = D(s*d);
  statistic_proposal = D(s*d);
  log_present = Minus_Energy(d,input_present,theta_present,heads,tails,dnedges,dn,dflag,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_present and compute exponential function of inner product <theta_present, statistic> */
  log_proposal = Minus_Energy(d,input_proposal,theta_proposal,heads,tails,dnedges,dn,dflag,bipartite,nterms,funnames,sonames,statistic); /* Compute statistic given input_proposal and compute exponential function of inner product <theta_proposal, statistic> */
  numerator = log_proposal - log_present; 
  sum_observed = log_present - log_proposal; 
  Set_D_D(d,statistic_present,statistic); /* Set elements 0..d-1 of statistic_present to statistic under data-generating value of input_proposal; note: vectors have not same dimension */
  Set_D_D(d,statistic_proposal,statistic); /* Set elements 0..d-1 of statistic_proposal to statistic under data-generating value of input_proposal; note: vectors have not same dimension */
  MCMC_wrapper_h(heads,tails,dnedges,maxpossibleedges,
                 dn,dflag,bipartite,nterms,funnames,sonames, 
                 MHproposaltype,MHproposalpackage,input_proposal,theta_proposal,samplesize, 
                 statistic_proposal,burnin,interval,newnetworkheads,newnetworktails, 
                 fVerbose,attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength, 
                 maxedges,mheads,mtails,mdnedges,input_present,statistic_present); /* Generate MCMC sample of networks given input_proposal and theta_proposal */
  log_ratio_n_const = Ratio_Partition_Functions_2(s,d,sum_observed,statistic_proposal,statistic_present,theta_proposal,theta_present);
  log_ratio = numerator + log_ratio_n_const; /* On log scale: ratio of ergm pmfs; note plus sign */
  /*
  vmaxset(vmax_k);
  */
  return log_ratio;
}

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
                        int number, int n_input, double **input, int d, double *theta, int generating, double *p)
/*
input: number of categories, number of input parameters, dimension of parameter, category indicator of given node to be used to generate importance sample 
output: full conditional of category indicator of given node (given node is specified by calling function and input is set accordingly) 
*/
{ 
  int k, l, s;
  double *input_k, *input_generating, log_k, log_generating, log_ratio, log_ratio_n_const, numerator, ratio_n_const, *statistic_0, *statistic_k, *statistic_generating, **statistic, sum, sum_observed;  
  char *vmax_k;
  /*
  vmax_k = vmaxget();
  */
  /* Full conditional of category indicator of given node:
  - node has been selected in Sample_Indicators_1 and input set in accordance
  - P(Z = k|others) = p_k / sum_l p_l
  - p_l can be written p_l = e_l / c_l, where e_l: exponential function, c_l: partition function
  - c_l = c_0 E_0, where E_0: expectation wrt p_0 and p_0 = e_0 / c_0: ERGM which generates importance sample
  - dividing numerator p_k and denominator sum_l p_l by p_0 = e_0 / c_0 gives p_l / p_0 = (e_l / e_0) (c_0 / c_l) = (e_l / e_0) / E_l
  Approximate full conditional:
  - by generating importance sample under p_0
  - by approximating (e_l / e_0) / E_l 
  - note: p_l must be normalized */
  input_generating = input[generating]; /* Set pointer input_generating to row generating of input */
  s = *samplesize; 
  statistic_0 = D(d);
  statistic_generating = D(s*d);
  statistic = DD(number,s*d);
  log_generating = Minus_Energy(d,input_generating,theta,heads,tails,dnedges,dn,dflag,bipartite,nterms,funnames,sonames,statistic_0); /* Compute statistic given input_generating and compute exponential function of inner product <theta, statistic_0> */
  Set_D_D(d,statistic_generating,statistic_0); /* Set elements 0..d-1 of statistic_generating to statistic under data-generating value of input_generating; note: vectors have not same dimension */
  for (k = 0; k < number; k++)
    {
    for (l = 0; l < d; l++)
      {
      statistic[k][l] = statistic_0[l]; /* Set elements 0..d-1 of statistic to statistic under data-generating value of 
input_generating; note: vectors have not same dimension */
      }
    }
  MCMC_wrapper_h_2(heads,tails,dnedges,maxpossibleedges,dn,dflag,bipartite,
                   nterms,funnames,sonames,MHproposaltype,MHproposalpackage,
                   input_generating,theta,samplesize,statistic_generating,burnin,interval, 
                   newnetworkheads,newnetworktails,fVerbose,attribs,maxout,maxin,minout,
                   minin,condAllDegExact,attriblength,maxedges,mheads,mtails,mdnedges,
                   number,n_input,input,statistic); /* Generate MCMC sample of graph given input_generating and theta */
  sum = 0;
  for (k = 0; k < number; k++) /* Full conditional probability of node i assigned to category k */
    {
    if (k == generating) log_ratio = 0; /* log(p_0 / p_0) = 0 */
    else
      {
      input_k = input[k]; /* Set pointer input_k to row k of input */
      log_k = Minus_Energy(d,input_k,theta,heads,tails,dnedges,dn,dflag,bipartite,nterms,funnames,sonames,statistic_0); /* Compute statistic given input_k and compute exponential function of inner product <theta, statistic_0>; statistic_0 is not used */
      numerator = log_k - log_generating; /* log(e_k / e_0) */
      sum_observed = log_k - log_generating; 
      statistic_k = statistic[k]; /* Set pointer statistic_k to row k of statistic */
      /*
      Rprintf("\nstatistic[0]: % -.4f",statistic_k[0]);
      */
      log_ratio_n_const = Ratio_Partition_Functions_2(s,d,sum_observed,statistic_generating,statistic_k,theta,theta);
      /*
      Rprintf("\nstatistic: % -.4f, % -.4f",statistic_k[0],statistic_k[1]); 
      Rprintf("\nlog numerator: % -.4f, log denominator: % -.4f",numerator,log_ratio_n_const);
      */
      log_ratio = numerator - log_ratio_n_const; /* On log scale: ratio of ergm pmfs; note minus sign */
      }
    p[k] = e(log_ratio); /* Full conditional probability of assigning node i to category k */
    sum = sum + p[k]; 
    /*
    Rprintf("\n sum: % -.4f",sum);
    */
    }
  /*
  Rprintf("\nFull conditional:");
  */
  for (k = 0; k < number; k++) /* Normalize */
    {
    p[k] = p[k] / sum;
    /*
    Rprintf(" % -.4f",p[k]);
    */
    }  
  /*
  vmaxset(vmax_k);
  */
}

