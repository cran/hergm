/*
 *  File ergm/src/netstats.h
 *  Part of the statnet package, http://statnetproject.org
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) in
 *    http://statnetproject.org/attribution
 *
 *  Copyright 2010 the statnet development team
 */
#ifndef NETSTATS_H
#define NETSTATS_H

#include "edgetree.h"
#include "model.h"

void network_stats_wrapper(int *heads, int *tails, int *dnedges, 
			   int *dn, int *dflag,  int *bipartite,
			   int *nterms, char **funnames,
			   char **sonames, double *inputs,  double *stats);
void SummStats(Edge n_edges, Vertex *heads, Vertex *tails,
	       Network *nwp, Model *m, double *stats);
#endif
