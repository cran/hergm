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

#include "h_ergm_utils.h"
#include "MCMC.h"
#include "netstats.h"

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

double* Set_Input_Block(int terms, int *hierarchical, int max_number, int n, int n_block, double *theta_block, double *input);
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

