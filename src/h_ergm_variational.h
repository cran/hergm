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

/*
Optimization algorithms based on evaluations of the objective function or its gradient:
#include <math.h>
#include <nlopt.h>
*/

double Expected_Density(int n, double **mu, int directed);
/*
input: number of nodes, expectations of edges under the assumption of independent edges, indicator of directed edges
output: expectation of number of edges under the assumption of independent edges
*/

double Expected_Transitivity(int n, double **mu, int directed);
/*
input: number of nodes, expectations of edges under the assumption of independent edges, indicator of directed edges
output: expected number of transitive triples under the assumption of independent edges
*/

double D_Expected_Transitivity(int n, int i, int j, double **mu, int directed);
/*
input: number nodes, nodes i and j, expectations of edges under the assumption of independent edges, indicator of directed edges
output: expectation of number of two-paths between nodes i and j under the assumption of independent edges
*/

double Expected_Stars(int n, double **mu);
/*
input: number of nodes, expectations of edges under the assumption of independent edges, indicator of directed edges
output: expected number of stars under the assumption of independent edges
*/


double D_Expected_Stars(int n, int i, int j, double **mu);
/*
input: number nodes, nodes i and j, expectations of edges under the assumption of independent edges, indicator of directed edges
output: expectation of number of two-paths between nodes i and j under the assumption of independent edges
*/

double Entropy(int n, double **mu, int directed);
/*
input: number of nodes, expectations of edges under the assumption of independence, indicator of directed edges
output: entropy of distribution under the assumption of independent edges
*/

double Expected_Energy(int n, int model, double *eta, double **mu, int directed);
/*
input: number of nodes, edge and triangle parameters, expectations of edges under the assumption of independent edges, indicator of directed edges, indicator of model
output: expectation of inner product of vector of parameters and vector of statistics under hte assumption of independent edges
*/

double Lower_Bound(int n, int model, double *eta, double **mu, int directed);
/*
input: number of nodes, edge and triangle parameters, expectations of edges under the assumption of independent edges, indicator of directed edges, indicator of model
output: lower bound on log partition function of triangle ERGM obtained by variational approximation
*/

double Update_Expectations(int n, int model, int i, int j, double *eta, double **mu, int directed);
/*
input: number of nodes, nodes i and j, parameters, expectations of edges under the assumption of independent edges, indicator of directed edges
output: updated expectations of edges under the assumption of independent edges
*/

double EM(int n, int model, double *eta, int directed, int verbose);
/*
input: number of nodes, parameters, indicator of directed edges
output: lower bound of log partition function obtained by variational approximation
*/

double Variational_EM(int n, int model, double *eta, int directed, int verbose);
/*
input: number of nodes, parameters, indicator of directed edges
output: best lower bound on log partition function obtained by variational approximation
*/

