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

#include "h_ergm_variational.h"

double Expected_Density(int n, double **mu, int directed)
/*
input: number of nodes, expectations of edges under the assumption of independent edges, indicator of directed edges
output: expectation of number of edges under the assumption of independent edges
*/
{
  int i, j;
  double sum;
  sum = 0.0;
  for (i = 0; i < n - 1; i++)
    {
    for (j = i + 1; j < n; j++)
      {
      sum = sum + mu[i][j];
      if (directed == 1) sum = sum + mu[j][i];
      }
    }
  return sum;
}

double Expected_Transitivity(int n, double **mu, int directed)
/*
input: number of nodes, expectations of edges under the assumption of independent edges, indicator of directed edges
output: expected number of transitive triples under the assumption of independent edges
*/
{
  int i, j, k;
  double sum;
  sum = 0.0;
  for (i = 0; i < n - 2; i++)
    {
    for (j = i + 1; j < n - 1; j++)
      {
      for (k = j + 1; k < n; k++)
        {
        sum = sum + (mu[i][j] * mu[j][k] * mu[i][k]); /* i, j, k */
        if (directed == 1) sum = sum + (mu[i][k] * mu[k][j] * mu[i][j]) /* i, k, j */
                                     + (mu[j][i] * mu[i][k] * mu[j][k]) /* j, i, k */
                                     + (mu[j][k] * mu[k][i] * mu[j][i]) /* j, k, i */
                                     + (mu[k][i] * mu[i][j] * mu[k][j]) /* k, i, j */
                                     + (mu[k][j] * mu[j][i] * mu[k][i]); /* k, j, i */
        }
      }
    }
  return sum;
}

double D_Expected_Transitivity(int n, int i, int j, double **mu, int directed)
/*
input: number nodes, nodes i and j, expectations of edges under the assumption of independent edges, indicator of directed edges
output: expectation of number of two-paths between nodes i and j under the assumption of independent edges
*/
{
  int k;
  double sum;
  sum = 0.0;
  for (k = 0; k < n; k++) /* Implicit assumption: expectations of self-edges must be 0 */
    {
    if ((k != i) && (k != j))
      {
      sum = sum + (mu[i][k] * mu[j][k]); /* First possible edge */
      if (directed == 1) sum = sum + (mu[k][i] * mu[k][j]) /* Second possible edge */
                                   + (mu[i][k] * mu[k][j]); /* Third possible edge */
      } 
    }
  return sum;
}

double Expected_Stars(int n, double **mu)
/*
input: number of nodes, expectations of edges under the assumption of independent edges, indicator of directed edges
output: expected number of stars under the assumption of independent edges
*/
{
  int i, j, k;
  double sum;
  sum = 0.0;
  for (i = 0; i < n; i++) /* i is the center of n - 1 2-stars */
    {
    for (j = 0; j < n - 1; j++)
      {
      for (k = j + 1; k < n; k++)
        {
        if ((k != i) && (k != j)) sum = sum + (mu[i][j] * mu[i][k]);
        }
      }
    }
  return sum;
}

double D_Expected_Stars(int n, int i, int j, double **mu)
/*
input: number nodes, nodes i and j, expectations of edges under the assumption of independent edges, indicator of directed edges
output: expectation of number of two-paths between nodes i and j under the assumption of independent edges
*/
{
  int k;
  double sum;
  sum = 0.0;
  for (k = 0; k < n; k++) /* Implicit assumption: expectations of self-edges must be 0 */
    {
    if ((k != i) && (k != j)) sum = sum + mu[i][k] + mu[j][k]; 
    }
  return sum;
}

double Entropy(int n, double **mu, int directed)
/*
input: number of nodes, expectations of edges under the assumption of independence, indicator of directed edges
output: entropy of distribution under the assumption of independent edges
*/
{
  int i, j;
  double entropy;
  entropy = 0.0;
  for (i = 0; i < n - 1; i++)
    {
    for (j = i + 1; j < n; j++)
      {
      entropy = entropy - (mu[i][j] * ln(mu[i][j])) - ((1.0 - mu[i][j]) * ln(1.0 - mu[i][j]));
      if (directed == 1) entropy = entropy - (mu[j][i] * ln(mu[j][i])) - ((1.0 - mu[j][i]) * ln(1.0 - mu[j][i]));
      }
    }
  return entropy;
}

double Expected_Energy(int n, int model, double *eta, double **mu, int directed)
/*
input: number of nodes, edge and triangle parameters, expectations of edges under the assumption of independent edges, indicator of directed edges, indicator of model
output: expectation of inner product of vector of parameters and vector of statistics under hte assumption of independent edges
*/
{
  double expected_energy, expected_statistic_0, expected_statistic_1;
  expected_statistic_0 = Expected_Density(n,mu,directed);
  expected_statistic_1 = 0.0;
  if (model == 1) expected_statistic_1 = Expected_Stars(n,mu);
  else if (model == 2) expected_statistic_1 = Expected_Transitivity(n,mu,directed);
  expected_energy = (eta[0] * expected_statistic_0) + (eta[1] * expected_statistic_1);
  return expected_energy;
}

double Lower_Bound(int n, int model, double *eta, double **mu, int directed)
/*
input: number of nodes, edge and triangle parameters, expectations of edges under the assumption of independent edges, indicator of directed edges, indicator of model
output: lower bound on log partition function of triangle ERGM obtained by variational approximation
*/
{
  double expected_energy, entropy, lower_bound;
  expected_energy = Expected_Energy(n,model,eta,mu,directed);
  entropy = Entropy(n,mu,directed);
  lower_bound = expected_energy + entropy;
  return lower_bound;
}

double Update_Expectations(int n, int model, int i, int j, double *eta, double **mu, int directed)
/*
input: number of nodes, nodes i and j, parameters, expectations of edges under the assumption of independent edges, indicator of directed edges
output: updated expectations of edges under the assumption of independent edges
*/
{
  double d_expected_statistic, log_odds, expectation;
  d_expected_statistic = 0.0;
  if (model == 1) d_expected_statistic = D_Expected_Stars(n,i,j,mu);
  else if (model == 2) d_expected_statistic = D_Expected_Transitivity(n,i,j,mu,directed); 
  log_odds = eta[0] + (eta[1] * d_expected_statistic);
  expectation = 1.0 / (1.0 + e(-log_odds));
  return expectation;
}

double EM(int n, int model, double *eta, int directed, int verbose)
/*
input: number of nodes, parameters, indicator of directed edges
output: lower bound of log partition function obtained by variational approximation
*/ 
{
  int i, j, iteration;
  double last_lower_bound, lower_bound, **mu;
  mu = (double**) calloc(n,sizeof(double*)); 
  if (mu == NULL) { Rprintf("\n\nEM: calloc failed...\n\n"); error("Error: out of memory"); }
  for (i = 0; i < n; i++) 
    {
    mu[i] = (double*) calloc(n,sizeof(double));
    if (mu[i] == NULL) { Rprintf("\n\nEM: calloc failed...\n\n"); error("Error: out of memory"); }
    }
  for (i = 0; i < n - 1; i++) /* Main diagonal of mu is 0, since it is initialized as 0 */
    {
    for (j = i + 1; j < n; j++)
      {
      mu[i][j] = unif_rand();
      if (directed == 0) mu[j][i] = mu[i][j];
      else mu[j][i] = unif_rand();
      }
    }
  if (verbose == 1) 
    {
    Rprintf("\nVariational EM:\n");
    Rprintf("\niteration   lower bound A(eta)");
    Rprintf("\n------------------------------");
    }
  lower_bound = -DBL_MAX;
  iteration = 0;
  do /* Experiments suggests that in most cases 3 iterations suffice */
    {
    iteration = iteration + 1;
    for (i = 0; i < n - 1; i++) /* Main diagonal of mu remains 0, since it is initialized at 0 and not updated */
      {
      for (j = i + 1; j < n; j++)
        {
        mu[i][j] = Update_Expectations(n,model,i,j,eta,mu,directed);
        if (directed == 0) mu[j][i] = mu[i][j];
        else mu[j][i] = Update_Expectations(n,model,j,i,eta,mu,directed);
        }
      }
    last_lower_bound = lower_bound;
    lower_bound = Lower_Bound(n,model,eta,mu,directed);
    if (verbose == 1) Rprintf("\n%i %8.4f",iteration,lower_bound);
    }
  while ((lower_bound - last_lower_bound) > 0.01);
  if (verbose == 1) Rprintf("\n------------------------------\n");
  for (i = 0; i < n; i++)
    {
    free(mu[i]);
    }
  free(mu);
  return lower_bound;
}

double Variational_EM(int n, int model, double *eta, int directed, int verbose)
/*
input: number of nodes, parameters, indicator of directed edges
output: best lower bound on log partition function obtained by variational approximation
*/
{
  int run;
  double lower_bound, best_lower_bound;
  best_lower_bound = -DBL_MAX;
  for (run = 0; run < 5; run++)
    {
    lower_bound = EM(n,model,eta,directed,verbose);
    if (lower_bound > best_lower_bound) best_lower_bound = lower_bound;
    }
  return best_lower_bound;
}

