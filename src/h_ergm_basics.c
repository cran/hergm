/***************************************************************************/
/* Copyright 2009 Michael Schweinberger                                    */
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

#include "h_ergm_basics.h"

double** Lower_Half_Matrix(int n)
/* 
input: number of rows of matrix
output: pointer to lower half of n x n matrix
*/ 
{
  int i, j;
  double **x;
  x = (double**) calloc(n+1,sizeof(double*)); /* Allocate memory for n rows 0..n (row 0 is redundant) */
  if (x == NULL) 
    { 
    Rprintf("\n\ncalloc failed...\n\n"); 
    exit(1); 
    }
  for (i = 0; i < n + 1; i++) /* For row i, allocate memory for elements 0..i (element 0 is redundant) */
    {
    x[i] = (double*) calloc(i,sizeof(double));
    if (x[i] == NULL) 
      { 
      Rprintf("\n\ncalloc failed...\n\n"); 
      exit(1); 
      }
    }
  return x;
}

void Set_I_I(int d, int *vector1, int *vector2)
/*
input: two vectors
output: vector1 = vector2
*/
{
  int i;
  for (i = 0; i < d; i++)
    {
    vector1[i] = vector2[i];
    }
}

void Set_II_II(int d1, int d2, int **matrix1, int **matrix2)
/*
input: two matrices
output: matrix 1 = matrix2
*/
{
  int i, j;
  for (i = 0; i < d1; i++)
    {
    for (j = 0; j < d2; j++)
      {
      matrix1[i][j] = matrix2[i][j];
      }
    }
}

void Set_D_D(int d, double *vector1, double *vector2)
/*
input: two vectors
output: vector1 = vector2
*/
{
  int i;
  for (i = 0; i < d; i++)
    {
    vector1[i] = vector2[i];
    }
}

void Set_DD_DD(int d1, int d2, double **matrix1, double **matrix2)
/*
input: two matrices
output: matrix 1 = matrix2
*/
{
  int i, j;
  for (i = 0; i < d1; i++)
    {
    for (j = 0; j < d2; j++)
      {
      matrix1[i][j] = matrix2[i][j];
      }
    }
}

void Get_Column(int d, double *vector, double **matrix, int column)
/*
input: dimension of vector, matrix, column 
output: vector
*/ 
{
  int i;
  for (i = 0; i < d; i++)
    {
    vector[i] = matrix[i][column];
    }
}

void Set_Column(int d, double **matrix, int column, double *vector)
/*
input: dimension of vector, vector, column 
output: vector
*/ 
{
  int i;
  for (i = 0; i < d; i++)
    {
    matrix[i][column] = vector[i];
    }
}

double** Scale(int d1, int d2, double **matrix, double scale)
/* 
input: order of matrix, matrix, scale
output: scaled matrix
*/ 
{
  int i, j;
  double **x;
  x = (double**) calloc(d1,sizeof(double*));
  if (x == NULL) { Rprintf("\n\ncalloc failed: Scale, x\n\n"); exit(1); }
  for (i = 0; i < d1; i++)
    {
    x[i] = (double*) calloc(d2,sizeof(double));
    if (x[i] == NULL) { Rprintf("\n\ncalloc failed: Scale, x[%i]\n\n",i); exit(1); }
    }
  for (i = 0; i < d1; i++)
    {
    for (j = 0; j < d2; j++)
      {
      x[i][j] = matrix[i][j] * scale;
      }
    }
  return x;
}

void Print_I(int d, int *vector)
{
  int i;
  Rprintf("\n");
  for (i = 0; i < d; i++) Rprintf(" %i",vector[i]);
  Rprintf("\n");
}

void Print_II(int d1, int d2, int **matrix)
{
  int i, j;
  Rprintf("\n");
  for (i = 0; i < d1; i++)
    {
    for (j = 0; j < d2; j++) Rprintf(" %i",matrix[i][j]);
    Rprintf("\n");
    }
}

void Print_D(int d, double *vector)
{
  int i;
  Rprintf("\n");
  for (i = 0; i < d; i++) Rprintf(" %8.4f",vector[i]);
  Rprintf("\n");
}

void Print_DD(int d1, int d2, double **matrix)
{
  int i, j;
  Rprintf("\n");
  for (i = 0; i < d1; i++)
    {
    for (j = 0; j < d2; j++) Rprintf(" %8.4f",matrix[i][j]);
    Rprintf("\n");
    }
}

int Minimizer(int n, int *vector)
/*
input: vector of integers
output: minimum of vector of integers
*/
{
  int i, minimum, minimizer;
  minimizer = INT_MAX;
  minimum = INT_MAX;
  for (i = 0; i < n; i++)
    {
    if (vector[i] < minimum) 
      {
      minimizer = i;
      minimum = vector[i];
      }
    }
  return minimizer;
}
 
void Get_Permutation(long int n, long int index, int *permutation)
/*
input: number of elements, index of permutation, permutation
output: permutation which follows given permutation in lexigraphical order
*/
{
  int i, j, permutation_i, permutation_j;
  for (i = 0; i < n; i++) 
    {
    permutation[index+i] = permutation[index-n+i];
    }
  i = n - 1;
  while (permutation[index+i-1] > permutation[index+i]) i = i - 1; 
  j = n;
  while (permutation[index+j-1] < permutation[index+i-1]) j = j - 1; 
  /* Swap values at positions i-1 and j-1 */
  permutation_i = permutation[index+i-1];
  permutation_j = permutation[index+j-1];
  permutation[index+i-1] = permutation_j;
  permutation[index+j-1] = permutation_i;
  i++; 
  j = n;
  while (i < j)
    {
    permutation_i = permutation[index+i-1];
    permutation_j = permutation[index+j-1];
    permutation[index+i-1] = permutation_j;
    permutation[index+j-1] = permutation_i;
    i++;
    j--;
    }
}

void Permutations(long int *n_elements, long int *n_permutations, int *permutation)
/*
input: number of elements, number of permutations of elements, elements in increasing order stored in first row of permutations
output: all possible permutations of elements in lexigraphical order
*/
{
  long int i, index, n;
  n = (int)*n_elements;
  index = 0;
  for (i = 1; i < *n_permutations; i++)
    {
    index = index + n;
    Get_Permutation(n,index,permutation);
    }
}

