#include "h_ergm_basics.h"

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

