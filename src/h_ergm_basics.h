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

#include <stdio.h>
#include <R.h>
#include <math.h>
#include <Rmath.h>

double** Lower_Half_Matrix(int n);
/* 
input: number of rows of matrix
output: pointer to the upper half of n x n matrix
note: 
- number of elements is n (n - 1) / 2
- elements are 0
*/ 

void Set_I_I(int d, int *vector1, int *vector2);
/*
input: two vectors
output: vector1 = vector2
*/

void Set_II_II(int d1, int d2, int **matrix1, int **matrix2);
/*
input: two matrices
output: matrix 1 = matrix2
*/

void Set_D_D(int d, double *vector1, double *vector2);
/*
input: two vectors
output: vector1 = vector2
*/

void Set_DD_DD(int d1, int d2, double **matrix1, double **matrix2);
/*
input: two matrices
output: matrix 1 = matrix2
*/

void Get_Column(int d, double *vector, double **matrix, int column);
/*
input: dimension of vector, matrix, column 
output: vector
*/ 

void Set_Column(int d, double **matrix, int column, double *vector);
/*
input: dimension of vector, vector, column 
output: vector
*/ 

double** Scale(int d1, int d2, double **matrix, double scale);
/* 
input: order of matrix, matrix, scale
output: scaled matrix
*/ 

void Print_I(int d, int *vector);

void Print_II(int d1, int d2, int **matrix);

void Print_D(int d, double *vector);

void Print_DD(int d1, int d2, double **matrix);

void Get_Permutation(long int n, long int index, int *permutation);
/*
input: number of elements, index of permutation, permutation
output: permutation which follows given permutation in lexigraphical order
*/

int Max(int n, int *vector);
/*
input: vector of integers
output: maximum of vector of integers
*/

void Permutations(long int *n_elements, long int *n_permutations, int *permutation);
/*
input: number of elements, number of permutations of elements, elements in increasing order stored in first row of permutations
output: all possible permutations of elements in lexigraphical order
*/

