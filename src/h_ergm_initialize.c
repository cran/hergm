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

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/*

SEXP R_subColSummarize_avg_log(SEXP RMatrix, SEXP R_rowIndexList){

  static SEXP(*fun)(SEXP, SEXP) = NULL;
  
  if (fun == NULL)
    fun =  (SEXP(*)(SEXP, SEXP))R_GetCCallable("ergm","network_stats_wrapper");
  
  return fun(RMatrix, R_rowIndexList);
}

CHM_DN(*compute_statistic)(compute_statistic,SEXP);

compute_statistic = (CHM_DN(*)(CHM_DN,SEXP)) R_GetCCallable("ergm", "network_stats_wrapper");

*/
