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
