# Make C/C++ functions of R package hergm available to C/C++ functions of other R packages by registering the functions

compute_statistic = R_GetCCallable("ergm", "network_stats_wrapper");

R_RegisterCCallable("hergm", "P_Independence", P_Independence);
DL_FUNC R_GetCCallable("hergm", "P_Independence");

R_RegisterCCallable("hergm", "Inference", Inference);
DL_FUNC R_GetCCallable("hergm", "Inference");

