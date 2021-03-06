# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

isTauSignificantlyChanged <- function(tauPrecision, tau, prevTau) {
    .Call('hergm_isTauSignificantlyChanged', PACKAGE = 'hergm', tauPrecision, tau, prevTau)
}

easy_E_step <- function(numOfVertices, numOfClasses, alpha, pi, stat, tau) {
    .Call('hergm_easy_E_step', PACKAGE = 'hergm', numOfVertices, numOfClasses, alpha, pi, stat, tau)
}

runFixedPointEstimationEStep <- function(numOfVertices, numOfClasses, alpha, pi, stat00, stat01, stat10, stat11, tau) {
    .Call('hergm_runFixedPointEstimationEStep', PACKAGE = 'hergm', numOfVertices, numOfClasses, alpha, pi, stat00, stat01, stat10, stat11, tau)
}

calculateStats <- function(network, stat00, stat01, stat10, stat11) {
    .Call('hergm_calculateStats', PACKAGE = 'hergm', network, stat00, stat01, stat10, stat11)
}

runFixedPointEstimationEStepMM <- function(numOfVertices, numOfClasses, alpha, pi, stat00, stat01, stat10, stat11, tau, network) {
    .Call('hergm_runFixedPointEstimationEStepMM', PACKAGE = 'hergm', numOfVertices, numOfClasses, alpha, pi, stat00, stat01, stat10, stat11, tau, network)
}

runModelEstimationMStep <- function(numOfVertices, numOfClasses, alpha, pi, stat00, stat01, stat10, stat11, tau) {
    .Call('hergm_runModelEstimationMStep', PACKAGE = 'hergm', numOfVertices, numOfClasses, alpha, pi, stat00, stat01, stat10, stat11, tau)
}

easy_M_Step <- function(numOfVertices, numOfClasses, alpha, pi, stat, tau) {
    .Call('hergm_easy_M_Step', PACKAGE = 'hergm', numOfVertices, numOfClasses, alpha, pi, stat, tau)
}