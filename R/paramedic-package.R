#' The 'paramedic' package.
#'
#' @description We estimate microbial abundances by combining compositional 16S data and observed absolute abundance obtained using quantitative PCR. We have two fully Bayesian hierarchical models allowing for varying efficiencies: the first uses a centered parameterization, while the second uses a noncentered parameterization. We use Stan to fit the chosen hierarchical model to the data: this may result in differences in speed and accuracy between the hierarchical models, depending on the dataset and initial values chosen. Please read any warning messages returned by the algorithm carefully, as these can help diagnose convergence issues. We return credible intervals and point estimates for the true concentrations, and point estimates and prediction intervals for the unobserved qPCR abundances. 
#'
#' @docType package
#' @name paramedic-package
#' @aliases paramedic
#' @useDynLib paramedic, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.17.4. http://mc-stan.org
#'
NULL
