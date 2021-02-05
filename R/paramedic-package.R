#' The 'paramedic' package.
#'
#' @description We estimate microbial abundances by combining compositional data and observed absolute abundance data and modeling efficiencies.
#' We use Stan to fit the chosen hierarchical model to the data: this may result in differences in speed and accuracy between the hierarchical models, depending on the dataset and initial values chosen. 
#' Please read any warning messages returned by the algorithm carefully, as these can help diagnose convergence issues. 
#' We return credible intervals and point estimates for the true concentrations, and point estimates and prediction intervals for the unobserved absolute abundances. 
#'
#' @section Author(s):
#' \bold{Maintainer}: Brian Williamson \url{http://bdwilliamson.github.io}
#' 
#' Methodology authors:
#' \itemize{
#'   \item{Brian D. Williamson}
#'   \item{James P. Hughes}
#'   \item{Amy D. Willis}
#' }
#' 
#' @section See also:
#' The preprint: \url{https://www.biorxiv.org/content/10.1101/761486v1}
#' Other useful links:
#' \itemize{
#'   \item{\url{http://bdwilliamson.github.io/paramedic}}
#'   \item{\url{http://github.com/statdivlab/paramedic}}
#'   \item{Report bugs at \url{http://github.com/statdivlab/paramedic/issues}}
#' }
#' 
#' @docType package
#' @name paramedic-package
#' @aliases paramedic
#' @useDynLib paramedic, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling extract summary
#' @importFrom magrittr %>%
#' @importFrom stats rnorm rmultinorm rnbinom rpois
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.17.4. http://mc-stan.org
#'
NULL
