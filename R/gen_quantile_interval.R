#' Generate a Poisson-quantile-based prediction interval for qPCR
#'
#' Return a Poisson-quantile-based prediction interval for qPCR values given Markov Chain Monte Carlo samples for the estimated concentrations.
#'
#' @param mu_quantiles the (alpha/2, 1 - alpha/2) quantiles from the MCMC sampling distribution for the true absolute concentration \eqn{\mu}; only supply if type = "credible_quantiles" (a matrix of dimension N (the sample size) x q (the number of taxa) by 2).
#' @param mu_samps the estimated concentrations [an array with dimension (number of MCMC samples) by N by q]; only supply if type = "sample_quantiles".
#' @param alpha the desired level (defaults to 0.05, corresponding to an interval using the 2.5\% and 97.5\% quantiles)
#' @param type the type of intervals desired, either "credible_quantiles" or "sample_quantiles" (please see Details for more information on the difference between these two types).
#' @param div_num the number to multiply by.
#'
#' @return A (1 - \eqn{\alpha})x100\% Poisson-quantile-based prediction interval for each qPCR
#'
#' @examples
#' # load the package, read in example data
#' library("paramedic")
#' data(example_16S_data)
#' data(example_qPCR_data)
#'
#' # run paramedic (with an extremely small number of iterations, for illustration only)
#' # on only the first 10 taxa
#' mod <- run_paramedic(W = example_16S_data[, 1:10], V = example_qPCR_data,
#' n_iter = 30, n_burnin = 25, 
#' n_chains = 1, stan_seed = 4747)
#' # get model summary
#' mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
#' # get samples
#' mod_samps <- rstan::extract(mod$stan_fit)
#' # extract relevant summaries
#' summs <- extract_posterior_summaries(stan_mod = mod_summ, stan_samps = mod_samps, 
#' taxa_of_interest = 1:3,
#' mult_num = 1, level = 0.95, interval_type = "quantile")
#'
#' @export
gen_quantile_interval <- function(mu_quantiles, mu_samps, alpha = 0.05, type = "credible_quantiles", div_num = 1) {
  if (type == "credible_quantiles") {
    if (is.null(dim(mu_quantiles))) mu_quantiles <- matrix(mu_quantiles, ncol = 2)

    lower_limits <- stats::qpois(p = alpha/2, lambda = mu_quantiles[, 1])
    upper_limits <- stats::qpois(p = 1 - alpha/2, lambda = mu_quantiles[, 2])
    pred_intervals <- cbind(lower_limits, upper_limits)
  } else if (type == "sample_quantiles") {
    # check to make sure it's an array of dimension 3
    if (length(dim(mu_samps)) != 3) mu_samps <- array(mu_samps, dim = c(dim(mu_samps), 1))

    # sample from Poisson
    v_samps <- apply(mu_samps, c(1, 2, 3), function(x) stats::rpois(1, x))

    # compute quantiles for each (i,j) pair
    quantiles <- apply(v_samps, c(2, 3), function(x) stats::quantile(x, probs = c(alpha/2, 1 - alpha/2)))

    pred_intervals <- matrix(quantiles, nrow = dim(quantiles)[2]*dim(quantiles)[3], ncol = 2, byrow = TRUE)
  }

  return(pred_intervals)
}
