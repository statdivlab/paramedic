#' Extract posterior summaries from a hierarchical model fit
#'
#' Return point estimates and credible intervals for the true concentration, and point estimates and prediction intervals for estimated qPCR obtained through a Stan fit.
#'
#' @param stan_mod the model summary object from Stan.
#' @param stan_samps the list of MCMC samples from Stan.
#' @param taxa_of_interest the taxa to return estimates and posterior summaries for.
#' @param mult_num the number to multiply the resulting estimates and standard deviations by (defaults to 1).
#' @param level the \code{alpha} level for prediction intervals (defaults to 0.95, for a nominal 95\% prediction interval).
#' @param interval_type the type of prediction interval desired (defaults to "wald", but "quantile" is also acceptable).
#'
#' @return An object of class \code{paramedic}. See Details for more information
#'
#' @details A \code{paramedic} object is a list containing the following elements:
#' \itemize{
#'  \item{estimates}{ - the point estimates of qPCR (a matrix with dimension sample size by number of taxa).}
#'  \item{pred_intervals}{ - predction intervals for qPCR (an array with dimension sample size by 2 by number of taxa).}
#'  \item{est_efficiency}{ - point estimates for estimated varying efficiency, if varying efficiency was modeled (a vector of length number of taxa); otherwise, NA.}
#'  \item{efficiency_intervals}{ - posterior level \code{level}\eqn{\times}100\% confidence intervals for the true efficiency, if efficiency was modeled (a matrix of dimension number of taxa by 2); otherwise, NA.}
#' }
#'
#' @examples
#' ## load the package, read in example data
#' library("paramedic")
#' data(simple_example_data)
#'
#' ## process the example data
#' q <- 3
#' q_obs <- 2
#' processed_data <- process_data(full_data = simple_example_data, rel_inds = 1:q,
#' abs_inds = (q + 1):(q + 1 + q_obs),
#' abs_plus_rel_inds = 1:q_obs,
#' regex_thr = NA, regex_abs = "", llod = 0,
#' m_min = 1000, div_num = 1)
#'
#' ## run paramedic (with an extremely small number of iterations, for illustration only)
#' mod <- run_paramedic(W = processed_data$relative, V = processed_data$absolute,
#' stan_model = "src/stan_files/variable_efficiency.stan", n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747,
#' params_to_save = c("mu", "Sigma", "beta", "e"))
#'
#' ## get summary, samples
#' mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
#' mod_samps <- rstan::extract(mod)
#'
#' ## extract relevant summaries
#' summs <- extract_posterior_summaries(stan_mod = mod_summ, stan_samps = mod_samps, q = q, taxa_of_interest = 1:3,
#' mult_num = 1, level = 0.95, interval_type = "wald")
#'
#' @export
extract_posterior_summaries <- function(stan_mod, stan_samps, taxa_of_interest, mult_num = 1, level = 0.95, interval_type = "wald") {
  ## get the posterior estimates
  mu_summ_lst <- lapply(as.list(taxa_of_interest), function(x) stan_mod[grepl("mu", rownames(stan_mod)) & !grepl("log", rownames(stan_mod)) & grepl(paste0(",", x, "]"), rownames(stan_mod), fixed = TRUE), c(1, 3, 4, 5)]*mult_num)
  est_lst <- lapply(mu_summ_lst, function(x) x[, 1])
  sd_lst <- lapply(mu_summ_lst, function(x) x[, 2])
  ci_lst <- lapply(mu_summ_lst, function(x) x[, c(3, 4)])
  estimates <- do.call(cbind, est_lst)
  sd <- do.call(cbind, sd_lst)
  rownames(estimates) <- as.character(1:dim(estimates)[1])
  rownames(sd) <- as.character(1:dim(sd)[1])

  ## get cis
  credible_intervals <- sapply(1:length(mu_summ_lst), function(x) mu_summ_lst[[x]][, c(3, 4)], simplify = "array")
  rownames(credible_intervals) <- as.character(1:dim(estimates)[1])

  ## get prediction intervals
  if (interval_type == "wald") {
    intervals <- sapply(1:length(taxa_of_interest), function(x) gen_wald_interval(estimates[, x], sd[, x], alpha = 1 - level), simplify = "array")
  } else if (interval_type == "quantile") {
    intervals <- sapply(1:length(taxa_of_interest), function(x) gen_quantile_interval(mu_quantiles = credible_intervals[, , x], mu_samps = stan_samps$mu[, , taxa_of_interest[x]], div_num = mult_num, alpha = 1 - level, type = "credible_quantiles"), simplify = "array")
  } else { ## next, add quantile
    stop("Unsupported prediction interval type. Please enter one of 'quantile' or 'wald'.")
  }

  ## extract summaries of varying efficiency, if they exist
  if (any(grepl("e", rownames(stan_mod)) & !grepl("beta", rownames(stan_mod)))) {
    e <- stan_mod[grepl("e", rownames(stan_mod)) & !grepl("beta", rownames(stan_mod)), 1]
    e_intervals <- stan_mod[grepl("e", rownames(stan_mod)) &!grepl("beta", rownames(stan_mod)), c(4, 5)]
  } else {
    e <- NA
    e_intervals <- NA
  }

  return(list(estimates = estimates, pred_intervals = intervals, cred_intervals = credible_intervals, est_efficiency = e, efficiency_intervals = e_intervals, sd = sd))
}
