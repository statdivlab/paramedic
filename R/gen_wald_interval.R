#' Generate a Wald-type prediction interval for qPCR
#'
#' Return a Wald-type prediction interval for qPCR values given the estimated concentrations.
#'
#' @param mu the estimated concentrations (a vector of length N times q).
#' @param sd the estimated standard deviation of the estimated concentrations (a vector of length N times q).
#' @param alpha the desired level (defaults to 0.05, corresponding to a 95\% interval).
#' @param truncate truncate negative lower limits at zero (defaults to \code{TRUE})
#'
#' @return A (1 - \eqn{\alpha})x100\% Wald-type prediction interval for each qPCR
#'
#' @examples
#' ## load the package, read in example data
#' library("paramedic")
#' data(example_16S_data)
#' data(example_qPCR_data)
#'
#' ## run paramedic (with an extremely small number of iterations, for illustration only)
#' mod <- run_paramedic(W = example_16S_data, V = example_qPCR_data,
#' stan_model = "src/stan_files/variable_efficiency.stan", n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747,
#' params_to_save = c("mu", "Sigma", "beta", "e"))
#' ## get model summary
#' mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
#' ## get samples
#' mod_samps <- rstan::extract(mod)
#' ## extract relevant summaries
#' summs <- extract_posterior_summaries(stan_mod = mod_summ, stan_samps = mod_samps, taxa_of_interest = 1:3,
#' mult_num = 1, level = 0.95, interval_type = "wald")
#'
#' @export
gen_wald_interval <- function(mu, sd, alpha = 0.05, truncate = TRUE) {
  ## check if it is a matrix first; if so, make it a vector
  if (is.matrix(mu)) mu <- as.vector(t(mu))
  if (is.matrix(sd)) sd <- as.vector(t(sd))
  ## estimated variance
  est_var <- mu + sd^2

  ## get intervals
  pred_intervals <- mu + sqrt(est_var) %o% qnorm(c(alpha/2, 1 - alpha/2))

  ## truncate
  if (truncate) {
    pred_intervals[, 1] <- ifelse(pred_intervals[, 1] < 0, 0, pred_intervals[, 1])
  }

  return(pred_intervals)
}
