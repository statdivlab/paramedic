#' Run a Stan algorithm to estimate concentrations from combined absolute and relative abundance data
#'
#' Estimate concentrations (and efficiencies, if applicable) using the specified Stan algorithm to combine absolute and relative abundance data.
#'
#' @param W The relative abundance data, e.g., from broad range 16S sequencing with "universal" primers.
#' @param V The absolute abundance data, e.g., from taxon-specific absolute primers.
#' @param stan_model The Stan algorithm to fit to the data. Expects an R object of class `stanmodel` corresponding to a pre-compiled Stan file (either \code{stanmodels$variable_efficiency} or \code{stanmodels$variable_efficiency_centered}).
#' @param n_iter The total number of iterations per chain to be run by the Stan algorithm. Defaults to 10500.
#' @param n_burnin The total number of warmups per chain to be run by the Stan algorithm. Defaults to 10000.
#' @param n_chains The total number of chains to be run by the Stan algorithm. Defaults to 4.
#' @param stan_seed The random number seed to initialize.
#' @param inits_lst An optional list of initial values of the parameters. Must be a named list; see \code{\link[rstan]{stan}}.
#' @param ... other arguments to pass to \code{\link[rstan]{sampling}} (e.g., control).
#'
#' @return An object of class \code{stanfit}.
#'
#' @details We fit a hierarchical model in Stan to the data, with goal to estimate true concentration for all taxa. There are two available hierarchical models. The two only differ in the way that the true concentration is parameterized. While these hierarchical models are mathematically identical, they are fit differently by Stan. The first hierarchical model is a "centered" model, where \deqn{\beta ~ N(0, \sqrt{50})} \deqn{\log \Sigma ~ N(0, \sqrt{50})} \deqn{\log \mu ~ N(\beta, \Sigma).} We call this model \code{variable_efficiency_centered.stan}. The second hierarchical model is a "noncentered" model, where \deqn{\beta ~ N(0, \sqrt{50})} \deqn{\log \Sigma ~ N(0, \sqrt{50})} \deqn{\log \gamma ~ N(0, I)} \deqn{\log \mu = \sqrt{\Sigma} \log \gamma + \beta.} We call this model \code{variable_efficiency.stan}.

#' In most cases, we suggest using the noncentered model. However, if the model does not converge in a reasonable amount of time using the noncentered model, consider trying the centered model.
#'
#' @examples
#' ## load the package, read in example data
#' library("paramedic")
#' data(example_16S_data)
#' data(example_qPCR_data)
#'
#' ## run paramedic (with an extremely small number of iterations, for illustration only)
#' ## on only the first 10 taxa
#' mod <- run_paramedic(W = example_16S_data[, 1:10], V = example_qPCR_data,
#' stan_model = stanmodels$variable_efficiency, n_iter = 30, n_burnin = 25, 
#' n_chains = 1, stan_seed = 4747)
#'
#' @seealso \code{\link[rstan]{stan}} and \code{\link[rstan]{sampling}} for specific usage of the \code{stan} and \code{sampling} functions.
#'
#' @export
run_paramedic <- function(W, V,
                      stan_model = stanmodels$variable_efficiency,
                      n_iter = 10500, n_burnin = 10000, n_chains = 4, stan_seed = 4747,
                      inits_lst = NULL,
                      ...) {
    N <- dim(W)[1]
    q <- dim(W)[2]
    q_obs <- dim(V)[2]
    ## error messages
    if (dim(W)[1] != dim(V)[1]) stop("The number of rows in W and V must match.")
    ## set up the data list
    data_lst <- list(W = W, V = V, N = N, q = q, q_obs = q_obs)
    ## get inits from the naive estimator
    naive_estimator <- function(idx, relatives, absolutes, known_absolute) {
      ## get the sum of the relative abundances
      sum_relative <- rowSums(relatives[, known_absolute, drop = FALSE])
      ## get the sum of the absolutes
      sum_absolute <- rowSums(absolutes)

      ## naive estimator
      absolute <- relatives[, idx]*sum_absolute/sum_relative
      ## if sum_relative was zero, predict zero (?)
      absolute <- ifelse (sum_relative == 0, 0, absolute)
      return(absolute)
    }
    naive_est <- cbind(as.matrix(V), apply(matrix((q_obs + 1):q), 1, naive_estimator, W, V, 1:q_obs))
    colnames(naive_est) <- colnames(W)
    log_naive <- ifelse(is.infinite(log(naive_est)), 0, log(naive_est))
    naive_beta <- colMeans(log_naive, na.rm = TRUE)
    naive_Sigma <- diag(stats::var(log_naive, na.rm = TRUE))
    tmp <- ifelse(naive_Sigma < 1e-3, 1e-3, naive_Sigma)
    naive_Sigma <- tmp
    log_naive_tilde <- sweep(sweep(log_naive, 2, naive_beta, FUN = "-"), 2, naive_Sigma, FUN = "/")
    tmp <- ifelse(is.na(log_naive_tilde), 0, log_naive_tilde)
    log_naive_tilde <- tmp
    if (is.null(inits_lst)) { # create inits if not passed in
        if (n_chains > 1) {
            if (grepl("centered", stan_model@model_name)) {
                inits_lst <- list(list(log_mu = log_naive), rep(list(init = "random"), n_chains - 1))
            } else {
                inits_lst <- list(list(log_mu_tilde = log_naive_tilde), rep(list(init = "random"), n_chains - 1))
            }
        } else {
            if (grepl("centered", stan_model@model_name)) {
                inits_lst <- list(list(log_mu = log_naive, beta = naive_beta, Sigma = naive_Sigma))
            } else {
                inits_lst <- list(list(log_mu_tilde = log_naive_tilde, beta = naive_beta, Sigma = naive_Sigma))
            }
        }
    }

    ## run the Stan algorithm
    mod <- rstan::sampling(stan_model, data = data_lst, pars = c("mu", "e", "beta", "Sigma"),
                           chains = n_chains, iter = n_iter, warmup = n_burnin, seed = stan_seed,
                           init = inits_lst, ...)
    return(mod)
}
