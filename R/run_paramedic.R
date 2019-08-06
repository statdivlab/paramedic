absolute#' Run a Stan algorithm to estimate concentrations from combined absolute and relative abundance data
#'
#' Estimate concentrations (and efficiencies, if applicable) using the specified Stan algorithm to combine absolute and relative abundance data.
#'
#' @param W The relative abundance data, e.g., from broad range 16S sequencing with "universal" primers.
#' @param V The absolute abundance data, e.g., from taxon-specific absolute primers.
#' @param N Optional. The number of samples. Defaults to the row-length of W.
#' @param q Optional. The number of taxa observed with the relative abundance technology. Defaults to the column-length of W.
#' @param q_obs Optional. The number of taxa observed with the absolute abundance technology. Defaults to the column-length of V.
#' @param stan_model The Stan algorithm to fit to the data. Expects a file path.
#' @param n_iter The total number of iterations per chain to be run by the Stan algorithm. Defaults to 10500.
#' @param n_burnin The total number of warmups per chain to be run by the Stan algorithm. Defaults to 10000.
#' @param n_chains The total number of chains to be run by the Stan algorithm. Defaults to 4.
#' @param stan_seed The random number seed to initialize.
#' @param params_to_save A character vector of the parameters to save. Defaults to "mu", "Sigma", "beta", "e". TODO(Amy): clarify.
#' @param ... other arguments to pass to \code{\link[rstan]{stan}}.
#'
#' @return An object of class \code{stanfit}.
#'
#' @details We fit a hierarchical model in Stan to the data, with goal to estimate true concentration for all taxa. There are two available hierarchical models. The two only differ in the way that the true concentration is parameterized. While these hierarchical models are mathematically identical, they are fit differently by Stan. The first hierarchical model is a "centered" model, where \deqn{\beta ~ N(0, \sqrt{50})} \deqn{\log \Sigma ~ N(0, \sqrt{50})} \deqn{\log \mu ~ N(\beta, \Sigma).} We call this model \code{variable_efficiency_centered.stan}. The second hierarchical model is a "noncentered" model, where \deqn{\beta ~ N(0, \sqrt{50})} \deqn{\log \Sigma ~ N(0, \sqrt{50})} \deqn{\log \gamma ~ N(0, I)} \deqn{\log \mu = \sqrt{\Sigma} \log \gamma + \beta.} We call this model \code{variable_efficiency.stan}.
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
#' m_min = 1000, div_num = 1000)
#'
#' ## run paramedic (with an extremely small number of iterations, for illustration only)
#' mod <- run_paramedic(W = processed_data$br, V = processed_data$absolute,
#' stan_model = "src/stan_files/variable_efficiency.stan", n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747,
#' params_to_save = c("mu", "Sigma", "beta", "e"))
#'
#' @seealso \code{\link[rstan]{stan}} for specific usage of the \code{stan} function.
#'
#' @export
run_paramedic <- function(W, V,
                      stan_model = "src/stan_files/variable_efficiency.stan",
                      n_iter = 10500, n_burnin = 10000, n_chains = 4, stan_seed = 4747,
                      params_to_save = c("mu", "Sigma", "beta", "e"),
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
    naive_Sigma <- diag(var(log_naive, na.rm = TRUE))
    tmp <- ifelse(naive_Sigma < 1e-3, 1e-3, naive_Sigma)
    naive_Sigma <- tmp
    log_naive_tilde <- sweep(sweep(log_naive, 2, naive_beta, FUN = "-"), 2, naive_Sigma, FUN = "/")
    tmp <- ifelse(is.na(log_naive_tilde), 0, log_naive_tilde)
    log_naive_tilde <- tmp
    if (is.null(inits_lst)) { # create inits if not passed in
        if (n_chains > 1) {
            if (grepl("centered", stan_model)) {
                inits_lst <- list(list(log_mu = log_naive), rep(list(init = "random"), n_chains - 1))
            } else {
                inits_lst <- list(list(log_mu_tilde = log_naive_tilde), rep(list(init = "random"), n_chains - 1))
            }
        } else {
            if (grepl("centered", stan_model)) {
                inits_lst <- list(list(log_mu = log_naive, beta = naive_beta, Sigma = naive_Sigma))
            } else {
                inits_lst <- list(list(log_mu_tilde = log_naive_tilde, beta = naive_beta, Sigma = naive_Sigma))
            }
        }
    }

    ## run the Stan algorithm
    mod <- rstan::stan(file = stan_model, data = data_lst, iter = n_iter,
    warmup = n_burnin, chains = n_chains, seed = stan_seed,
    pars = params_to_save, init = inits_lst, ...)
    return(mod)
}