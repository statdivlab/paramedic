#' Run a Stan algorithm to estimate concentrations from combined qPCR and br16S data
#'
#' Estimate concentrations (and efficiencies, if applicable) using the specified Stan algorithm to combine qPCR and br16S data.
#'
#' @param W the br16S data.
#' @param V the qPCR data.
#' @param N the sample size.
#' @param q the total number of taxa.
#' @param q_obs the number of taxa with observed qPCR.
#' @param stan_model the Stan algorithm to fit to the data.
#' @param n_iter the total number of iterations per chain.
#' @param n_burnin the total number of warmups per chain.
#' @param n_chains the total number of chains.
#' @param stan_seed the random number seed to initialize.
#' @param params_to_save a character vector of the parameters to save.
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
#' processed_data <- process_data(full_data = simple_example_data, br_inds = 1:q,
#' qpcr_inds = (q + 1):(q + 1 + q_obs),
#' pcr_plus_br_inds = 1:q_obs,
#' regex_thr = NA, regex_cps = "", llod = 0,
#' m_min = 1000, div_num = 1000)
#'
#' ## run paramedic (with an extremely small number of iterations, for illustration only)
#' mod <- paramedic(W = processed_data$br, V = processed_data$qpcr, q = q, q_obs = q_obs,
#' stan_model = "src/stan_files/variable_efficiency.stan", n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747,
#' params_to_save = c("mu", "Sigma", "beta", "e"))
#'
#' @seealso \code{\link[rstan]{stan}} for specific usage of the \code{stan} function.
#'
#' @export
paramedic <- function(W, V, N = dim(W)[1], q = dim(W)[2], q_obs = dim(V)[2],
                      stan_model = "src/stan_files/variable_efficiency.stan",
                      n_iter = 10500, n_burnin = 10000, n_chains = 4, stan_seed = 4747,
                      params_to_save = c("mu", "Sigma", "beta", "e"),
                      inits_lst = NULL,
                      ...) {
    ## error messages
    if (dim(W)[1] != dim(V)[1]) stop("The number of rows in W and V must match.")
    if (dim(W)[1] != N) stop("The number of rows in W and V must be equal to the sample size N.")
    if (dim(W)[2] != q) stop("The number of columns in W must be equal to the total number of taxa q.")
    if (dim(V)[2] != q_obs) stop("The number of columns in V must be equal to the number of taxa with observed qPCR, q_obs.")
    ## set up the data list
    data_lst <- list(W = W, V = V, N = N, q = q, q_obs = q_obs)
    ## get inits from the naive estimator
    naive_estimator <- function(idx, brs, qpcrs, known_qpcr) {
      ## get the sum of the br16s's
      sum_br16s <- rowSums(brs[, known_qpcr, drop = FALSE])
      ## get the sum of the qpcr's
      sum_qpcr <- rowSums(qpcrs)

      ## naive estimator
      qpcr <- brs[, idx]*sum_qpcr/sum_br16s
      ## if sum_br16s was zero, predict zero (?)
      qpcr <- ifelse (sum_br16s == 0, 0, qpcr)
      return(qpcr)
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
