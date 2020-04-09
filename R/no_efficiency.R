#' Run a Stan algorithm to estimate concentrations from combined absolute and relative abundance data without modeling efficiency
#'
#' Estimate concentrations (without modeling efficiency) by combining absolute and relative abundance data.
#'
#' @param W The relative abundance data, e.g., from broad range 16S sequencing with "universal" primers. Expects data (e.g., matrix, data.frame, tibble) with sample identifiers in the first column. Sample identifiers must be the same between W and V, and the column must have the same name in W and V.
#' @param V The absolute abundance data, e.g., from taxon-specific absolute primers. Expects data (e.g., matrix, data.frame, tibble) with sample identifiers in the first column. Sample identifiers must be the same between W and V, and the column must have the same name in W and V.
#' @param X The covariate data. Expects data (e.g., matrix, data.frame, tibble) with sample identifiers in the first column. Sample identifiers must be the same between W, V, and X, and the column must have the same name in W, V, and X. If X only consists of the subject identifiers, then no covariates are used.
#' @param n_iter The total number of iterations per chain to be run by the Stan algorithm. Defaults to 10500.
#' @param n_burnin The total number of warmups per chain to be run by the Stan algorithm. Defaults to 10000.
#' @param n_chains The total number of chains to be run by the Stan algorithm. Defaults to 4.
#' @param stan_seed The random number seed to initialize.
#' @param inits_lst An optional list of initial values of the parameters. Must be a named list; see \code{\link[rstan]{stan}}.
#' @param sigma_beta Hyperparameter specifying the prior variance on \eqn{\beta_0}. Defaults to \eqn{\sqrt{50}}.
#' @param sigma_Sigma Hyperparameter specifying the prior variance on \eqn{\Sigma}. Defaults to \eqn{\sqrt{50}}.
#' @param ... other arguments to pass to \code{\link[rstan]{sampling}} (e.g., control).
#'
#' @return An object of class \code{stanfit}.
#'
#' @details We fit a hierarchical model in Stan to the data, with goal to estimate true concentration for all taxa. This model, in contrast with \code{run_paramedic} and \code{run_paramedic_centered}, does not allow for varying efficiency. **We always recommend using the model that allows for varying efficiency**, unless it is known *a priori* that there is not varying efficiency in the relative abundance data. We include this model to facilitate replicating the results of the manuscript.

#'
#' @examples
#' ## load the package, read in example data
#' library("paramedic")
#' data(example_16S_data)
#' data(example_qPCR_data)
#'
#' ## run no efficiency (with an extremely small number of iterations, for illustration only)
#' ## on only the first 10 taxa
#' mod <- no_efficiency(W = example_16S_data[, 1:10], V = example_qPCR_data,
#' n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747)
#'
#' @seealso \code{\link[rstan]{stan}} and \code{\link[rstan]{sampling}} for specific usage of the \code{stan} and \code{sampling} functions.
#'
#' @export
no_efficiency <- function(W, V, X = V[, 1, drop = FALSE],
                      n_iter = 10500, n_burnin = 10000, n_chains = 4, stan_seed = 4747,
                      inits_lst = NULL, sigma_beta = sqrt(50), sigma_Sigma = sqrt(50),
                      ...) {
    N <- dim(W)[1]
    q <- dim(W)[2] - 1
    q_obs <- dim(V)[2] - 1
    p <- dim(X)[2]
    ## --------------
    ## error messages
    ## --------------
    ## error if there aren't the same number of rows in W and V
    if (dim(W)[1] != dim(V)[1]) stop("The number of rows in W and V must match.")
    ## error if there aren't the same number of rows in V and X
    if (dim(V)[1] != dim(X)[1]) stop("The number of rows in V and X must match.")
    ## error if there is any difference between the first column of W and V
    if (length(setdiff(as.numeric(unlist(W[, 1])), as.numeric(unlist(V[, 1])))) != 0) stop("W and V must have the same samples.")
    ## error if there is any difference between the first column of X and V
    if (length(setdiff(as.numeric(unlist(V[, 1])), as.numeric(unlist(X[, 1])))) != 0) stop("V and X must have the same samples.")
    ## error if the id columns are named different things
    if ("data.frame" %in% class(W) | "tbl" %in% class(W)) {
      if (names(W)[1] != names(V)[1]) stop("W and V must have the same name.")
      if (names(W)[1] != names(X)[1]) stop("W and X must have the same name.")
    } else {
      if (colnames(W)[1] != colnames(V)[1]) stop("W and V must have the same name.")
      if (colnames(W)[1] != colnames(X)[1]) stop("W and X must have the same name.")
    }
    ## error if q < q_obs
    if (q < q_obs) stop("V must have fewer taxa than W (or the same number of taxa).")

    ## ---------------------------
    ## pre-processing and warnings
    ## ---------------------------
    ## make tibbles out of W and V, if they aren't already
    W_tbl <- tibble::as_tibble(W)
    V_tbl <- tibble::as_tibble(V)
    X_tbl <- tibble::as_tibble(X)
    ## if the rows are scrambled between W and V, change W to match V
    if (sum(W_tbl[, 1] == V_tbl[, 1]) != dim(V_tbl)[1]) {
      warning("Re-ordering samples so that the rows of W match the rows of V. The results will be in terms of the rows of V.")
      combined_tbl <- dplyr::left_join(V_tbl, W_tbl, by = names(V_tbl)[1])
      tmp <- combined_tbl %>%
        dplyr::select(-dplyr::ends_with(".x")) %>%
        dplyr::rename_at(.vars = dplyr::vars(dplyr::ends_with(".y")),
                  .funs = list(~sub("[.]y$", "", .)))
      W_tbl <- tmp
    }
    ## if the rows are scrambled between X and V, change X to match V
    if (sum(X_tbl[, 1] == V_tbl[, 1]) != dim(V_tbl)[1]) {
        warning("Re-ordering samples so that the rows of X match the rows of V. The results will be in terms of the rows of V.")
        combined_tbl <- dplyr::left_join(V_tbl, X_tbl, by = names(V_tbl)[1])
        tmp <- combined_tbl %>%
          dplyr::select(-dplyr::ends_with(".x")) %>%
          dplyr::rename_at(.vars = dplyr::vars(dplyr::ends_with(".y")),
                    .funs = list(~sub("[.]y$", "", .)))
        X_tbl <- tmp
    }
    ## if the absolute abundance-observed columns are scrambled between W and V, change W to match V
    if (sum(names(W_tbl)[-1][1:q_obs] == names(V_tbl)[-1]) != q_obs) {
      warning("Re-ordering columns so that the first q_obs columns of W are in the same order as V.")
      tmp <- W_tbl %>%
        dplyr::select(names(V_tbl), names(W_tbl)[(q_obs + 1):q])
      W_tbl <- tmp
    }

    ## make matrix version of W and V and X, if they aren't already; remove first column
    W_mat <- as.matrix(W)[, -1, drop = FALSE]
    V_mat <- as.matrix(V)[, -1, drop = FALSE]
    X_mat <- as.matrix(X)[, -1, drop = FALSE]
    ## ----------------------------------------
    ## set up the data and initial values lists
    ## ----------------------------------------
    if (dim(X_mat)[2] == 0) {
        # don't use covariate data
        data_lst <- list(W = W_mat, V = V_mat, N = N, q = q, q_obs = q_obs, sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma)
    } else {
        data_lst <- list(W = W_mat, V = V_mat, N = N, q = q, q_obs = q_obs, p = dim(X_mat)[2], X = X_mat, sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma)
    }

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
    ## if q == q_obs, the naive estimator is V
    if (q == q_obs) {
        naive_est <- as.matrix(V_mat)
    } else { ## otherwise, we have to do something
        naive_est <- cbind(as.matrix(V_mat), apply(matrix((q_obs + 1):q), 1, naive_estimator, W_mat, V_mat, 1:q_obs))
    }
    colnames(naive_est) <- colnames(W_mat)
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
            inits_lst <- list(list(log_mu_tilde = log_naive_tilde), rep(list(init = "random"), n_chains - 1))
        } else {
            inits_lst <- list(list(log_mu_tilde = log_naive_tilde, beta_0 = naive_beta, log_Sigma = log(naive_Sigma)))
        }
    }
    ## ----------------------
    ## run the Stan algorithm
    ## ----------------------
    if (dim(X_mat)[2] == 0) {
        mod <- rstan::sampling(stanmodels$no_efficiency, data = data_lst, pars = c("mu", "beta_0", "Sigma"),
                               chains = n_chains, iter = n_iter, warmup = n_burnin, seed = stan_seed,
                               init = inits_lst, ...)
    } else {
        mod <- rstan::sampling(stanmodels$no_efficiency_covariates, data = data_lst, pars = c("mu", "beta_0", "beta_1", "Sigma"), chains = n_chains, iter = n_iter, warmup = n_burnin, seed = stan_seed, init = inits_lst, ...)
    }
    return(mod)
}
