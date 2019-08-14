#' Run a Stan algorithm to estimate concentrations from combined absolute and relative abundance data
#'
#' Estimate concentrations (and efficiencies, if applicable) by combining absolute and relative abundance data.
#'
#' @param W The relative abundance data, e.g., from broad range 16S sequencing with "universal" primers. Expects data (e.g., matrix, data.frame, tibble) with sample identifiers in the first column. Sample identifiers must be the same between W and V, and the column must have the same name in W and V.
#' @param V The absolute abundance data, e.g., from taxon-specific absolute primers. Expects data (e.g., matrix, data.frame, tibble) with sample identifiers in the first column. Sample identifiers must be the same between W and V, and the column must have the same name in W and V.
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
#' n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747)
#'
#' @seealso \code{\link[rstan]{stan}} and \code{\link[rstan]{sampling}} for specific usage of the \code{stan} and \code{sampling} functions.
#'
#' @export
run_paramedic <- function(W, V,
                      n_iter = 10500, n_burnin = 10000, n_chains = 4, stan_seed = 4747,
                      inits_lst = NULL,
                      ...) {
    N <- dim(W)[1]
    q <- dim(W)[2] - 1
    q_obs <- dim(V)[2] - 1
    ## --------------
    ## error messages
    ## --------------
    ## error if there aren't the same number of rows in W and V
    if (dim(W)[1] != dim(V)[1]) stop("The number of rows in W and V must match.")
    ## error if there is any difference between the first column of W and V
    if (length(setdiff(as.numeric(unlist(W[, 1])), as.numeric(unlist(V[, 1])))) != 0) stop("W and V must have the same samples.")
    ## error if the id columns are named different things
    if ("data.frame" %in% class(W) | "tbl" %in% class(W)) {
      if (names(W)[1] != names(V)[1]) stop("W and V must have the same name ")  
    } else {
      if (colnames(W)[1] != colnames(V)[1]) stop("W and V must have the same name ")  
    }
    
    ## ---------------------------
    ## pre-processing and warnings
    ## ---------------------------
    ## make tibbles out of W and V, if they aren't already
    W_tbl <- tibble::as_tibble(W)
    V_tbl <- tibble::as_tibble(V)
    ## if the rows are scrambled between W and V, change W to match V
    if (sum(W_tbl[, 1] == V_tbl[, 1]) != dim(V_tbl)[1]) {
      warning("Re-ordering samples so that the rows of W match the rows of V. The results will be in terms of the rows of V.")
      combined_tbl <- dplyr::left_join(V_tbl, W_tbl, by = names(V_tbl)[1])
      tmp <- combined_tbl %>% 
        select(-ends_with(".x")) %>% 
        rename_at(.vars = vars(ends_with(".y")),
                  .funs = funs(sub("[.]y$", "", .)))
      W_tbl <- tmp
    }
    ## if the absolute abundance-observed columns are scrambled between W and V, change W to match V
    if (sum(names(W_tbl)[-1][1:q_obs] == names(V_tbl)[-1]) != q_obs) {
      warning("Re-ordering columns so that the first q_obs columns of W are in the same order as V.")
      tmp <- W_tbl %>% 
        select(names(V_tbl), names(W_tbl)[(q_obs + 1):q])
      W_tbl <- tmp
    }
    
    ## make matrix version of W and V, if they aren't already; remove first column
    W_mat <- as.matrix(W)[, -1]
    V_mat <- as.matrix(V)[, -1]
    ## ----------------------------------------
    ## set up the data and initial values lists
    ## ----------------------------------------
    data_lst <- list(W = W_mat, V = V_mat, N = N, q = q, q_obs = q_obs)
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
    naive_est <- cbind(as.matrix(V_mat), apply(matrix((q_obs + 1):q), 1, naive_estimator, W_mat, V_mat, 1:q_obs))
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
            inits_lst <- list(list(log_mu_tilde = log_naive_tilde, beta = naive_beta, Sigma = naive_Sigma))
        }
    }
    ## ----------------------
    ## run the Stan algorithm
    ## ----------------------
    mod <- rstan::sampling(stanmodels$variable_efficiency, data = data_lst, pars = c("mu", "e", "beta", "Sigma"),
                           chains = n_chains, iter = n_iter, warmup = n_burnin, seed = stan_seed,
                           init = inits_lst, ...)
    return(mod)
}
