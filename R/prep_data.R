#' Prepare Data for Paramedic
#' 
#' Functions in \code{paramedic} expect data (W, V, and X) in matrix, data.frame
#' , or tibble form, with the first column being sample identifiers. 
#' Sample identifiers must be the same between W, V, and X and the column must 
#' have the same name in W and V. This function checks the incoming data, 
#' transforms it into the data list that is passed to the Stan functions, and 
#' computes any initial values of model parameters.
#' 
#' @param W The relative abundance data, e.g., from broad range 16S sequencing 
#' with "universal" primers. Expects data (e.g., matrix, data.frame, tibble) 
#' with sample identifiers in the first column. Sample identifiers must be the 
#' same between W and V, and the column must have the same name in W and V.
#' @param V The absolute abundance data, e.g., from taxon-specific absolute 
#' primers. Expects data (e.g., matrix, data.frame, tibble) with sample 
#' identifiers in the first column. Sample identifiers must be the same 
#' between W and V, and the column must have the same name in W and V.
#' @param X The covariate data. Expects data (e.g., matrix, data.frame, 
#' tibble) with sample identifiers in the first column. Sample identifiers 
#' must be the same between W, V, and X, and the column must have the same 
#' name in W, V, and X. If X only consists of the subject identifiers, then
#'  no covariates are used.
#' @param k the number of batches; if > 0, expects a column named "batch" 
#' as part of each tibble
#' @param inits_lst An optional list of initial values of the parameters. 
#' Must be a named list; see \code{\link[rstan]{stan}}.
#' @param sigma_beta Hyperparameter specifying the prior variance on 
#' \eqn{\beta_0}. Defaults to \eqn{\sqrt{50}}.
#' @param sigma_Sigma Hyperparameter specifying the prior variance on 
#' \eqn{\Sigma}. Defaults to \eqn{\sqrt{50}}.
#' @param alpha_sigma Hyperparameter specifying the shape parameter of the 
#' prior distribution on \eqn{\sigma_e}. Defaults to 2.
#' @param kappa_sigma Hyperparameter specifying the scale parameter of the
#'  prior distribution on \eqn{\sigma_e}. Defaults to 1.
check_entered_data <- function(W, V, X, k, inits_lst, sigma_beta, sigma_Sigma, 
                               alpha_sigma, kappa_sigma) {
    # error if there aren't the same number of rows in W and V
    if (dim(W)[1] != dim(V)[1]) stop("The number of rows in W and V must match.")
    # error if there aren't the same number of rows in V and X
    if (dim(V)[1] != dim(X)[1]) stop("The number of rows in V and X must match.")
    # error if there is any difference between the first column of W and V
    if (length(setdiff(as.numeric(unlist(W[, 1])), as.numeric(unlist(V[, 1])))) != 0) stop("W and V must have the same samples.")
    # error if there is any difference between the first column of X and V
    if (length(setdiff(as.numeric(unlist(V[, 1])), as.numeric(unlist(X[, 1])))) != 0) stop("V and X must have the same samples.")
    # error if the id columns are named different things
    if ("data.frame" %in% class(W) | "tbl" %in% class(W)) {
        if (names(W)[1] != names(V)[1]) stop("W and V must have the same name ")
        if (names(W)[1] != names(X)[1]) stop("W and X must have the same name.")
    } else {
        if (colnames(W)[1] != colnames(V)[1]) stop("W and V must have the same name ")
        if (colnames(W)[1] != colnames(X)[1]) stop("W and X must have the same name.")
    }
    q <- dim(W)[2] - 1 - ifelse(k > 0, 1, 0)
    q_obs <- dim(V)[2] - 1 - ifelse(k > 0, 1, 0)
    # error if q < q_obs
    if (q < q_obs) stop("V must have fewer taxa than W (or the same number of taxa).")
    return(invisible(NULL))
}

#' @describeIn check_entered_data Make tibbles from entered data
#' @importFrom tibble as_tibble 
#' @importFrom dplyr select rename_at ends_with .data left_join
#' @importFrom rlang !! sym
make_paramedic_tibbles <- function(W, V, X, k, inits_lst, 
                                   sigma_beta, sigma_Sigma, 
                                   alpha_sigma, kappa_sigma) {
    # make tibbles out of W and V, if they aren't already
    W_tbl <- tibble::as_tibble(W)
    V_tbl <- tibble::as_tibble(V)
    X_tbl <- tibble::as_tibble(X)
    q <- dim(W_tbl)[2] - 1 - ifelse(k > 0, 1, 0)
    q_obs <- dim(V_tbl)[2] - 1 - ifelse(k > 0, 1, 0)
    # if the rows are scrambled between W and V, change W to match V
    if (sum(W_tbl[, 1] == V_tbl[, 1]) != dim(V_tbl)[1]) {
        warning("Re-ordering samples so that the rows of W match the rows of V. The results will be in terms of the rows of V.")
        combined_tbl <- dplyr::left_join(V_tbl, W_tbl, by = names(V_tbl)[1])
        tmp <- combined_tbl %>%
            dplyr::select(-dplyr::ends_with(".x")) %>%
            dplyr::rename_at(.vars = dplyr::vars(dplyr::ends_with(".y")),
                             .funs = list(~sub("[.]y$", "", .)))
        W_tbl <- tmp
    }
    # if the rows are scrambled between X and V, change X to match V
    if (sum(X_tbl[, 1] == V_tbl[, 1]) != dim(V_tbl)[1]) {
        warning("Re-ordering samples so that the rows of X match the rows of V. The results will be in terms of the rows of V.")
        combined_tbl <- dplyr::left_join(V_tbl, X_tbl, by = names(V_tbl)[1])
        tmp <- combined_tbl %>%
            dplyr::select(-dplyr::ends_with(".x")) %>%
            dplyr::rename_at(.vars = dplyr::vars(dplyr::ends_with(".y")),
                             .funs = list(~sub("[.]y$", "", .)))
        X_tbl <- tmp
    }
    # if the absolute abundance-observed columns are scrambled between W and V,
    # change W to match V
    if (sum(names(W_tbl)[-1][1:q_obs] == names(V_tbl)[-1]) != q_obs) {
        warning("Re-ordering columns so that the first q_obs columns of W are in the same order as V.")
        tmp <- W_tbl %>%
            dplyr::select(names(V_tbl), names(W_tbl)[(q_obs + 1):q])
        W_tbl <- tmp
    }
    
    # make matrix version of W and V, if they aren't already; remove first column
    # if k > 0, make an array with K in the first position
    W_mat <- array(NA, dim = c(k, nrow(W), q))
    V_mat <- array(NA, dim = c(k, nrow(V), q_obs))
    X_mat <- array(NA, dim = c(nrow(X), ncol(X) - 1))
    if (k > 0) {
        for (j in 1:k) {
            W_mat[j, , ] <- (W %>% 
                                 dplyr::filter(!! rlang::sym("batch") == j) %>% 
                                 dplyr::select(- !! rlang::sym("batch")) %>% 
                                 as.matrix())[, -1, drop = FALSE]
            V_mat[j, , ] <- (V %>% 
                                 dplyr::filter(!! rlang::sym("batch") == j) %>%
                                 dplyr::select(- !! rlang::sym("batch")) %>% 
                                 as.matrix())[, -1, drop = FALSE]
        }
        X_mat <- as.matrix(X)[, -1, drop = FALSE]
    } else {
        W_mat <- as.matrix(W)[, -1, drop = FALSE]
        V_mat <- as.matrix(V)[, -1, drop = FALSE]
        X_mat <- as.matrix(X)[, -1, drop = FALSE]    
    }
    
    return(list(w_tbl = W_tbl, v_tbl = V_tbl, x_tbl = X_tbl, 
                w_mat = W_mat, v_mat = V_mat, x_mat = X_mat))
}

#' @param k the number of batches (0 if a single batch)
#' @param W_mat the pre-processed matrix W
#' @param V_mat the pre-processed matrix V
#' @param X_mat the pre-processed matrix X
#' @param inits_lst An optional list of initial values of the parameters. 
#' Must be a named list; see \code{\link[rstan]{stan}}.
#' @param sigma_beta Hyperparameter specifying the prior variance on 
#' \eqn{\beta_0}. Defaults to \eqn{\sqrt{50}}.
#' @param sigma_Sigma Hyperparameter specifying the prior variance on 
#' \eqn{\Sigma}. Defaults to \eqn{\sqrt{50}}.
#' @param alpha_sigma Hyperparameter specifying the shape parameter of the 
#' prior distribution on \eqn{\sigma_e}. Defaults to 2.
#' @param kappa_sigma Hyperparameter specifying the scale parameter of the
#'  prior distribution on \eqn{\sigma_e}. Defaults to 1.
#' @param alpha_phi Hyperparameter specifying the shape parameter of the 
#' prior distribution on \eqn{phi} (the optional Negative Binomial 
#' dispersion parameter). Defaults to 0 (in which case a Poisson distribution
#' on V is used).
#' @param beta_phi Hyperparameter specifying the scale parameter of the 
#' prior distribution on \eqn{phi} (the optional Negative Binomial 
#' dispersion parameter). Defaults to 0 (in which case a Poisson distribution
#' on V is used).
#' @param n_chains the number of chains to run
#' @param centered whether or not to use the centered parameterization of the 
#' stan algorithm. Defaults to \code{FALSE}.
#' @describeIn check_entered_data Make data and initial values lists to
#'  pass to stan
#' @importFrom stats var
make_paramedic_stan_data <- function(k, W_mat, V_mat, X_mat, inits_lst, 
                                     sigma_beta, sigma_Sigma, 
                                     alpha_sigma, kappa_sigma, 
                                     alpha_phi, beta_phi,
                                     n_chains, centered = FALSE) {
    N <- ifelse(k > 0, dim(W_mat)[2], dim(W_mat)[1])
    q <- ifelse(k > 0, dim(W_mat)[3], dim(W_mat)[2])
    q_obs <- ifelse(k > 0, dim(V_mat)[3], dim(V_mat)[2])
    d <- dim(X_mat)[2]
    data_lst <- list(W = W_mat, V = V_mat, X = X_mat,
                     N = N, q = q, q_obs = q_obs, d = d,
                     K = k,
                     sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                     alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                     alpha_phi = alpha_phi, beta_phi = beta_phi)
    # get inits from the naive estimator
    naive_estimator <- function(idx, relatives, absolutes, known_absolute) {
        # get the sum of the relative abundances
        sum_relative <- rowSums(relatives[, known_absolute, drop = FALSE])
        # get the sum of the absolutes
        sum_absolute <- rowSums(absolutes)
        
        # naive estimator
        absolute <- relatives[, idx]*sum_absolute/sum_relative
        # if sum_relative was zero, predict zero (?)
        absolute <- ifelse (sum_relative == 0, 0, absolute)
        return(absolute)
    }
    # if q == q_obs, the naive estimator is V
    if (q == q_obs) {
        if (k > 0) {
            naive_est <- apply(V_mat, c(2, 3), mean)
        } else {
            naive_est <- as.matrix(V_mat)
        }
    } else { # otherwise, we have to do something
        if (k > 0) {
            naive_v <- apply(V_mat, c(2, 3), mean)
            naive_ests_arr <- array(NA, dim = c(k, dim(W_mat)[2], 
                                                length((q_obs + 1):q)))
            for (l in 1:k) {
                naive_ests_arr[l, , ] <- apply(matrix((q_obs + 1):q), 1, 
                                               naive_estimator, W_mat[l, , ], 
                                               V_mat[l, , ], 1:q_obs)   
            }
            naive_ests <- apply(naive_ests_arr, c(2, 3), mean)
        } else {
            naive_v <- as.matrix(V_mat)
            naive_ests <- apply(matrix((q_obs + 1):q), 1, naive_estimator, 
                                W_mat, V_mat, 1:q_obs)   
        }
        if (N == 1) {
          naive_ests <- matrix(naive_ests, nrow = 1)
        }
        naive_est <- cbind(naive_v, naive_ests)
    }
    colnames(naive_est) <- colnames(W_mat)
    log_naive <- ifelse(is.infinite(log(naive_est)), 0, log(naive_est))
    naive_beta <- colMeans(log_naive, na.rm = TRUE)
    naive_Sigma <- diag(stats::var(log_naive, na.rm = TRUE))
    tmp <- ifelse(naive_Sigma < 1e-3, 1e-3, naive_Sigma)
    naive_Sigma <- tmp
    log_naive_tilde <- sweep(sweep(log_naive, 2, naive_beta, FUN = "-"), 2, 
                             naive_Sigma, FUN = "/")
    tmp <- ifelse(is.na(log_naive_tilde), 0, log_naive_tilde)
    log_naive_tilde <- tmp
    if (is.null(inits_lst)) { # create inits if not passed in
        if (centered) {
            inits_lst <- list(list(mu = ifelse(naive_est == 0, 1e-4, 
                                               naive_est)))
        } else {
            inits_lst <- list(list(log_mu_tilde = log_naive_tilde))
        }
        if (n_chains > 1) {
            inits_lst <- c(inits_lst, replicate(n_chains - 1, 
                                                list(init = "random"), 
                                                simplify = FALSE))
        } else {
            inits_lst <- list(c(inits_lst[[1]], list(log_Sigma = 
                                                         log(naive_Sigma))))
        }
    }
    return(list(data_lst = data_lst, inits_lst = inits_lst))
}