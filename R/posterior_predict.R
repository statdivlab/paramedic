#' Draw from the posterior predictive distribution
#' 
#' Predict concentrations (and efficiencies, if applicable), 
#' absolute abundances, and relative abundances based on the posterior 
#' distributions from a previously fitted model resulting from a call to 
#' \code{run_paramedic}.
#' 
#' @param object An object of class \code{"paramedic"}, 
#'   resulting from a call to \code{run_paramedic}.
#' @param W The new relative abundance data, e.g., 
#'   from broad range 16S sequencing with "universal" primers. 
#'   Expects data (e.g., matrix, data.frame, tibble) with sample identifiers 
#'   in the first column. Sample identifiers must be the same between W and V, 
#'   and the column must have the same name in W and V.
#' @param V The new absolute abundance data, e.g., from taxon-specific 
#'   absolute primers. Expects data (e.g., matrix, data.frame, tibble) with 
#'   sample identifiers in the first column. Sample identifiers must be the 
#'   same between W and V, and the column must have the same name in W and V.
#' @param X The new covariate data. Expects data (e.g., matrix, data.frame, 
#'   tibble) with sample identifiers in the first column. Sample identifiers
#'    must be the same between W, V, and X, and the column must have the same 
#'    name in W, V, and X. If X only consists of the subject identifiers, 
#'    then no covariates are used.
#' @param draws the number of draws to return. The default and maximum number 
#'   of draws is the size of the posterior sample.
#' @param alpha_sigma Hyperparameter specifying the shape parameter of the 
#'   prior distribution on \eqn{\sigma_e}. Defaults to 2.
#' @param kappa_sigma Hyperparameter specifying the scale parameter of the 
#'   prior distribution on \eqn{\sigma_e}. Defaults to 1.
#' @param alpha_phi Hyperparameter specifying the shape parameter of the 
#'   prior distribution on \eqn{\phi}. Defaults to 0; a negative binomial 
#'   model can be specified if both \code{alpha_phi} and \code{beta_phi} are 
#'   nonzero.
#' @param beta_phi Hyperparameter specifying the rate parameter of the prior 
#'   distribution on \eqn{\phi}. Defaults to 0; a negative binomial model can 
#'   be specified if both \code{alpha_phi} and \code{beta_phi} are nonzero.
#' @param k the number of batches that the relative abundance data W were 
#'   analyzed in. If k = 0 (the default), then batch effects are not considered. 
#'   (currently not used)
#' @param sigma_xi Hyperparameters specifying the variance of efficiencies 
#'   over batches. Only used if \code{k} is greater than zero. 
#'   Defaults to 1. (currently not used)
#' @param ... Ignored
#'   
#' @return A list of \code{draws} by \code{ncol(W)} matrices (for taxon-level 
#'   parameters) or \code{draws} by \code{nrow(W)} by \code{ncol(W)} arrays 
#'   (for individual-level parameters and data) of simulations from the 
#'   posterior predictive distribution. Each row of the matrices, and each 
#'   of the first dimension of the arrays, denotes the predictions generated 
#'   using a single draw of the model parameters from the posterior distribution.
#'  
#' @aliases posterior_predict
#' @export  
posterior_predict.paramedic <- function(object, W = NULL, V = NULL, 
                                        X = V[, 1, drop = FALSE], 
                                        draws = NULL, 
                                        alpha_sigma = 2, kappa_sigma = 1, 
                                        alpha_phi = 0, beta_phi = 0, 
                                        k = 0, sigma_xi = 1,
                                        ...) {
    # --------------
    # error messages
    # --------------
    check_entered_data(W, V, X, k, alpha_sigma = alpha_sigma, 
                       kappa_sigma = kappa_sigma)
    # ---------------------------
    # pre-processing and warnings
    # ---------------------------
    pre_processed_lst <- make_paramedic_tibbles(W, V, X, k, 
                                                alpha_sigma = alpha_sigma, 
                                                kappa_sigma = kappa_sigma)
    W_mat <- pre_processed_lst$w_mat
    V_mat <- pre_processed_lst$v_mat
    X_mat <- pre_processed_lst$x_mat
    N <- nrow(W_mat)
    q <- ifelse(k > 0, dim(W_mat)[3], dim(W_mat)[2])
    q_obs <- ifelse(k > 0, dim(V_mat)[3], dim(V_mat)[2])
    d <- dim(X_mat)[2]
    M <- rowSums(W_mat)
    # ---------------------------
    # extract relevant quantities 
    # ---------------------------
    posterior_params <- get_pp_params_paramedic(object$stan_fit, draws = draws)
    
    draws <- get_pp_draws_paramedic(posterior_params, N = N, q = q, 
                                    q_obs = q_obs, d = d,
                                    X = X_mat, M = M,
                                    alpha_sigma = alpha_sigma, 
                                    kappa_sigma = kappa_sigma, 
                                    alpha_phi = alpha_phi, 
                                    beta_phi = beta_phi)
    structure(draws, class = c("ppd", class(draws)))
}

# extract relevant quantities
# @param object A stanfit object
# @param draws Number of draws
# @return the matrix of posterior draws for each parameter
get_pp_params_paramedic <- function(object, draws = NULL) {
    ext_obj <- rstan::extract(object)
    S <- nrow(ext_obj$beta_0)
    if (is.null(draws)) {
        draws <- S
    }
    if (draws > S) {
        stop(paste0("'draws' should be <= posterior sample size (", S, ")."))
    }
    samp <- switch((draws < S) + 1, 1:S, sample(1:S, draws))
    
    beta_0 <- ext_obj$beta_0[samp, , drop = FALSE]
    Sigma <- ext_obj$Sigma[samp, , drop = FALSE]
    if (any(grepl("sigma_e", names(ext_obj)))) {
        sigma_e <- ext_obj$sigma_e[samp]
    } else {
        sigma_e <- NULL
    }
    if (any(grepl("beta_1", names(ext_obj)))) {
        beta_1 <- ext_obj$beta_1[samp, , drop = FALSE]
    } else {
        beta_1 <- NULL
    }
    if (any(grepl("phi", names(ext_obj)))) {
        phi <- ext_obj$phi[samp, , drop = FALSE]
    } else {
        phi <- NULL
    }
    list(beta_0 = beta_0, Sigma = Sigma, sigma_e = sigma_e, beta_1 = beta_1,
         phi = phi)
}

# Draw from posterior predictive distribution
# @param params A list of posterior parameters
# @param N The sample size
# @param q The number of taxa with observed relative abundance
# @param q_obs The number of taxa with observed absolute abundance
# @param d The number of covariates
# @param X The covariates (observed)
# @param M The read numbers (observed)
# @param alpha_sigma Hyperparameter specifying the shape parameter of the 
#   prior distribution on \eqn{\sigma_e}. Defaults to 2.
# @param kappa_sigma Hyperparameter specifying the scale parameter of the 
#   prior distribution on \eqn{\sigma_e}. Defaults to 1.
# @param alpha_phi Hyperparameter specifying the shape parameter of the 
#   prior distribution on \eqn{\phi}. Defaults to 0; a negative binomial 
#   model can be specified if both \code{alpha_phi} and \code{beta_phi} are 
#   nonzero.
# @param beta_phi Hyperparameter specifying the rate parameter of the prior 
#   distribution on \eqn{\phi}. Defaults to 0; a negative binomial model can 
#   be specified if both \code{alpha_phi} and \code{beta_phi} are nonzero.
get_pp_draws_paramedic <- function(params, N = NULL, q = NULL, q_obs = NULL, 
                                   d = 0,
                                   X = NULL, M = NULL,
                                   alpha_sigma = 2, kappa_sigma = 1, 
                                   alpha_phi = 0, beta_phi = 0) {
    # taxon-level
    if (alpha_sigma > 0 && kappa_sigma > 0) {
        e <- exp(.pp_gaussian(matrix(0, nrow = length(params$sigma_e),
                                     ncol = q), 
                              sqrt(params$sigma_e))) 
    } else {
        e <- matrix(1, nrow = nrow(params$beta_0), ncol = q)
    }
    if (d > 0) {
        mu_mean <- params$beta_0 + params$beta_1 %*% X
    } else {
        mu_mean <- params$beta_0
    }
    # indidivual- and taxon-level; note first dimension is the draw,
    # second dimension is observation, third dimension is taxon
    mu <- aperm(replicate(N, exp(.pp_gaussian(mu_mean, sqrt(params$Sigma)))), 
                c(3, 2, 1))
    if (alpha_phi > 0 && beta_phi > 0) {
        V <- aperm(sapply(1:N, function(i) {
            .pp_neg_binomial(t(mu[i, , ]), params$phi)
        }, simplify = "array"), c(3, 2, 1))
    } else {
        V <- aperm(sapply(1:N, function(i) {
            .pp_poisson(t(mu[i, , ]))
        }, simplify = "array"), c(3, 2, 1))
    }
    W <- aperm(sapply(1:N, function(i) {
        .pp_multinomial(t(mu[i, , ]), M[i], e)
    }, simplify = "array"), c(3, 2, 1))
    list(mu = mu, e = e, V = V, W = W)
}

# functions to draw from posterior distributions
.pp_gaussian <- function(mu, sigma) {
    t(sapply(1:nrow(mu), function(s) {
        rnorm(ncol(mu), mu[s, ], sigma[s])
    }))
}
.pp_neg_binomial <- function(mu, phi) {
    t(sapply(1:nrow(mu), function(s) {
        rnbinom(ncol(mu), size = phi[s], mu = mu[s, ])
    }))
}
.pp_poisson <- function(lambda) {
    t(sapply(1:nrow(lambda), function(s) {
        rpois(ncol(lambda), lambda = lambda[s, ])
    }))
}
.pp_multinomial <- function(mu, M, e) {
    t(sapply(1:nrow(mu), function(s) {
        p <- softmax(log(mu[s, ]) + log(e[s, ]))
        rmultinom(1, size = M, prob = p)
    }))
}
softmax <- function(z) {
    exp(z) / sum(exp(z))
}