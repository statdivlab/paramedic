#' Generate predictions based on a paramedic model fit
#'
#' Predict concentrations (and efficiencies, if applicable), absolute abundances, and relative abundances based on the posterior distributions from a previously fitted model resulting from a call to \code{run_paramedic.
#'
#' @param W The new relative abundance data, e.g., from broad range 16S sequencing with "universal" primers. Expects data (e.g., matrix, data.frame, tibble) with sample identifiers in the first column. Sample identifiers must be the same between W and V, and the column must have the same name in W and V.
#' @param V The new absolute abundance data, e.g., from taxon-specific absolute primers. Expects data (e.g., matrix, data.frame, tibble) with sample identifiers in the first column. Sample identifiers must be the same between W and V, and the column must have the same name in W and V.
#' @param X The new covariate data. Expects data (e.g., matrix, data.frame, tibble) with sample identifiers in the first column. Sample identifiers must be the same between W, V, and X, and the column must have the same name in W, V, and X. If X only consists of the subject identifiers, then no covariates are used.
#' @param beta_0 The posterior distribution on \eqn{\beta_0} (the intercept in the distribution on concentrations \eqn{\mu}), resulting from a call to \code{run_paramedic}.
#' @param beta_1 The posterior distribution on \eqn{\beta_1} (the slope in the distribution on concentrations \eqn{\mu}; only used if X consists of more than just sample identifiers), resulting from a call to \code{run_paramedic}.
#' @param Sigma The posterior distribution on \eqn{\Sigma} (the variance in the distribution on concentrations \eqn{\mu}).
#' @param sigma_e The posterior distribution on \eqn{\sigma_e} (the variance of the efficiencies); only used if both \code{alpha_sigma} and \code{kappa_sigma} below are nonzero.
#' @param phi The posterior distribution on \eqn{\phi} (the overdispersion parameter on V); only used if both \code{alpha_phi} and \code{beta_phi} below are nonzero.
#' @param k the number of batches that the relative abundance data W were analyzed in. If k = 0 (the default), then batch effects are not considered. (currently not used)
#' @param n_iter The total number of iterations per chain to be run by the Stan algorithm. Defaults to 10500.
#' @param n_burnin The total number of warmups per chain to be run by the Stan algorithm. Defaults to 10000.
#' @param n_chains The total number of chains to be run by the Stan algorithm. Defaults to 4.
#' @param stan_seed The random number seed to initialize.
#' @param alpha_sigma Hyperparameter specifying the shape parameter of the prior distribution on \eqn{\sigma_e}. Defaults to 2.
#' @param kappa_sigma Hyperparameter specifying the scale parameter of the prior distribution on \eqn{\sigma_e}. Defaults to 1.
#' @param alpha_phi Hyperparameter specifying the shape parameter of the prior distribution on \eqn{\phi}. Defaults to 0; a negative binomial model can be specified if both \code{alpha_phi} and \code{beta_phi} are nonzero.
#' @param beta_phi Hyperparameter specifying the rate parameter of the prior distribution on \eqn{\phi}. Defaults to 0; a negative binomial model can be specified if both \code{alpha_phi} and \code{beta_phi} are nonzero.
#' @param sigma_xi Hyperparameters specifying the variance of efficiencies over batches. Only used if \code{k} is greater than zero. Defaults to 1. (currently not used)
#' @param ... other arguments to pass to \code{\link[rstan]{sampling}} (e.g., control).
#'
#' @return An object of class \code{stanfit}.
#'
#' @details Using the posterior distributions from a call to \code{run_paramedic}, generate predicted values.
#'
#' @examples
#' ## load the package, read in example data
#' library("paramedic")
#' data(example_16S_data)
#' data(example_qPCR_data)
#' ## generate a train/test split
#' set.seed(1234)
#' folds <- sample(1:2, size = nrow(example_16S_data), replace = TRUE)
#'
#' ## run paramedic (with an extremely small number of iterations, for illustration only)
#' ## on only the first 10 taxa, and on a subset of the data
#' mod <- run_paramedic(W = example_16S_data[folds == 1, 1:10], 
#' V = example_qPCR_data[folds == 1, ], n_iter = 30, n_burnin = 25, 
#' n_chains = 1, stan_seed = 4747)
#' 
#' # generate predictions on test data
#' ext_mod <- rstan::extract(mod)
#' beta0_post <- ext_mod$beta_0
#' sigma_post <- ext_mod$Sigma
#' sigmae_post <- ext_mod$sigma_e
#' pred_mod <- predict_paramedic(beta_0 = beta0_post, Sigma = sigma_post,
#' sigma_e = sigmae_post, ) 
#'
#' @seealso \code{\link[rstan]{stan}} and \code{\link[rstan]{sampling}} for specific usage of the \code{stan} and \code{sampling} functions; \code{\link[paramedic]{run_paramedic}} for usage of the \code{run_paramedic} function.
#'
#' @export
predict_paramedic <- function(W, V, X = V[, 1, drop = FALSE], beta_0, 
                          beta_1 = array(0, dim = c(nrow(beta_0), ncol(X) - 1,
                                                    ncol(W) - 1)),
                          Sigma, sigma_e = matrix(0, nrow = nrow(Sigma), 
                                           ncol = ncol(Sigma)),
                          phi = matrix(1, nrow = nrow(Sigma), ncol = nrow(X)),
                          k = 0, n_iter = 10500, n_burnin = 10000, n_chains = 4, 
                          stan_seed = 4747, alpha_sigma = 2, kappa_sigma = 1, 
                          alpha_phi = 0, beta_phi = 0, sigma_xi = 1,
                          ...) {
    ## --------------
    ## error messages
    ## --------------
    check_entered_data(W, V, X, k, alpha_sigma = alpha_sigma, 
                       kappa_sigma = kappa_sigma)
    ## ---------------------------
    ## pre-processing and warnings
    ## ---------------------------
    pre_processed_lst <- make_paramedic_tibbles(W, V, X, k, 
                                                alpha_sigma = alpha_sigma, 
                                                kappa_sigma = kappa_sigma)
    W_mat <- pre_processed_lst$w_mat
    V_mat <- pre_processed_lst$v_mat
    X_mat <- pre_processed_lst$x_mat
    ## ----------------------------------------
    ## set up the data list
    ## ----------------------------------------
    N <- ifelse(k > 0, dim(W_mat)[2], dim(W_mat)[1])
    N_samples <- nrow(Sigma)
    q <- ifelse(k > 0, dim(W_mat)[3], dim(W_mat)[2])
    q_obs <- ifelse(k > 0, dim(V_mat)[3], dim(V_mat)[2])
    d <- dim(X_mat)[2]
    M <- rowSums(W_mat)
    data_lst <- list(N = N, N_samples = N_samples, q_obs = q_obs, q = q,
                     M = M, d = d, X = X_mat, sigma_e = sigma_e, 
                     beta_0 = beta_0, beta_1 = beta_1, Sigma = Sigma,
                     phi = phi, alpha_sigma = alpha_sigma, 
                     kappa_sigma = kappa_sigma, alpha_phi = alpha_phi,
                     beta_phi = beta_phi)
    ## ----------------------
    ## run the Stan algorithm
    ## ----------------------
    pars <- c("mu",
              switch((alpha_sigma > 0 & kappa_sigma > 0) + 1, NULL, "e"),
              "V", "W")
    stan_model <- stanmodels$predict_variable_efficiency
    mod <- rstan::sampling(stan_model, data = data_lst, pars = pars,
                           chains = n_chains, iter = n_iter, warmup = n_burnin, 
                           seed = stan_seed, algorithm = "Fixed_param", ...)
    return(mod)
}
