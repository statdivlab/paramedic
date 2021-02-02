#' Run a Stan algorithm to estimate concentrations from combined absolute and relative abundance data
#'
#' Estimate concentrations (and efficiencies, if applicable) by combining absolute and relative abundance data.
#'
#' @param W The relative abundance data, e.g., from broad range 16S sequencing with "universal" primers. Expects data (e.g., matrix, data.frame, tibble) with sample identifiers in the first column. Sample identifiers must be the same between W and V, and the column must have the same name in W and V.
#' @param V The absolute abundance data, e.g., from taxon-specific absolute primers. Expects data (e.g., matrix, data.frame, tibble) with sample identifiers in the first column. Sample identifiers must be the same between W and V, and the column must have the same name in W and V.
#' @param X The covariate data. Expects data (e.g., matrix, data.frame, tibble) with sample identifiers in the first column. Sample identifiers must be the same between W, V, and X, and the column must have the same name in W, V, and X. If X only consists of the subject identifiers, then no covariates are used.
#' @param k the number of batches that the relative abundance data W were analyzed in. If k = 0 (the default), then batch effects are not considered.
#' @param n_iter The total number of iterations per chain to be run by the Stan algorithm. Defaults to 10500.
#' @param n_burnin The total number of warmups per chain to be run by the Stan algorithm. Defaults to 10000.
#' @param n_chains The total number of chains to be run by the Stan algorithm. Defaults to 4.
#' @param stan_seed The random number seed to initialize.
#' @param centered If \code{TRUE}, uses a centered parameterization on the true concentrations \eqn{\mu}; if \code{FALSE} (the default), uses a non-centered parameterization.
#' @param inits_lst An optional list of initial values of the parameters. Must be a named list; see \code{\link[rstan]{stan}}.
#' @param sigma_beta Hyperparameter specifying the prior variance on \eqn{\beta_0}. Defaults to \eqn{\sqrt{50}}.
#' @param sigma_Sigma Hyperparameter specifying the prior variance on \eqn{\Sigma}. Defaults to \eqn{\sqrt{50}}.
#' @param alpha_sigma Hyperparameter specifying the shape parameter of the prior distribution on \eqn{\sigma_e}. Defaults to 2.
#' @param kappa_sigma Hyperparameter specifying the scale parameter of the prior distribution on \eqn{\sigma_e}. Defaults to 1.
#' @param alpha_phi Hyperparameter specifying the shape parameter of the prior distribution on \eqn{\phi}. Defaults to 0; a negative binomial model can be specified if both \code{alpha_phi} and \code{beta_phi} are nonzero.
#' @param beta_phi Hyperparameter specifying the rate parameter of the prior distribution on \eqn{\phi}. Defaults to 0; a negative binomial model can be specified if both \code{alpha_phi} and \code{beta_phi} are nonzero.
#' @param sigma_xi Hyperparameters specifying the variance of efficiencies over batches. Only used if \code{k} is greater than zero. Defaults to 1.
#' @param ... other arguments to pass to \code{\link[rstan]{sampling}} (e.g., control).
#'
#' @return An object of class \code{stanfit}.
#'
#' @details We fit a hierarchical model in Stan to the data, with goal to estimate true concentration for all taxa. There are two available hierarchical models. The two only differ in the way that the true concentration is parameterized. While these hierarchical models are mathematically identical, they are fit differently by Stan. The first hierarchical model is a "centered" model, where \deqn{\beta_0 ~ N(0, \sqrt{50})} \deqn{\log \Sigma ~ N(0, \sqrt{50})} \deqn{\log \mu ~ N(\beta_0 + X\beta_1, \Sigma).} We call this model \code{variable_efficiency_centered.stan}. The second hierarchical model is a "noncentered" model, where \deqn{\beta_0 ~ N(0, \sqrt{50})} \deqn{\log \Sigma ~ N(0, \sqrt{50})} \deqn{\log \gamma ~ N(0, I)} \deqn{\log \mu = \sqrt{\Sigma} \log \gamma + \beta_0 + X\beta_1.} We call this model \code{variable_efficiency.stan}.

#' In most cases, we suggest using the noncentered model. However, if the model does not converge in a reasonable amount of time using the noncentered model, consider trying the centered model.
#'
#' @examples
#' # load the package, read in example data
#' library("paramedic")
#' data(example_16S_data)
#' data(example_qPCR_data)
#'
#' # run paramedic (with an extremely small number of iterations, for illustration only)
#' # on only the first 10 taxa
#' mod <- run_paramedic(W = example_16S_data[, 1:10], V = example_qPCR_data,
#' n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747)
#'
#' @seealso \code{\link[rstan]{stan}} and \code{\link[rstan]{sampling}} for specific usage of the \code{stan} and \code{sampling} functions.
#'
#' @export
run_paramedic <- function(W, V, X = V[, 1, drop = FALSE], k = 0,
                      n_iter = 10500, n_burnin = 10000, n_chains = 4, stan_seed = 4747,
                      centered = FALSE, inits_lst = NULL,
                      sigma_beta = sqrt(50), sigma_Sigma = sqrt(50), alpha_sigma = 2, kappa_sigma = 1, 
                      alpha_phi = 0, beta_phi = 0, sigma_xi = 1,
                      ...) {
    # --------------
    # error messages
    # --------------
    check_entered_data(W, V, X, k, inits_lst, sigma_beta, sigma_Sigma, alpha_sigma, kappa_sigma)
    # ---------------------------
    # pre-processing and warnings
    # ---------------------------
    pre_processed_lst <- make_paramedic_tibbles(W, V, X, k, inits_lst, sigma_beta, sigma_Sigma, alpha_sigma, kappa_sigma)
    W_mat <- pre_processed_lst$w_mat
    V_mat <- pre_processed_lst$v_mat
    X_mat <- pre_processed_lst$x_mat
    # ----------------------------------------
    # set up the data and initial values lists
    # ----------------------------------------
    data_inits_lst <- make_paramedic_stan_data(k, W_mat, V_mat, X_mat, inits_lst, 
                                               sigma_beta, sigma_Sigma, 
                                               alpha_sigma, kappa_sigma, 
                                               alpha_phi, beta_phi,
                                               n_chains, centered)
    data_lst <- data_inits_lst$data_lst
    inits_lst <- data_inits_lst$inits_lst
    # ----------------------
    # run the Stan algorithm
    # ----------------------
    pars <- c("mu",
              if (alpha_sigma > 0 & kappa_sigma > 0) "e",
              "beta_0",
              if (ncol(X_mat) > 0) "beta_1",
              "Sigma",
              if (alpha_phi > 0 & beta_phi > 0) "phi",
              if (alpha_sigma > 0 & kappa_sigma > 0) "sigma_e")
    stan_model <- switch((k > 0) + 1, 
                         switch(centered + 1,
                         stanmodels$variable_efficiency,
                         stanmodels$variable_efficiency_centered),
                         stanmodels$variable_efficiency_batches)
    mod <- rstan::sampling(stan_model, data = data_lst, pars = pars,
                           chains = n_chains, iter = n_iter, warmup = n_burnin, seed = stan_seed,
                           init = inits_lst, ...)
    output <- list(
        stan_fit = mod,
        stan_data = data_lst,
        stan_inits = inits_lst,
        stan_model = stan_model,
        summary = rstan::summary(mod, probs = c(0.025, 0.975))$summary
    )
    structure(output, class = "paramedic")
}
