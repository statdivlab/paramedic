#' Run a Stan algorithm to estimate concentrations from combined absolute and relative abundance data
#'
#' Estimate concentrations (and efficiencies, if applicable) by combining absolute and relative abundance data.
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
#' @param alpha_sigma Hyperparameter specifying the shape parameter of the prior distribution on \eqn{\sigma_e}. Defaults to 2.
#' @param kappa_sigma Hyperparameter specifying the scale parameter of the prior distribution on \eqn{\sigma_e}. Defaults to 1.
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
#' ## run paramedic centered (with an extremely small number of iterations, for illustration only)
#' ## on only the first 10 taxa
#' mod <- run_paramedic_centered(W = example_16S_data[, 1:10], V = example_qPCR_data,
#' n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747)
#'
#' @seealso \code{\link[rstan]{stan}} and \code{\link[rstan]{sampling}} for specific usage of the \code{stan} and \code{sampling} functions.
#' 
#' @importFrom rstan sampling
#' @export
run_paramedic_centered <- function(W, V, X = V[, 1, drop = FALSE],
                          n_iter = 10500, n_burnin = 10000, n_chains = 4, stan_seed = 4747,
                          inits_lst = NULL, sigma_beta = sqrt(50), sigma_Sigma = sqrt(50), alpha_sigma = 2, kappa_sigma = 1,
                          ...) {
    ## --------------
    ## error messages
    ## --------------
    check_entered_data(W, V, X, inits_lst, sigma_beta, sigma_Sigma, alpha_sigma, kappa_sigma)
    ## ---------------------------
    ## pre-processing and warnings
    ## ---------------------------
    pre_processed_lst <- make_paramedic_tibbles(W, V, X, inits_lst, sigma_beta, sigma_Sigma, alpha_sigma, kappa_sigma)
    W_mat <- pre_processed_lst$w_mat
    V_mat <- pre_processed_lst$v_mat
    X_mat <- pre_processed_lst$x_mat
    ## ----------------------------------------
    ## set up the data and initial values lists
    ## ----------------------------------------
    data_inits_lst <- make_paramedic_stan_data(W_mat, V_mat, X_mat, inits_lst, sigma_beta, sigma_Sigma, alpha_sigma, kappa_sigma, n_chains, centered = TRUE)
    data_lst <- data_inits_lst$data_lst
    inits_lst <- data_inits_lst$inits_lst
    ## run the Stan algorithm
    if (dim(X_mat)[2] == 0) {
        mod <- rstan::sampling(stanmodels$variable_efficiency_centered, data = data_lst, pars = c("mu", "e", "beta_0", "Sigma", "sigma_e"),
                               chains = n_chains, iter = n_iter, warmup = n_burnin, seed = stan_seed,
                               init = inits_lst, ...)
    } else {
        mod <- rstan::sampling(stanmodels$variable_efficiency_centered_covariates, data = data_lst, pars = c("mu", "e", "beta_0", "beta_1", "Sigma", "sigma_e"),
                               chains = n_chains, iter = n_iter, warmup = n_burnin, seed = stan_seed,
                               init = inits_lst, ...)
    }
    return(mod)
}
