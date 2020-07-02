data{
    // declares N, q_obs, q, d, V, W, X
    // also declares hyperparameters:
    // sigma_beta, sigma_Sigma,
    // alpha_sigma, kappa_sigma,
    // alpha_phi, beta_phi
#include /data/data_paramedic.stan
}
parameters{
    // first-level parameters
    vector[q] log_mu[N];
    // declares shared parameters log_e,
    // beta_0, beta_1, log_Sigma,
    // sigma_e, phi
#include /parameters/parameters_paramedic.stan
}
transformed parameters{
    // declares p, log_mu_v
#include /tparameters/tparameters_paramedic_centered.stan
}
model {
    // specifies priors on beta_0, beta_1, Sigma,
    // sigma_e, log_e, and phi
#include /model/shared_model_paramedic.stan

    for (i in 1:N){
        if (d > 0) {
            log_mu[i] ~ normal(beta_0 + (X[i] * beta_1)', exp(log_Sigma));
        }
        else {
            log_mu[i] ~ normal(beta_0, exp(log_Sigma));
        }
        if (alpha_phi > 0 && beta_phi > 0) {
            V[i] ~ neg_binomial_2_log(log_mu_v[i], phi[i]);
        }
        else {
            V[i] ~ poisson_log(log_mu_v[i]);
        }
        W[i] ~ multinomial(p[i]);
    }
}
generated quantities{
    vector[q] mu[N];
    vector[q] e;
    vector[q] Sigma;

    mu = exp(log_mu);
    if (alpha_sigma > 0 && kappa_sigma > 0)
        e = exp(log_e);
    else
        e = rep_vector(1, q);

    Sigma = exp(log_Sigma);
}
