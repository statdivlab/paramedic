data{
    // declares N, q_obs, q, d, K, V, W, X
    // also declares hyperparameters:
    // sigma_beta, sigma_Sigma,
    // alpha_sigma, kappa_sigma,
    // alpha_phi, beta_phi
#include /data/data_paramedic_batches.stan
}
parameters{
    // first-level parameters
    vector[q] log_mu_tilde[N];
    // declares shared parameters log_e,
    // beta_0, beta_1, log_Sigma,
    // sigma_e, phi
#include /parameters/parameters_paramedic_batches.stan
}
transformed parameters{
    // declares p, log_mu_v, log_mu
#include /tparameters/tparameters_paramedic_batches.stan
}
model {
    // specifies priors on beta_0, beta_1, Sigma,
    // sigma_e, log_e, and phi
#include /model/shared_model_paramedic_batches.stan

    for (k in 1:K) {
        for (i in 1:N) {
            log_mu_tilde[i] ~ std_normal();
            if (alpha_phi > 0 && beta_phi > 0) {
                V[k][i] ~ neg_binomial_2_log(log_mu_v[i], phi[i]);
            }
            else {
                V[k][i] ~ poisson_log(log_mu_v[i]);
            }
            W[k][i] ~ multinomial(p[i,k]);
        }
    }
}
generated quantities{
    vector[q] mu[N];
    vector[q] e[K];
    vector[q] Sigma;

    mu = exp(log_mu);
    if (alpha_sigma > 0 && kappa_sigma > 0) {
        e = exp(log_e);
    }
    else {
        for (k in 1:K) {
            e[k] = rep_vector(1, q);
        }
    }
    Sigma = exp(log_Sigma);
}
