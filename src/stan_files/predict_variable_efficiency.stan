data{
    // declares N, q_obs, q, d, X
    // also declares posterior means and hyperparameters:
    // beta_0, beta_1, Sigma, sigma_e
    // alpha_sigma, kappa_sigma,
    // alpha_phi, beta_phi
#include /data/data_paramedic_predict.stan
}
// note that parameters, transformed parameters, and model are empty
// because we use the posterior distributions from the fit on training data
parameters{
}
transformed parameters{
}
model {
}
generated quantities{
    vector[q] mu[N];
    vector[q] e;
    vector[q] V[N];
    vector[q] W[N];

    // predicted values for e
    if (alpha_sigma > 0 && kappa_sigma > 0) {
        e = exp(normal_rng(0, sqrt(sigma_e)));
    }
    else {
        e = rep_vector(1, q);
    }
    for (i in 1:N) {
        // predicted values for mu
        if (d > 0) {
            mu[i] = exp(normal_rng(beta_0 + (X[i] * beta_1)', Sigma));
        } else {
            mu[i] = exp(normal_rng(beta_0, Sigma));
        }
        // predicted values for V
        if (alpha_phi > 0 && beta_phi > 0) {
            V[i] = neg_binomial_rng(mu[i], phi[i]);
        }
        else {
            V[i] = poisson_rng(mu[i]);
        }
        // predicted values for W
        if (alpha_sigma > 0 && kappa_sigma > 0) {
            W[i] = multinomial_rng(softmax(log(mu[i]) + log(e)), M[i]);
        }
        else {
            W[i] = multinomial_rng(softmax(log(mu[i])), M[i]);
        }
    }
}
