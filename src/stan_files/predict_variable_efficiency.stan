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
    vector[q] mu[N_samples, N];
    vector[q] p[N_samples, N];
    vector[q] e[N_samples];
    real V[N_samples, N, q];
    int W[N_samples, N, q];

    for (l in 1:N_samples) {
        // predicted values for e
        if (alpha_sigma > 0 && kappa_sigma > 0) {
            for (j in 1:q) {
                e[l, j] = exp(normal_rng(0, sqrt(sigma_e[l])));
            }
        }
        else {
            for (j in 1:q) {
                e[l, j] = 1;
            }
        }
        for (i in 1:N) {
            // predicted values for mu
            if (d > 0) {
                mu[l, i] = to_vector(exp(normal_rng(beta_0[l] + (X[i] * beta_1[l])', Sigma[l])));
            } else {
                mu[l, i] = to_vector(exp(normal_rng(beta_0[l], Sigma[l])));
            }
            // predicted values for V
            if (alpha_phi > 0 && beta_phi > 0) {
                V[l, i] = neg_binomial_rng(mu[l, i], phi[l, i]);
            }
            else {
                // if mu is *super* large, use normal approximation
                for (k in 1:q) {
                    if (mu[l, i, k] > 1e09) {
                        V[l, i, k] = normal_rng(mu[l, i, k], mu[l, i, k]);
                    } else {
                        V[l, i, k] = poisson_rng(mu[l, i, k]);
                    }
                }
            }
            // predicted values for W
            p[l, i] = softmax(log(mu[l, i]) + log(e[l]));
            W[l, i] = multinomial_rng(p[l, i], M[i]);
        }
    }
}
