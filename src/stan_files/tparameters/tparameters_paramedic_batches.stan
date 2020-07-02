    //
    vector[q] p[N,K];
    vector[q_obs] log_mu_v[N];
    vector[q] log_mu[N];

    for (i in 1:N) {
        if (d > 0) {
            log_mu[i] = beta_0 + (X[i] * beta_1)' + exp(log_Sigma) .* log_mu_tilde[i];
        }
        else {
            log_mu[i] = beta_0 + exp(log_Sigma) .* log_mu_tilde[i];
        }
        if (alpha_sigma > 0 && kappa_sigma > 0) {
            for (k in 1:K) {
                p[i,k] = softmax(log_mu[i] + log_e[k]);
            }
        }
        else {
            for (k in 1:K) {
                p[i,k] = softmax(log_mu[i]);
            }
        }
        log_mu_v[i] = head(log_mu[i], q_obs);
    }
