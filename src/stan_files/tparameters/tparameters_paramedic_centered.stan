    // parameters for centered model
    vector[q] p[N];
    vector[q_obs] log_mu_v[N];

    for (i in 1:N) {
        if (alpha_sigma > 0 && kappa_sigma > 0) {
            p[i] = softmax(log_mu[i] + log_e);
        }
        else {
            p[i] = softmax(log_mu[i]);
        }
        log_mu_v[i] = head(log_mu[i], q_obs);
    }
