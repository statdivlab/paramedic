data{
    int<lower=1> N;
    int<lower=1> q_obs;
    int<lower=1> q;
    int<lower=0> d;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
    matrix[N,d] X;
    // hyperparameters
    real sigma_beta;
    real sigma_Sigma;
    // 0 for both = efficiency-naive model
    // otherwise, fit varying-efficiency model
    real<lower=0> alpha_sigma;
    real<lower=0> kappa_sigma;
    // 0 for both = Poisson model
    // otherwise, fit negative binomial
    real<lower=0> alpha_phi;
    real<lower=0> beta_phi;
}
parameters{
    // first-level parameters
    vector[q] log_mu_tilde[N];
    vector[q] log_e;
    // second-level hyperparameters
    vector[q] beta_0;
    matrix[d,q] beta_1;
    vector[q] log_Sigma;
    real<lower=0> sigma_e;
    vector<lower=0>[N] phi;
}
transformed parameters{
    vector[q] p[N];
    vector[q_obs] log_mu_v[N];
    vector[q] log_mu[N];

    for (i in 1:N) {
        if (d > 0)
            log_mu[i] = beta_0 + (X[i] * beta_1)' + exp(log_Sigma) .* log_mu_tilde[i];
        else
            log_mu[i] = beta_0 + exp(log_Sigma) .* log_mu_tilde[i];

        if (alpha_sigma > 0 && kappa_sigma > 0)
            p[i] = softmax(log_mu[i] + log_e);
        else
            p[i] = softmax(log_mu[i]);

        log_mu_v[i] = head(log_mu[i], q_obs);
    }
}
model {
    // hierarchical model
    beta_0 ~ normal(0, sigma_beta);
    log_Sigma ~ normal(0, sigma_Sigma);

    if (alpha_sigma > 0 && kappa_sigma > 0)
        sigma_e ~ inv_gamma(alpha_sigma, kappa_sigma);
        log_e ~ normal(0, sqrt(sigma_e));

    if (d > 0)
        for (j in 1:q) {
            beta_1[:,j] ~ std_normal();
        }

    if (alpha_phi > 0 && beta_phi > 0)
        phi ~ gamma(alpha_phi, beta_phi);

    for (i in 1:N) {
        log_mu_tilde[i] ~ std_normal();
        if (alpha_phi > 0 && beta_phi > 0)
            V[i] ~ neg_binomial_2_log(log_mu_v[i], phi[i]);
        else
            V[i] ~ poisson_log(log_mu_v[i]);

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
