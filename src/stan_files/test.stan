data{
    int<lower=1> N;
    int<lower=1> q_obs;
    int<lower=1> q;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
    // hyperparameters
    real sigma_beta;
    real sigma_Sigma;
}
parameters{
    vector[q] log_mu_tilde[N];
    vector[q] beta_0;
    vector<lower=0>[q] Sigma;
}
transformed parameters{
    vector[q] log_mu[N];
    for (i in 1:N){
        for (j in 1:q) {
            log_mu[i,j] = beta_0[j] + Sigma[j] * log_mu_tilde[i,j];
        }
    }
}
model {
    // hierarchical model
    beta_0 ~ normal(0, sigma_beta);
    Sigma ~ lognormal(0, sigma_Sigma);

    for (i in 1:N) {
        log_mu_tilde[i] ~ std_normal();
        V[i] ~ poisson_log(head(log_mu[i],q_obs));
        W[i] ~ multinomial(softmax(log_mu[i]));
    }
}
generated quantities{
    vector[q] mu[N];
    for (i in 1:N){
        for (j in 1:q) {
            mu[i,j] = exp(log_mu[i,j]);
        }
    }
}
