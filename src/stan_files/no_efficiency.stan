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
    vector<lower=0>[q] mu[N];
    for (i in 1:N){
        for (j in 1:q) {
            mu[i,j] = exp(beta_0[j] + Sigma[j] * log_mu_tilde[i,j]);
        }
    }
}
model {
    // hierarchical model
    beta_0 ~ normal(0, sigma_beta);
    Sigma ~ lognormal(0, sigma_Sigma);

    for (i in 1:N){
        log_mu_tilde[i] ~ std_normal();
        V[i] ~ poisson(head(mu[i],q_obs));
        W[i] ~ multinomial((mu[i])/sum(mu[i]));
    }
}
