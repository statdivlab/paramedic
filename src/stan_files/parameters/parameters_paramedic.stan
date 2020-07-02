    // shared parameters
    vector[q] log_e[K];
    // second-level hyperparameters
    // control the mus
    vector[q] beta_0;
    matrix[d,q] beta_1;
    vector[q] log_Sigma;
    // control the es
    real<lower=0> sigma_e;
    vector<lower=0>[q] xi;
    // control overdispersion in V
    vector<lower=0>[N] phi;
