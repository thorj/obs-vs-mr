data {
  int<lower=1> N;                 // number of protein-phenotype pairs
  int<lower=1> P;                 // number of phenotypes
  int<lower=1> S;                 // number of somamers
  int<lower=1> G;                 // number of phenotype groups
  int<lower=1> K;                 // number of unique proteins targeted by SOMAmers

  int<lower=1, upper=P> phenotype[N]; // Phenotype at i-th index
  int<lower=1, upper=S> somamer[N];   // SOMAmer at i-th index
  // Maps
  int<lower=1, upper=G> group_map[P];     // Map phenotype at i-th index to its group
  int<lower=1, upper=K> protein_map[S]; // maps SOMAmer s to protein k
  vector[N] f_stat;   // harmonic f stat
  vector[N] log_tc;   // log(1 + total_coding_variants)
  
  matrix[S, S] L;   // Cholesky of correlation matrix between SOMAmers

  int<lower=0, upper=1> y[N];
}

parameters {
  real alpha;                      // global intercept

  // Non-centered parameterizations for the hierarchical effects
  vector[G] gamma_raw;             // raw group effects
  real<lower=0> sigma_gamma;

  vector[P] alpha_raw;             // raw phenotype effects
  real<lower=0> sigma_alpha;
  
  vector[S] z_beta;               // raw somamer effects
  real<lower=0> sigma_eps;
  
  vector[K] zeta_raw;            // raw protein effects
  real<lower=0> sigma_zeta;
  
  real theta_f;   // Harmonic F 
  real theta_tc;  // Total coding variants behind phenotype-protein pair
}

transformed parameters {
  vector[G] gamma_g = sigma_gamma * gamma_raw;
  vector[P] alpha_p = sigma_alpha * alpha_raw;
  vector[K] zeta = sigma_zeta * zeta_raw;
  
  vector[S] beta_mu; // Mean somamer effect from zeta
  vector[S] beta_s;
  
  for (s in 1:S) {
    beta_mu[s] = zeta[protein_map[s]];
  }
  beta_s = beta_mu + sigma_eps * (L * z_beta); // Transformed protein effects
}

model {
  // Priors
  alpha ~ normal(0, 1);

  gamma_raw ~ normal(0, 1);
  sigma_gamma ~ exponential(1);

  alpha_raw ~ normal(0, 1);
  sigma_alpha ~ exponential(1);

  zeta_raw ~ normal(0, 1);
  sigma_zeta ~ exponential(1);
  
  z_beta ~ normal(0, 1);
  sigma_eps ~ exponential(1);
 
  theta_f ~ normal(0, 1);
  theta_tc ~ normal(0, 0.5);

  // Likelihood: vectorized the linear predictor
  {
    vector[N] linpred;
    for (i in 1:N) {
      int p = phenotype[i];
      int g = group_map[p];
      int s = somamer[i];
      linpred[i] = alpha + gamma_g[g] + alpha_p[p] + beta_s[s] + 
      theta_f * f_stat[i] + 
      theta_tc * log_tc[i];
    }
    y ~ bernoulli_logit(linpred);
  }
}

generated quantities {
  vector[N] log_lik;

  for (i in 1:N) {
    int p = phenotype[i];
    int g = group_map[p];
    int s = somamer[i];
    log_lik[i] = bernoulli_logit_lpmf(y[i] |
      alpha + gamma_g[g] +
      alpha_p[p] +
      beta_s[s] +
      theta_f * f_stat[i] +
      theta_tc * log_tc[i]
    );
  }
}
