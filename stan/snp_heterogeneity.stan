data {
  int<lower=1> N;                  // Number of SNPs
  vector[N] beta_hat;              // Wald ratios
  matrix[N, N] Sigma;              // Covariance matrix of Wald ratios (LD-adjusted)
}

parameters {
  real theta;                     // Overall causal effect
  real<lower=0> tau;              // Between-SNP heterogeneity (SD)
}

model {
  matrix[N, N] Sigma_total;
  matrix[N, N] L_Sigma_total;

  // Priors
  theta ~ normal(0, 1);
  tau ~ cauchy(0, 0.1);

  // Total covariance = measurement error + heterogeneity
  Sigma_total = Sigma + diag_matrix(rep_vector(tau^2, N));
  L_Sigma_total = cholesky_decompose(Sigma_total);

  // Multivariate normal likelihood
  beta_hat ~ multi_normal_cholesky(rep_vector(theta, N), L_Sigma_total);
}
