// Trangucci Marginal Likelihood

data {

  int<lower=1> N;
  vector[N] y;
  real x[N];

}

transformed data {

  vector[N] zeros;
  zeros = rep_vector(0, N);

}

parameters {

  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;

}

model {

matrix[N, N] L_cov;
  {
    matrix[N, N] cov;
    cov = cov_exp_quad(x, alpha, rho);
    for (n in 1:N)
      cov[n, n] = cov[n, n] + square(sigma);
    L_cov = cholesky_decompose(cov);
  }
  // rho ~ gamma(2, 20);
  // alpha ~ normal(0, 1);
  // sigma ~ normal(0, 1);
  y ~ multi_normal_cholesky(zeros, L_cov);

}

