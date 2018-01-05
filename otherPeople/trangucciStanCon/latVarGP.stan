// Trangucci Latent Variable

data {

  int<lower=1> N;
  vector[N] y;
  real x[N];

}

parameters {

  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector[N] eta;

}

transformed parameters {

vector[N] f;
  {
    matrix[N, N] L;
    matrix[N, N] K;
    K = cov_exp_quad(x, alpha, rho);
    for (n in 1:N)
      K[n, n] = K[n, n] + 1e-12;
    L = cholesky_decompose(K);
    f = L * eta;
  }

}

model {

  rho ~ gamma(2, 20);
  alpha ~ normal(0, 1);
  sigma ~ normal(0, 1);
  eta ~ normal(0, 1);
  y ~ normal(f, sigma);

}

// generated quantities {
// 
// vector[N_pred] f_pred;
// vector[N_pred] y_pred;
// f_pred = gp_pred_rng(x_pred, y, x, alpha, rho, sigma);
// for (n in 1:N_pred)
// y_pred[n] = normal_rng(f_pred[n], sigma);
// 
// }

