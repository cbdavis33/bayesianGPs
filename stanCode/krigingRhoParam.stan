functions {

  // See pg. 256 in Stan User Guide for version 2.17.0 for a better way to create this matrix.
  matrix getR(real[] x,
              vector rho) {
                
    int N = size(x);
    matrix[N, N] R;
    for (i in 1:(N-1)) {
      R[i, i] = 1;
      for (j in (i + 1):N) {
        R[i, j] = rho^(16 * (x[i] - x[j])^2);
        R[j, i] = R[i, j];
      }
    }
    R[N, N] = 1;
    //return cholesky_decompose(R);
    return R;
  }
  
  // matrix getR(int N,
  //             real[] x_is,
  //             real rho){
  // 
  //   matrix[N, N] R;
  // 
  //   for(j in 1:N){
  //     for(i in 1:N){
  //       if(i >= j){
  //         R[i, j] = rho^(16 * (x_is[i] - x_is[j])^2);
  //       }else{
  //         R[i, j] = R[j, i];
  //       }
  //     }
  //   }
  //   return R;                
  // }
  
  matrix getC(matrix R,
              real sigma,
              real sigmaEps){
    
    int N = rows(R);
    matrix[N, N] C = sigma^2 * R;
    
    for (n in 1:N)
        C[n, n] = C[n, n] + sigmaEps^2;
    // C += diag_matrix(rep_vector(sigmaEps^2,N));
    
    return C;              
  
  }

  vector gp_pred_rng(real[] x_pred,
                     vector y_is,
                     real[] x_is,
                     real sigma,
                     real rho,
                     real sigmaEps) {
                       
    vector[size(x_pred)] f_pred;
    int N_pred = size(x_pred);
    int N = rows(y_is);
    
    {
      matrix[N, N] R = getR(x_is, rho);
      matrix[N, N] C = getC(R, sigma, sigmaEps);
      matrix[N, N] L_C = cholesky_decompose(C);
      matrix[N_pred, N_pred] CStar;
      
      vector[N] K_div_y_is_left = mdivide_left_tri_low(L_C, y_is);
      vector[N] K_div_y_is = mdivide_right_tri_low(K_div_y_is_left',L_C)';
      
      matrix[N, N_pred] k_x_is_x_pred;
      matrix[N, N_pred] v_pred;
      vector[N_pred] f_pred_mu;
      
      matrix[N_pred, N_pred] nug_pred;
    
      
      
      K_div_y_is = mdivide_right_tri_low(K_div_y_is',L_C)';
      k_x_is_x_pred = cov_exp_quad(x_is, x_pred, sigma, rho);
      f_pred_mu = (k_x_is_x_pred' * K_div_y_is);
      v_pred = mdivide_left_tri_low(L_sigmaEps, k_x_is_x_pred);
      cov_f_pred = cov_exp_quad(x_pred, sigma, rho) - v_pred' * v_pred;
      nug_pred = diag_matrix(rep_vector(1e-12,N_pred));
      f_pred = multi_normal_rng(f_pred_mu, cov_f_pred + nug_pred);
    }
    return f_pred;
    
    // {
    //   matrix[N, N] R = getR(x_is, rho);
    //   matrix[N, N] C = getC(R, sigma, sigmaEps);
    //   matrix[N, N] L_C = cholesky_decompose(C);
    //   
    //   
    //   vector[N] K_div_y_is_left = mdivide_left_tri_low(L_C, y_is);
    //   vector[N] K_div_y_is = mdivide_right_tri_low(K_div_y_is_left',L_C)';
    //   
    //   matrix[N, N_pred] k_x_is_x_pred;
    //   matrix[N, N_pred] v_pred;
    //   vector[N_pred] f_pred_mu;
    //   matrix[N_pred, N_pred] cov_f_pred;
    //   matrix[N_pred, N_pred] nug_pred;
    // 
    //   
    //   
    //   K_div_y_is = mdivide_right_tri_low(K_div_y_is',L_C)';
    //   k_x_is_x_pred = cov_exp_quad(x_is, x_pred, sigma, rho);
    //   f_pred_mu = (k_x_is_x_pred' * K_div_y_is);
    //   v_pred = mdivide_left_tri_low(L_sigmaEps, k_x_is_x_pred);
    //   cov_f_pred = cov_exp_quad(x_pred, sigma, rho) - v_pred' * v_pred;
    //   nug_pred = diag_matrix(rep_vector(1e-12,N_pred));
    //   f_pred = multi_normal_rng(f_pred_mu, cov_f_pred + nug_pred);
    // }
    // return f_pred;
  }

}
data {
  
  int<lower=1> N;
  int<lower=1> N_pred;
  vector[N] y;
  real x[N];
  real x_pred[N_pred];

}

parameters {
  
  real<lower=0> rho;
  real<lower=0> sigma;
  real<lower=0> sigmaEps;
  vector[N] eta;
}

transformed parameters {
  
  vector[N] f;
  
  {
    matrix[N, N] L;
    matrix[N, N] K;
    K = cov_exp_quad(x, sigma, rho);
    for (n in 1:N)
      K[n, n] = K[n, n] + 1e-12;
    L = cholesky_decompose(K);
    f = L * eta;
  }

}
model {

  rho ~ gamma(2, 20);
  sigma ~ normal(0, 1);
  sigmaEps ~ normal(0, 1);
  eta ~ normal(0, 1);
  y ~ normal(f, sigmaEps);

}

generated quantities {
  
  vector[N_pred] f_pred;
  vector[N_pred] y_pred;
  f_pred = gp_pred_rng(x_pred, y, x, sigma, rho, sigmaEps);
  for (n in 1:N_pred)
    y_pred[n] = normal_rng(f_pred[n], sigmaEps);

}

