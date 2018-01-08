functions{
  
  matrix getR(real[] x,
              real rho) {
                
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

  matrix getK(matrix R, 
              real sigma){
                
    int N = rows(R);
    matrix[N, N] K = sigma^2 * R;
  
    return K;
  }

  // fine for d = 1. Will need some modifications for d > 1
  matrix getC(real[] x,
              real rho,
              real sigma,
              real sigmaEps){
    
    int N = size(x); //maybe will be rows(x) if d > 1 (x will be a matrix or vector of arrays)
    matrix[N, N] R = getR(x, rho);
    matrix[N, N] K = getK(R, sigma);
    
    matrix[N, N] C = K + diag_matrix(rep_vector(sigmaEps^2,N));
    // 
    // matrix[N, N] C = sigma^2 * R; // This is matrix K in my handwritten work
    // 
    // for (n in 1:N)
    //     C[n, n] = C[n, n] + sigmaEps^2;
    // // C += diag_matrix(rep_vector(sigmaEps^2,N)); // Is looping or adding sig2eps*I faster?
    
    return C;              
  
  }
}
data {
  
  int<lower=1> n;
  real mu;
  real<lower=0> rho;
  real<lower=0> sigma;
  real<lower=0> sigmaEps;
  
}
transformed data {

  vector[n] muVec = rep_vector(mu, n);

}
model {}
generated quantities {
  
  real x[n];
  vector[n] y;
  vector[n] f;
  for (i in 1:n)
    x[i] = uniform_rng(0,1);
  
  {
    matrix[n, n] R = getR(x, rho);
    matrix[n, n] K = getK(R, sigma);
    matrix[n, n] L_K;
    for (i in 1:n)
      K[i, i] = K[i, i] + 1e-12;
    L_K = cholesky_decompose(K);
    f = multi_normal_cholesky_rng(muVec, L_K);
  }
  
  for (j in 1:n)
    y[j] = normal_rng(f[j], sigmaEps);
  
}

