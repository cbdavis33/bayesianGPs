functions{
  
  matrix getR(matrix x,
              vector rho) {
                
    int N = rows(x);
    int D = cols(x);
    vector[D] dist4;
    matrix[N, N] R;
    
    for (i in 1:(N-1)) {
      R[i, i] = 1.0;
      for (j in (i + 1):N) {
        dist4 = 4.0 * (x[i] - x[j])';
        R[i, j] = 1.0;
        for(k in 1:D){
          R[i,j] = R[i,j] * rho[k] ^ (dist4[k]^2);
        }
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
  matrix getC(matrix x,
              vector rho,
              real sigma,
              real sigmaEps){
    
    int N = rows(x); 
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
  
  int<lower=1> n;                    // number of observations
  int<lower=1> D;                    // number of input dimensions
  real mu;                           // GP mean
  vector<lower=0, upper=1>[D] rho;   // correlation parameters for each dim
  real<lower=0> sigma;               // GP standard deviation
  real<lower=0> sigmaEps;            // Noise standard deviation
  
}
transformed data {

  vector[n] muVec = rep_vector(mu, n);

}
model {}
generated quantities {
  
  matrix[n, D] X;
  vector[n] y;
  vector[n] f;
  for(i in 1:D){   
    for (j in 1:n){
      X[j, i] = uniform_rng(0, 1);
    }
  }
  
  {
    matrix[n, n] R = getR(X, rho);
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

