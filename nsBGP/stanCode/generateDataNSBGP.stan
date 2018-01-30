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
              vector sig2X){

    int N = rows(R);
    vector[N] rootV = sqrt(sig2X);
    matrix[N, N] K = quad_form_diag(R, rootV);

    return K;
  }

  matrix getC(matrix x,
              vector rho,
              vector sig2X,
              real sigmaEps){
    
    int N = rows(x);
    matrix[N, N] R = getR(x, rho);
    matrix[N, N] K = getK(R, sig2X);
    
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
  real<lower=0> sigmaEps;            // Noise standard deviation
  real muV;                          // Mean of the log GP for the variance function
  vector<lower=0, upper=1>[D] rhoV;  // correlation parameters for log GP for the variance function
  real<lower = 0> sigmaV;            // Standard deviation for log GP for the variance function
  
}
transformed data {

  vector[n] muVec = rep_vector(mu, n);
  vector[n] muVVec = rep_vector(muV, n);

}
model {}
generated quantities {

  matrix[n, D] X;
  vector[n] y;
  vector[n] f;
  vector[n] logSig2X;
  vector[n] sig2X;
  for(i in 1:D){
    for (j in 1:n){
      X[j, i] = uniform_rng(0, 1);
    }
  }

  {
    matrix[n, n] R = getR(X, rho);
    matrix[n, n] RV = getR(X, rhoV);
    matrix[n, n] KV = sigmaV^2 * RV;
    matrix[n, n] K;
    matrix[n, n] L_KV;
    matrix[n, n] L_K;
    for (i in 1:n)
      KV[i, i] = KV[i, i] + 1e-12;
    L_KV = cholesky_decompose(KV);
    logSig2X = multi_normal_cholesky_rng(muVVec, L_KV);
    sig2X = exp(logSig2X);
    K = getK(R, sig2X);

    for (i in 1:n)
      K[i, i] = K[i, i] + 1e-12;
    L_K = cholesky_decompose(K);
    f = multi_normal_cholesky_rng(muVec, L_K);
  }

  for (j in 1:n)
    y[j] = normal_rng(f[j], sigmaEps);

}

