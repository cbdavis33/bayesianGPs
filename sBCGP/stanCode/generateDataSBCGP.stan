functions{
  
  matrix getG(matrix x,
              vector rho) {
                
    int N = rows(x);
    int D = cols(x);
    vector[D] dist4;
    matrix[N, N] G;
    
    for (i in 1:(N-1)) {
      G[i, i] = 1.0;
      for (j in (i + 1):N) {
        dist4 = 4.0 * (x[i] - x[j])';
        G[i, j] = 1.0;
        for(k in 1:D){
          G[i,j] = G[i,j] * rho[k] ^ (dist4[k]^2);
        }
        G[j, i] = G[i, j];
      }
    }
    G[N, N] = 1;
    //return cholesky_decompose(R);
    return G;
  }
  
  matrix getR(real w,
              matrix G,
              matrix L){
   
    int N = rows(G); 
    matrix[N, N] R = w*G + (1 - w)*L;
    
    return R;
  }

  matrix getK(matrix R, 
              real sigma){
                
    int N = rows(R);
    matrix[N, N] K = sigma^2 * R;
  
    return K;
  }

  matrix getC(matrix x,
              real w,
              vector rhoG,
              vector rhoL,
              real sigma,
              real sigmaEps){
    
    int N = rows(x);
    matrix[N, N] G = getG(x, rhoG);
    matrix[N, N] L = getG(x, rhoL);
    matrix[N, N] R = getR(w, G, L);
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
  real<lower = 0.5, upper = 1> w;    // weight on global vs. local 
  vector<lower=0, upper=1>[D] rhoG;  // global correlation parameters for each dim
  vector<lower=0, upper=1>[D] rhoL;  // local correlation parameters for each dim
                                     // leave it on the user to ensure rho_L_i < rho_G_i
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
    matrix[n, n] G = getG(X, rhoG);
    matrix[n, n] L = getG(X, rhoL);
    matrix[n, n] R = getR(w, G, L);
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

