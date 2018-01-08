functions {
  
  // Look at pg. 256 in Stan Version 2.17.0 User Guide for a reference on 
  // creating this R matrix
  
  // getR for d = 1
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
  
  // getRop for d = 1
  matrix getRop(real[] x,
                real[] xP,
                real rho) {
                
    int n = size(x);
    int nP = size(xP);
    
    matrix[n, nP] Rop;
    for(i in 1:nP){
      for(j in 1:n){
        Rop[j,i] = rho^(16 * (x[j] - xP[i])^2);
      }
    }
    
    return Rop;
  }

  matrix getKop(matrix Rop, 
                real sigma){
                
    int n = rows(Rop);
    int nP = cols(Rop);
    
    matrix[n, nP] Kop = sigma^2 * Rop;
  
    return Kop;
  }
  
  matrix getDop(real[] x,
                real[] xP){
    
    int n = size(x);
    int nP = size(xP);
    
    matrix[n, nP] Dop;
    for(i in 1:nP){
      for(j in 1:n){
        Dop[j,i] = (x[j] == xP[i]);
      }
    }
    
    return Dop;
    
  }

  // fine for d = 1. Will need some modifications for d > 1
  matrix getCop(real[] x,
                real[] xP,
                real rho,
                real sigma,
                real sigmaEps){
    
    int n = size(x);
    int nP = size(xP); //maybe will be rows(x) if d > 1 (x will be a matrix or vector of arrays)
    
    matrix[n, nP] Rop = getRop(x, xP, rho);
    matrix[n, nP] Kop = getKop(Rop, sigma);
    matrix[n, nP] Dop = getDop(x, xP);
    matrix[n, nP] cop = Kop + Dop;
    
    matrix[n, nP] Cop = Kop + sigmaEps^2 * Dop;
    
    return Cop;              
  
  }

  // // getR for d > 1
  // matrix getR(real[] x,
  //             vector rho) {
  //               
  //   int N = size(x);
  //   matrix[N, N] R;
  //   for (i in 1:(N-1)) {
  //     R[i, i] = 1;
  //     for (j in (i + 1):N) {
  //       R[i, j] = rho^(16 * (x[i] - x[j])^2);
  //       R[j, i] = R[i, j];
  //     }
  //   }
  //   R[N, N] = 1;
  //   //return cholesky_decompose(R);
  //   return R;
  // }
  
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
  


  // for d = 1
  vector gp_pred_rng(real[] xPred,
                     vector yObs,
                     real[] xObs,
                     real mu,
                     real rho,
                     real sigma,
                     real sigmaEps) {
    
    int nPred = size(xPred);  // this might need to be changed to rows(xPred) when d > 1                
    vector[nPred] fPred;
    int n = rows(yObs);
    vector[nPred] yPred;
    
    {
      // matrix[n, n] R = getR(xObs, rho);
      // matrix[n, n] C = getC(R, sigma, sigmaEps);
      matrix[n, n] C = getC(xObs, rho, sigma, sigmaEps);
      matrix[n, n] L_C = cholesky_decompose(C);
      
      // matrix[nPred, nPred] Rp = getR(x_pred, rho);
      // matrix[nPred, nPred] Cp = getC(Rp, sigma, sigmaEps);
      matrix[nPred, nPred] Cp = getC(xPred, rho, sigma, sigmaEps);

      matrix[n, nPred] Cop = getCop(xObs, xPred, rho, sigma, sigmaEps);
      
      matrix[n, nPred] L_C_Inv_Cop = mdivide_left_tri_low(L_C, Cop);
      vector[n] yMinusMu = yObs - mu;
      vector[n] L_C_Inv_yMinusMu = mdivide_left_tri_low(L_C, yMinusMu);
      
      vector[nPred] meanPred = mu + L_C_Inv_Cop * L_C_Inv_yMinusMu;
      matrix[nPred, nPred] varPred = Cp - L_C_Inv_Cop' * L_C_Inv_Cop;
      
      // meanPred = mu + Cop' * inverse_spd(C) * (yObs - mu);
      // varPred = Cp - Cop' * inverse_spd(C) * Cop;
      
      yPred = multi_normal_rng(meanPred, varPred);
      
    }
    
    return yPred;
    
  }

}
data {
  
  int<lower=1> n;
  int<lower=1> nPred;
  vector[n] y;
  real x[n];
  real xPred[nPred];

}

parameters {
  
  real mu;
  real<lower=0> rho;
  real<lower=0> sigma;
  real<lower=0> sigmaEps;
  vector[n] eta;
  
}

transformed parameters {
  
  vector[n] f;
  
  {
    matrix[n, n] R = getR(x, rho);
    matrix[n, n] K = getK(R, sigma);
    matrix[n, n] L_K;
    for (i in 1:n)
      K[i, i] = K[i, i] + 1e-12;
    L_K = cholesky_decompose(K);
  
    f = mu + L_K * eta;
  }

}
model {

  mu ~ normal(0, 1e6);
  rho ~ uniform(0, 1);
  sigma ~ normal(0, 1);
  sigmaEps ~ normal(0, 1);
  eta ~ normal(0, 1);
  y ~ normal(f, sigmaEps);

}

// generated quantities {
//   
//   vector[nPred] f_pred;
//   vector[nPred] y_pred;
//   
//   // f_pred = gp_pred_rng(x_pred, y, x, sigma, rho, sigmaEps);
//   // for (n in 1:N_pred)
//   //   y_pred[n] = normal_rng(f_pred[n], sigmaEps);
// 
// }

