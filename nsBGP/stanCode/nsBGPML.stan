functions {
  
  // This might need to be redone as an array or array of vectors (pg 323 in Stan User Guide)
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

  matrix getRop(matrix x,
                matrix xP,
                vector rho) {

    int n = rows(x);
    int nP = rows(xP);
    int d = cols(x);
    vector[d] dist4;

    matrix[n, nP] Rop;
    for(i in 1:nP){
      for(j in 1:n){
        dist4 = 4.0 * (xP[i] - x[j])';
        Rop[j, i] = 1.0;
        for(k in 1:d){
          Rop[j,i] = Rop[j,i] * rho[k] ^ (dist4[k]^2);
        }
      }
    }

    return Rop;

  }

  matrix getKop(matrix Rop,
                vector sig2X,
                vector sig2Xp){

    int n = rows(Rop);
    int nP = cols(Rop);

    vector[n] rootV = sqrt(sig2X);
    vector[nP] rootVp = sqrt(sig2Xp);
    matrix[n, nP] Kop = (rootV * rootVp') .* Rop;

    return Kop;
  }

  matrix getDop(matrix x,
                matrix xP){

    int n = rows(x);
    int nP = rows(xP);
    int d = cols(x);

    matrix[n, nP] Dop;
    for(i in 1:nP){
      for(j in 1:n){
        Dop[j,i] = 1;
        for(k in 1:d){
          Dop[j,i] = Dop[j,i] * (x[j,k] == xP[i,k]);
        }
      }
    }

    return Dop;

  }

  matrix getCop(matrix x,
                matrix xP,
                vector rho,
                vector sig2X,
                vector sig2Xp,
                real sigmaEps){

    int n = rows(x);
    int nP = rows(xP); 

    matrix[n, nP] Rop = getRop(x, xP, rho);
    matrix[n, nP] Kop = getKop(Rop, sig2X, sig2Xp);
    matrix[n, nP] Dop = getDop(x, xP);

    matrix[n, nP] Cop = Kop + sigmaEps^2 * Dop;

    return Cop;

  }

  vector gp_f_rng(matrix xPred,
                  vector yObs,
                  matrix xObs,
                  real mu,
                  vector rho,
                  vector sig2X,
                  vector sig2Xp,
                  real sigmaEps) {

    int nPred = rows(xPred);  
    int n = rows(yObs);

    vector[nPred] fPred;
    // vector[nPred] yPred;

    {
      // matrix[n, n] R = getR(xObs, rho);
      // matrix[n, n] C = getC(R, sigma, sigmaEps);
      matrix[n, n] C = getC(xObs, rho, sig2X, sigmaEps);
      matrix[n, n] L_C = cholesky_decompose(C);
      
      matrix[nPred, nPred] Rp = getR(xPred, rho);
      matrix[nPred, nPred] Kp = getK(Rp, sig2Xp) + diag_matrix(rep_vector(1e-12, nPred));
      // matrix[nPred, nPred] Cp = getC(xPred, rho, sigma, sigmaEps);

      matrix[n, nPred] Rop = getRop(xObs, xPred, rho);
      matrix[n, nPred] Kop = getKop(Rop, sig2X, sig2Xp);
      // matrix[n, nPred] Cop = getCop(xObs, xPred, rho, sigma, sigmaEps);

      matrix[n, nPred] L_C_Inv_Kop = mdivide_left_tri_low(L_C, Kop);
      vector[n] yMinusMu = yObs - mu;
      vector[n] L_C_Inv_yMinusMu = mdivide_left_tri_low(L_C, yMinusMu);

      vector[nPred] meanPred = mu + L_C_Inv_Kop' * L_C_Inv_yMinusMu;
      matrix[nPred, nPred] varPred = Kp - crossprod(L_C_Inv_Kop) + // L_C_Inv_Kop' * L_C_Inv_Kop +
                                     diag_matrix(rep_vector(1e-9, nPred));

      // meanPred = mu + Kop' * inverse_spd(C) * (yObs - mu);
      // varPred = Kp - Kop' * inverse_spd(C) * Kop;

      fPred = multi_normal_rng(meanPred, varPred);

    }

    return fPred;

  }

  // NEEDS MODIFICATION AFTER REWRITING getKop
  vector gp_y_rng(matrix xPred,
                  vector yObs,
                  matrix xObs,
                  real mu,
                  vector rho,
                  vector sig2X,
                  vector sig2Xp,
                  real sigmaEps) {

    int nPred = rows(xPred);  // this might need to be changed to rows(xPred) when d > 1
    int n = rows(yObs);

    vector[nPred] yPred;

    {
      // matrix[n, n] R = getR(xObs, rho);
      // matrix[n, n] C = getC(R, sigma, sigmaEps);
      matrix[n, n] C = getC(xObs, rho, sig2X, sigmaEps);
      matrix[n, n] L_C = cholesky_decompose(C);

      // matrix[nPred, nPred] Rp = getR(x_pred, rho);
      // matrix[nPred, nPred] Cp = getC(Rp, sigma, sigmaEps);
      matrix[nPred, nPred] Cp = getC(xPred, rho, sig2Xp, sigmaEps);

      matrix[n, nPred] Cop = getCop(xObs, xPred, rho, sig2X, sig2Xp, sigmaEps);

      matrix[n, nPred] L_C_Inv_Cop = mdivide_left_tri_low(L_C, Cop);
      vector[n] yMinusMu = yObs - mu;
      vector[n] L_C_Inv_yMinusMu = mdivide_left_tri_low(L_C, yMinusMu);

      vector[nPred] meanPred = mu + L_C_Inv_Cop' * L_C_Inv_yMinusMu;
      matrix[nPred, nPred] varPred = Cp - crossprod(L_C_Inv_Cop); // L_C_Inv_Cop' * L_C_Inv_Cop;

      // meanPred = mu + Cop' * inverse_spd(C) * (yObs - mu);
      // varPred = Cp - Cop' * inverse_spd(C) * Cop;

      yPred = multi_normal_rng(meanPred, varPred);

    }

    return yPred;

  }

  vector gp_V_rng(matrix xPred,
                  vector logVObs,
                  matrix xObs,
                  real muV,
                  real sigmaV,
                  vector rhoV) {

    int nPred = rows(xPred);  
    int n = rows(logVObs);

    vector[nPred] vPred;

    {
      
      matrix[nPred, nPred] KV = sigmaV^2 * getR(xPred, rhoV)  + diag_matrix(rep_vector(1e-9, nPred));
      matrix[n, n] K = sigmaV^2 * getR(xObs, rhoV)  + diag_matrix(rep_vector(1e-9, n));
      matrix[n, n] L_K = cholesky_decompose(K);
      matrix[n, nPred] Kop = sigmaV^2 * getRop(xObs, xPred, rhoV);   
      
      matrix[n, nPred] L_K_Inv_Kop = mdivide_left_tri_low(L_K, Kop);
      vector[n] yMinusMu = logVObs - muV;
      vector[n] L_K_Inv_yMinusMu = mdivide_left_tri_low(L_K, yMinusMu);
      
      vector[nPred] meanVPred = muV + L_K_Inv_Kop' * L_K_Inv_yMinusMu;
      matrix[nPred, nPred] varPred = KV - crossprod(L_K_Inv_Kop); // L_K_Inv_Kop' * L_K_Inv_Kop;
      
      matrix[nPred, nPred] L_varPred = cholesky_decompose(varPred);
      
      vector[nPred] logVPred = multi_normal_cholesky_rng(meanVPred, L_varPred);
      
      vPred = exp(logVPred);

    }

    return vPred;

  }

  
}
data {
  
  int<lower=1> n;
  int<lower=1> d;
  int<lower=1> nPred;
  vector[n] y;
  matrix[n,d] x;
  matrix[nPred,d] xPred;
  
}

parameters {
  
  real mu;
  vector<lower = 0, upper = 1>[d] rho;
  real<lower=0> sigmaEps;
  vector[n] logSig2X;
  real muV;
  vector<lower = 0, upper = 1>[d] rhoV;
  real<lower=0> sigmaV;
  
}

transformed parameters {
  
  vector<lower = 0>[n] sig2X = exp(logSig2X); 
  
}
model {
  
  matrix[n, n] L_C;
  vector[n] muVec;
  
  matrix[n, n] L_KV;
  vector[n] muVVec;
  
  {
   
    matrix[n, n] C = getC(x, rho, sig2X, sigmaEps);
    matrix[n, n] KV = sigmaV^2 * getR(x, rhoV) + diag_matrix(rep_vector(1e-9, n));
    L_C = cholesky_decompose(C);
    L_KV = cholesky_decompose(KV);
    muVec = rep_vector(mu, n);
    muVVec = rep_vector(muV, n);

  }
  
  mu ~ normal(0, 1e6);
  rho ~ beta(1, 1);
  sigmaEps ~ normal(0, 1);
  muV ~ normal(0, 1);
  rhoV ~ beta(1,1);
  sigmaV ~ normal(0, 1);
  logSig2X ~ multi_normal_cholesky(muVVec, L_KV);
  
  y ~ multi_normal_cholesky(muVec, L_C);
  
}

generated quantities {

  vector[nPred] fPred;
  vector[nPred] yPred;
  vector[nPred] yPred2;
  vector[nPred] sig2Xp;

  {
    
    sig2Xp = gp_V_rng(xPred, logSig2X, x, muV, sigmaV, rhoV);
    
  }

  fPred = gp_f_rng(xPred, y, x, mu, rho, sig2X, sig2Xp, sigmaEps);
  for (i in 1:nPred)
    yPred[i] = normal_rng(fPred[i], sigmaEps);
  yPred2 = gp_y_rng(xPred, y, x, mu, rho, sig2X, sig2Xp, sigmaEps);


}

