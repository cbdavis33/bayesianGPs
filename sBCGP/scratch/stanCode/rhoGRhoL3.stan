// functions{
//   
//     real truncated_beta_lpdf(real x1, 
//                              real x2,
//                              vector alpha,
//                              vector beta){
//       
//       real prob = x1^(alpha-1)*(x2 - x1)^(beta - 1)/x2^(alpha + beta - 1)
//       
//       real lprob = log(prob); 
//       
//     }
//   
// }

data{
  
  int<lower = 1> d;
  real<lower = 0> alphaW;
  real<lower = 0> betaW;
  vector<lower = 0>[d] alphaRhoG;
  vector<lower = 0>[d] betaRhoG;
  vector<lower = 0>[d] alphaRhoL;
  vector<lower = 0>[d] betaRhoL;
  
}

parameters{
  
  real<lower = 0, upper = 1> wRaw;
  vector<lower = 0, upper = 1>[d] rhoG;
  vector<lower = 0, upper = 1>[d] rhoLRaw;
  
}

transformed parameters{
  
  real<lower = 0, upper = 1> w = 0.5 + 0.5*wRaw;
  vector<lower = 0, upper = 1>[d] rhoL = rhoG .* rhoLRaw;
  
}

model{
  
  wRaw ~ beta(alphaW, betaW);
  rhoG ~ beta(alphaRhoG, betaRhoG);
  rhoLRaw ~ beta(alphaRhoL, betaRhoL);
  target += sum(log(rhoLRaw));
  
}


