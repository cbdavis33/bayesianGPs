data{
  int<lower = 1> d;
  vector<lower = 0>[d] alphaRhoG;
  vector<lower = 0>[d] betaRhoG;
  vector<lower = 0>[d] alphaRhoL;
  vector<lower = 0>[d] betaRhoL;
}
parameters{
  
  vector<lower = 0, upper = 1>[d] rhoG; 
  vector<lower = 0, upper = 1>[d] rhoL; 
  
}
transformed parameters{
  
}
model{
  
}