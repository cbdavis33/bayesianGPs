sumAlphaBeta <- function(alpha, beta){
  
  return(alpha + beta)
  
}

prodAlphaBeta <- function(alpha, beta){
  
  return(alpha * beta)
  
}

meanSDRhoL <- function(alphaRhoG, betaRhoG, alphaRhoL, betaRhoL){
  
  alphaPlusBetaL <- sumAlphaBeta(alphaRhoL, betaRhoL)
  alphaPlusBetaG <- sumAlphaBeta(alphaRhoG, betaRhoG)
  
  alphaTimesBetaL <- prodAlphaBeta(alphaRhoL, betaRhoL)
  alphaTimesBetaG <- prodAlphaBeta(alphaRhoG, betaRhoG)
  
  means <- (alphaRhoL/(alphaPlusBetaL)) * (alphaRhoG/(alphaPlusBetaG))
  
  vars <- alphaTimesBetaL/((alphaPlusBetaL)^2 * (alphaPlusBetaL + 1)) * 
    (alphaTimesBetaG/((alphaPlusBetaG)^2 * (alphaPlusBetaG + 1)) + (alphaRhoG/alphaPlusBetaG)^2) +
    ((alphaRhoL/alphaPlusBetaL)^2) * alphaTimesBetaG/((alphaPlusBetaG)^2 * (alphaPlusBetaG + 1))
    
  
  toReturn <- list(means = means, sds = sqrt(vars))
    
  return(toReturn)
  
}

meanSDRhoL(alphaRhoG, betaRhoG, alphaRhoL, betaRhoL)
