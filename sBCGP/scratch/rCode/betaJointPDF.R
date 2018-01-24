betaJointPDF <- function(alpha1, beta1, alpha2, beta2){
  
  x1 <- seq(0, 1, 0.1)
  x2 <- seq(0, 1, 0.1)
  
  myGrid <- expand.grid(x1,x2)
  names(myGrid) = c("x1", "x2")
  
  dens1 <- 1/beta(alpha2, beta2) * x2^(alpha2 - 1) * (x1 - x2)^(beta2 - 1) / 
    (x1^(alpha2 + beta2 - 1)) * (x2 < x1)
  
  dens2 <- dbeta(alpha1, beta1)
  
  density <- dens1 * dens2
  
  toReturn <- list(x1 = x1, x2 = x2, density = density)

  return(toReturn)
  
}

