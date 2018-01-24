betaJointPDF <- function(alpha1, beta1, alpha2, beta2){
  
  x1 <- seq(0.01, 0.99, 0.01)
  x2 <- seq(0.01, 0.99, 0.01)
  
  myGrid <- expand.grid(x1,x2)
  names(myGrid) = c("x1", "x2")
  
  
  dens1 <- dbeta(myGrid$x1, alpha1, beta1)
  
  dens2 <- 1/beta(alpha2, beta2) * myGrid$x2^(alpha2 - 1) * (myGrid$x1 - myGrid$x2)^(beta2 - 1) / 
    (myGrid$x1^(alpha2 + beta2 - 1)) * (myGrid$x2 < myGrid$x1)
  
  density <- dens1 * dens2
  
  toReturn <- list(x1 = myGrid$x1, x2 = myGrid$x2, density = density)

  return(toReturn)
  
}

alphaRhoG <- 1
betaRhoG <- 1
alphaRhoL <- 1
betaRhoL <- 1


density <- betaJointPDF(alphaRhoG, betaRhoG, alphaRhoL, betaRhoL)
density <- data.frame(rhoG = density$x1, rhoL = density$x2, density = density$density)
head(density)

ggplot(data = density, aes(x = rhoG, y = rhoL, z = density)) + 
  geom_contour() +
  geom_abline(slope = 1, intercept = 0) +
  xlim(0,1) + 
  ylim(0,1)

###########################################################


tmp <- mesh(seq(0.01, 0.99, 0.01), seq(0.01, 0.99, 0.01))
rhoG <- tmp$x
rhoL <- tmp$y


dens1 <- dbeta(rhoG, alphaRhoG, betaRhoG)

dens2 <- 1/beta(alphaRhoL, betaRhoL) * rhoL^(alphaRhoL - 1) * (rhoG - rhoL)^(betaRhoL - 1) / 
  (rhoG^(alphaRhoL + betaRhoL - 1)) * (rhoL < rhoG)

density <- dens1 * dens2

surf3D(rhoG, rhoL, density, theta = 40)

#######################################################
# 
# kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))
# p <- plotly::plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% plotly::add_surface()

rhoG <- seq(0.01, 0.99, 0.01)
rhoL <- seq(0.01, 0.99, 0.01)
density <- density
p <- plotly::plot_ly(x = rhoG, y = rhoL, z = density) %>% plotly::add_surface()
p
