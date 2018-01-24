rm(list = ls())
cat("\014")

library(plyr)
library(tidyverse)
library(stringr)
library(rstan)
library(reshape2)
library(bayesplot)

rstan_options(auto_write = TRUE)
numCores <- 4                     # Find the number of cores on your machine
options(mc.cores = numCores)      # Then use one less than that for MCMC sampling

source("sBCGP/scratch/rCode/betaMeanVar.R")

alphaW <- 6
betaW <- 4

# alphaRhoG <- c(1, 1, 1)
# betaRhoG <- c(1, 1, 1)
# alphaRhoL <- c(1, 8, 3)
# betaRhoL <- c(1, 2, 6)

alphaRhoG <- 0.5
betaRhoG <- 0.5
alphaRhoL <- 0.5
betaRhoL <- 0.5

D <- length(alphaRhoG)

stanData <- list(d = D, alphaW = alphaW, betaW = betaW, alphaRhoG = as.array(alphaRhoG), 
                 betaRhoG =  as.array(betaRhoG),
                 alphaRhoL = as.array(alphaRhoL), betaRhoL = as.array(betaRhoL))

# stanModel1 <- stan_model(file = 'sBCGP/scratch/stanCode/rhoGRhoL1.stan')
# stanModel2 <- stan_model(file = 'sBCGP/scratch/stanCode/rhoGRhoL2.stan')
# stanModel3 <- stan_model(file = 'sBCGP/scratch/stanCode/rhoGRhoL3.stan')

samples1 <- sampling(stanModel1, 
                     data = stanData, 
                     cores = 1, 
                     chains = 4, 
                     iter = 3000, 
                     control = list(adapt_delta = 0.90))

samples2 <- sampling(stanModel2, 
                     data = stanData, 
                     cores = 1, 
                     chains = 4, 
                     iter = 3000, 
                     control = list(adapt_delta = 0.90))

samples3 <- sampling(stanModel3, 
                     data = stanData, 
                     cores = 1, 
                     chains = 4, 
                     iter = 3000, 
                     control = list(adapt_delta = 0.90))


parsToMonitor <- c("w", "wRaw", "rhoG", "rhoL", "rhoLRaw")
print(samples1, pars = parsToMonitor, digits = 4)
print(samples2, pars = parsToMonitor, digits = 4)
print(samples3, pars = parsToMonitor, digits = 4)

meanSDRhoL(alphaRhoG, betaRhoG, alphaRhoL, betaRhoL)

samps1 <- extract(samples1)
samps2 <- extract(samples2)
samps3 <- extract(samples3)

samps1Mat <- as.matrix(samples1)
samps2Mat <- as.matrix(samples2)
samps3Mat <- as.matrix(samples3)

rhoG1 <- samps1Mat[,"rhoG[1]"]
rhoL1 <- samps1Mat[,"rhoL[1]"]

rhoG2 <- samps2Mat[,"rhoG[1]"]
rhoL2 <- samps2Mat[,"rhoL[1]"]

rhoG3 <- samps3Mat[,"rhoG[1]"]
rhoL3 <- samps3Mat[,"rhoL[1]"]

contourPlot1 <- ggplot(data = data.frame(rhoG = rhoG1, rhoL = rhoL1), aes(x = rhoG, y = rhoL)) +
  geom_density2d(aes(color = ..level..)) +
  theme_classic() +
  xlim(0, 1) + 
  ylim(0, 1) +
  ggtitle("No Gradient Adjustment")
contourPlot1

contourPlot2 <- ggplot(data = data.frame(rhoG = rhoG2, rhoL = rhoL2), aes(x = rhoG, y = rhoL)) +
  geom_density2d(aes(color = ..level..)) +
  theme_classic() +
  xlim(0, 1) + 
  ylim(0, 1) +
  ggtitle("With Gradient Adjustment")
contourPlot2

contourPlot3 <- ggplot(data = data.frame(rhoG = rhoG3, rhoL = rhoL3), aes(x = rhoG, y = rhoL)) +
  geom_density2d(aes(color = ..level..)) +
  theme_classic() +
  xlim(0, 1) + 
  ylim(0, 1) +
  ggtitle("With Gradient Adjustment RhoLRaw")
contourPlot3


wRaw <- seq(0, 1, .005)
fwRaw <- dbeta(wRaw, alphaW, betaW)

par(mfrow = c(3,1))
hist(samps1$w)
hist(samps1$wRaw)
plot(w, fw, type = 'l')

par(mfrow = c(3,1))
hist(samps2$w)
hist(samps2$wRaw)
plot(w, fw, type = 'l')

print(round(alphaW/(alphaW + betaW),4))
print(round(sqrt(alphaW*betaW/((alphaW + betaW)^2 * (alphaW + betaW + 1))), 4))
      