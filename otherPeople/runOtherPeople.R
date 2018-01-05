rm(list = ls())
cat("\014")

library(plyr)
library(tidyverse)
library(stringr)
library(rstan)
library(reshape2)
library(bayesplot)

# minor edit
# minor edit from cbdMerck
# setwd("~/Documents/bayesianGPs/")

rstan_options(auto_write = TRUE)
numCores <- 4                     # Find the number of cores on your machine
options(mc.cores = numCores)      # Then use one less than that for MCMC sampling

sim_data_model <- stan_model(file = 'stanCode/generateKrigingData.stan')

alpha <- 1
rho <- 0.1
sigma <- sqrt(0.1)


dat_list <- list(N = 100, alpha = alpha, length_scale = rho, sigma = sigma)
set <- sample(1:dat_list$N, size = 60, replace = F)
# draw <- sampling(sim_data_model, iter = 1, algorithm = 'Fixed_param', chains = 1, data = dat_list,
#                  seed = 363360090)
draw <- sampling(sim_data_model, iter = 1, algorithm = 'Fixed_param', chains = 1, data = dat_list)

samps <- rstan::extract(draw)
plt_df = with(samps,data.frame(x = x[1,], y = y[1,], f = f[1,]))

ggplot(data = plt_df[set,], aes(x=x, y=y)) +
  geom_point(aes(colour = 'Realized data')) +
  geom_line(data = plt_df, aes(x = x, y = f, colour = 'Latent mean function')) +
  theme_bw() + theme(legend.position="bottom") +
  scale_color_manual(name = '', values = c('Realized data'='black','Latent mean function'='red')) +
  xlab('X') +
  ylab('y') +
  ggtitle(str_c('N = ',length(set),' from length-scale = ', rho, ', alpha = ', alpha, ', sigma = ', round(sigma,2)))


stan_data <- list(N = length(set), N_pred = dat_list$N - length(set),
                  zeros = rep(0,length(set)), x = samps$x[1,set], y = samps$y[1,set],
                  x_pred = samps$x[1,-set], f_pred = samps$f[1,-set])

stanMLGP <- stan_model(file = 'otherPeople/stanUserGuide/margLikeGP.stan')
stanLVGP <- stan_model(file = 'otherPeople/stanUserGuide/latVarGP.stan')
trangucciMLGP <- stan_model(file = 'otherPeople/trangucciStanCon/margLikeGP.stan')
trangucciLVGP <- stan_model(file = 'otherPeople/trangucciStanCon/latVarGP.stan')


stanMLGPFit <- sampling(stanMLGP, 
                        data = stan_data, 
                        cores = 4, 
                        chains = 4, 
                        iter = 1000, 
                        control = list(adapt_delta = 0.999))

stanLVGPFit <- sampling(stanLVGP, 
                        data = stan_data, 
                        cores = 4, 
                        chains = 4, 
                        iter = 1000, 
                        control = list(adapt_delta = 0.999))

trangucciMLGPFit <- sampling(stanMLGP, 
                             data = stan_data, 
                             cores = 4, 
                             chains = 4, 
                             iter = 1000, 
                             control = list(adapt_delta = 0.999))

trangucciLVGPFit <- sampling(stanLVGP, 
                             data = stan_data, 
                             cores = 4, 
                             chains = 4, 
                             iter = 1000, 
                             control = list(adapt_delta = 0.999))


parsToFollow <- c("rho", "alpha", "sigma")
print(stanMLGPFit, pars = parsToFollow)
print(stanLVGPFit, pars = parsToFollow)
print(trangucciMLGPFit, pars = parsToFollow)
print(trangucciLVGPFit, pars = parsToFollow)
