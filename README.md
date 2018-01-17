# Bayesian Gaussian Processes


## This Github Repository

This particular Github repo, *bayesianGPs* is a rough draft of putting together some software that will do the stuff in my dissertation, *A Bayesian Approach to Prediction adn Variable Selection Using Non-Stationary Gaussian Processes*. This one will use Stan, a Bayeian programming language/software. I might also duplicate in R using Rcpp what I did when I originally wrote this software in MATLAB, but that'll probably only be for random practice. 

I'll probably build up from the most basic stationary model and work towards the stuff in my dissertation. Eventually, I'll make an actual R package (Python, too?) in a different repo.

## A couple different Formulations

The marginal likelihood formulation is faster for Gaussian data, but can't be extended to non-Gaussian data. It seems like it does not need a nugget, as the residual error basically does the job of keeping the matrix positive definite.

The latent variable formulation is slower for Gaussian data, but can be extended to non-Gaussian data. It needs a nugget for positive definiteness of the covariance matrix.

## Setup

We have observed data, \mathbf{Y}. Conditional on model parameters, $\Lambda$, the data are treated as a realization of a Gaussian process, Y(x). $$ Y(\mathbf(x)) | \Lambda \sim GP\left( \mu, C\left(\cdot, \cdot \right) ) \right)$$

Let \mathbf{R} be an $n \times n$ correlation matrix for the observed data: 

## It turns out this README file doesn't support LaTeX. See bayesianGPs.html and/or bayesianGPs.pdf for more info.

