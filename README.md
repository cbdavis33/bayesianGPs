# bayesianGPs
Bayesian Gaussian Processes

The marginal likelihood formulation is faster for Gaussian data, but can't be extended to non-Gaussian data. It seems like it does not need a nugget, as the residual error basically does the job of keeping the matrix positive definite.

The latent variable formulation is slower for Gaussian data, but can be extended to non-Gaussian data. It needs a nugget for positive definiteness of the covariance matrix.

## Setup

We have observed data, \mathbf{Y}. Conditional on model parameters, $\Lambda$, the data are treated as a realization of a Gaussian process, Y(x). $$ Y(\mathbf(x)) | \Lambda \sim GP\left( \mu, C\left(\cdot, \cdot \right) ) \right)$$

Let \mathbf{R} be an $n \times n$ correlation matrix for the observed data: $$ R_{ij} = \rho^{16} $$

