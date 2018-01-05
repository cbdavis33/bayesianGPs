Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

C <- Posdef(n=5, ev=1:5)
eigen(C)$val

L_C <- t(chol(C))

L_C

L_C %*% t(L_C)
C
solve(t(L_C)) %*% solve(L_C)
solve(C)
solve(t(L_C)) %*% solve(L_C) - solve(C)