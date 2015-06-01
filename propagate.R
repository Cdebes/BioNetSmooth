propagate_algebra <- function(W, z, alpha, iter)
{
  F <- z
  W_p <- W * alpha
  diag(W_p) <- 1 - alpha
  for (i in c(1:iter))
  {
    F <- (W_p %*% F)
  }
  F
}


