# Script 1 - Gamma count distribution --------------------------------
# Author: Wagner Hugo Bonat LEG/UFPR ---------------------------------
# Date: 24/02/2017 ---------------------------------------------------

# Probability mass function ------------------------------------------
dgc <- function(y, beta, alpha) {
  p <- pgamma(q = 1,
              shape = y * alpha,
              rate = alpha * beta) -
    pgamma(q = 1,
           shape = (y + 1) * alpha,
           rate = alpha * beta)
  return(p)
}

# Mean and variance using numerical integration ----------------------
moments_gc <- function(beta, alpha) {
  E_gc <- function(y, beta, alpha) {
    EE <- y*dgc(y, beta, alpha)
    return(EE)
  }
  E2_gc <- function(y, beta, alpha) {
    E2 <- (y^2)*dgc(y, beta, alpha)
    return(E2)
  }
  Exp <- integrate(E_gc, lower = 0, upper = Inf, beta = beta, alpha = alpha)
  Exp2 <- integrate(E2_gc, lower = 0, upper = Inf, beta = beta, alpha = alpha)
  VV <- Exp2$value - Exp$value^2
  return(c("Expectation" = Exp$value, "Variance" = VV))
}
moments_gc <- Vectorize(moments_gc, c("beta"))

# Dispersion index ---------------------------------------------------
disp_index_gc <- function(beta, alpha) {
 TEMP <- moments_gc(beta = beta, alpha = alpha)
 DI <- TEMP[2]/TEMP[1]
 return(DI)
}
disp_index_gc <- Vectorize(disp_index_gc, c("beta"))

# Zero-inflation index -----------------------------------------------
zero_inflation_gc <- function(beta, alpha) {
  PX0 <- dgc(y = 0, beta = beta, alpha = alpha)
  mu <- moments_gc(beta = beta, alpha = alpha)[1]
  ZI <- 1 + log(PX0)/mu
  return(ZI)
}
zero_inflation_gc <- Vectorize(zero_inflation_gc, c("beta"))

# Heavy-tail index ---------------------------------------------------
heavy_tail_gc <- function(x, beta, alpha) {
  Px1 <- dgc(y = c(x, x+1), beta = beta, alpha = alpha)
  LTI <- Px1[2]/Px1[1]
  return(LTI)
}
heavy_tail_gc <- Vectorize(heavy_tail_gc, c("x"))

# END ------------------------------------------------------------------