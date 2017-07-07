# Script 3: COM-Poisson distribution -----------------------------------
# Author: Wagner Hugo Bonat LEG/UFPR -----------------------------------
# Date: 27/02/2017 -----------------------------------------------------

# Loading extra package
require(compoisson)

# Convergency of Z(lambda, nu) constant --------------------------------
# devtools::install_github("JrEduardo/cmpreg")
grid <- expand.grid(
    lambda = c(0.5, 1, 5, 10, 30, 50),
    nu = seq(0, 1, length.out = 11))
grid$z <- apply(grid, 1, function(x) {
    # com.compute.z(x[1], x[2]) # There is no a max iteration criterion
    cmpreg::computez(loglambda = log(x[1]), phi = log(x[2]),
                     tol = 1e-5, maxit = 1000)
})
format(xtabs(z ~ nu + lambda, data = grid), digits = 3)

# Computing mean and variance using numerical integration --------------
moments_cp <- function(lambda, nu, order = 1, upper = 2*lambda) {
    aux_fc <- function(x, lambda, nu, order) {
        output <- (x^order)*dcom(x, lambda = lambda, nu = nu)
        return(output)
    }
    output <- integrate(aux_fc, lower = 0, upper = upper, 
                        lambda = lambda, nu = nu, order = order)
    return(output$value)
}
moments_cp <- Vectorize(moments_cp, c("lambda"))

# Approximated mean and variance ---------------------------------------
ap_moments_cp <- function(lambda, nu) {
    expec <- lambda^(1/nu) - (nu - 1)/(2*nu)
    variance <- (1/nu)*(lambda^(1/nu))
    return(c("Expectation" = expec, "Variance" = variance))
}
ap_moments_cp <- Vectorize(ap_moments_cp, c("lambda"))

# Dispersion index -----------------------------------------------------
disp_index_cp <- function(lambda, nu) {
    TEMP <- ap_moments_cp(lambda = lambda, nu = nu)
    DI <- TEMP[2]/TEMP[1]
    return(DI)
}
disp_index_cp <- Vectorize(disp_index_cp, c("lambda"))

# Zero-inflation index -------------------------------------------------
zero_inflation_cp <- function(lambda, nu) {
    PX0 <- dcom(x = 0, lambda = lambda, nu = nu)
    mu <- ap_moments_cp(lambda = lambda, nu = nu)[1]
    ZI <- 1 + log(PX0)/mu
    return(ZI)
}
zero_inflation_cp <- Vectorize(zero_inflation_cp, c("lambda"))

# Heavy-tail index -----------------------------------------------------
heavy_tail_cp <- function(x, lambda, nu) {
    Px1 <- dcom(x = c(x, x+1), lambda = lambda, nu = nu)
    LTI <- Px1[2]/Px1[1]
    return(LTI)
}
heavy_tail_cp <- Vectorize(heavy_tail_cp, c("x"))

# Plot COM-Poisson mass function ---------------------------------------
plot_cmp <- function(lambda, nu, title, grid_y = 0:100) {
    dens <- dcom(x = grid_y, lambda = lambda, nu = nu)
    plot(dens, ylab = "Mass function", xlab = "y", type = "h",
         main = title)
}

# END-------------------------------------------------------------------
