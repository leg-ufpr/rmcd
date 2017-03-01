# Script 3: COM-Poisson distribution -----------------------------------
# Author: Wagner Hugo Bonat LEG/UFPR -----------------------------------
# Date: 27/02/2017 -----------------------------------------------------

# Loading extra package
require(compoisson)

# Computing mean and variance using numerical integration --------------
moments_cp <- function(lambda, nu, order = 1) {
    aux_fc <- function(x, lambda, nu, order) {
        output <- (x^order)*dcom(x, lambda = lambda, nu = nu)
        return(output)
    }
    output <- integrate(aux_fc, lower = 0, upper = lambda*2, 
                        lambda = lambda, nu = nu, order = order)
    return(output)
}

# Approximated mean and variance ---------------------------------------
ap_moments_cp <- function(lambda, nu) {
    expec <- lambda^(1/nu) - (nu - 1)/(2*nu)
    variance <- (1/nu)*(lambda^(1/nu))
    return(c(expec, variance))
}

di_cmp <- function(lambda, nu) {
    tt <- ap_moments_cp(lambda = lambda, nu = nu)
    return(tt[2]/tt[1])
}
# Plot COM-Poisson mass function ---------------------------------------
plot_cmp <- function(lambda, nu, title, grid_y = 0:100) {
    dens <- dcom(x = grid_y, lambda = lambda, nu = nu)
    plot(dens, ylab = "Mass function", xlab = "y", type = "h",
         main = title)
}

# Find parameter values ------------------------------------------------
system_equation <- function(param, Expec, DI) {
    lambda <- param[1]
    nu <- param[2]
    eq1 <- ap_moments_cp(lambda, nu)[1] -  Expec
    eq2 <- di_cmp(lambda, nu) - DI
    return(c(eq1, eq2))
}

# Using rootSolve package
#require(rootSolve)
#multiroot(system_equation, Expec = 10, DI = 0.5, start = c(10, 0.5))
#ap_moments_cp(lambda = 118.51, nu = 2.05)

#multiroot(system_equation, Expec = 10, DI = 2, start = c(2, 0.5))
#ap_moments_cp(lambda = 2.88, nu = 0.47)

#multiroot(system_equation, Expec = 10, DI = 10, start = c(1, 0.01))
#ap_moments_cp(lambda = 1.30, nu = 0.13)

#multiroot(system_equation, Expec = 10, DI = 6, start = c(1.2, 0.13))

