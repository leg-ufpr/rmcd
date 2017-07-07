# Script 4: Find parameters values -------------------------------------
# Author: Wagner Hugo Bonat LEG/UFPR -----------------------------------
# Date: 27/02/2017 -----------------------------------------------------

# Loading extra functions ----------------------------------------------
source("Script1.R")
source("Script3.R")

# System equation to find parameters to fixed Dispersion Index ---------
system_equation <- function(param, Expec, DI, moments_function) {
    par1 <- param[1]
    par2 <- param[2]
    print(param)
    moments <- moments_function(par1, par2)
    eq1 <- moments[1] -  Expec
    eq2 <- moments[2] / moments[1] - DI
    return(c(eq1, eq2))
}

# Solving the system using rootSolve package ---------------------------
require(rootSolve)

#-------------------------------------------
# For Gamma-Count model
out_gc1 <- multiroot(system_equation, start = c(10, 2), Expec = 10,
                     DI = 0.5, moments_function = moments_gc)
moments_gc(beta = out_gc1$root[1], alpha = out_gc1$root[2])

out_gc2 <- multiroot(system_equation, start = c(10, 0.5), Expec = 10,
                     DI = 2, moments_function = moments_gc)
moments_gc(beta = out_gc2$root[1], alpha = out_gc2$root[2])

out_gc3 <- multiroot(system_equation, start = c(10, 0.2), Expec = 10,
                     DI = 5, moments_function = moments_gc)
moments_gc(beta = out_gc3$root[1], alpha = out_gc3$root[2])

# # Could not find parameters to have Expectation = 10 and DI = 20
# out_gc4 <- multiroot(system_equation, start = c(10, 0.001),
#                      Expec = 10, DI = 20,
#                      moments_function = moments_gc)

#-------------------------------------------
# For COM-Poisson model
out_cp1 <- multiroot(system_equation, start = c(10, 2), Expec = 10,
                     DI = 0.5, moments_function = ap_moments_cp)
ap_moments_cp(lambda = out_cp1$root[1], nu = out_cp1$root[2])

out_cp2 <- multiroot(system_equation, start = c(10, 0.5), Expec = 10,
                     DI = 2, moments_function = ap_moments_cp)
ap_moments_cp(lambda = out_cp2$root[1], nu = out_cp2$root[2])

out_cp3 <- multiroot(system_equation, start = c(10, 0.2), Expec = 10,
                     DI = 5, moments_function = ap_moments_cp)
ap_moments_cp(lambda = out_cp3$root[1], nu = out_cp3$root[2])

# # Could not find parameters to have Expectation = 10 and DI = 20
# out_cp4 <- multiroot(system_equation, start = c(0.1, 0.01),
#                      Expec = 10, DI = 20,
#                      moments_function = ap_moments_cp)

#-------------------------------------------
# Compare parameters from each distribution combined with DI
pars <- sapply(list(out_gc1, out_gc2, out_gc3, out_cp1, out_cp2,
                    out_cp3), "[[", "root")
aux <- data.frame(
    model = rep(c("Gamma-Count", "COM-Poisson"), each = 6),
    param = rep(c("par1", "par2"), 6),
    DI = rep(rep(c(0.5, 2, 5), each = 2), 2),
    value = c(pars))
ftable(xtabs(value ~ param + model + DI, data = aux))

# END-------------------------------------------------------------------
