# Script 2: Poisson-Tweedie distribution -------------------------------
# Author: Wagner Hugo Bonat LEG/UFPR -----------------------------------
# Date: 27/02/2017 -----------------------------------------------------

# Integrand ------------------------------------------------------------
integrand <- function(x, y, mu, phi, power) {
    int = dpois(y, lambda = x)*dtweedie(x, mu = mu, 
                                        phi = phi, power = power)
    return(int)
}

# Numerical integration using Gauss-Laguerre method --------------------
gauss_laguerre <- function(integrand, n_pts, y, mu, phi, power) {
    pts <- gauss.quad(n_pts, kind="laguerre")
    integral <- sum(pts$weights*integrand(pts$nodes, y = y, mu = mu, 
                                          phi = phi, power = power)/
                        exp(-pts$nodes))
    return(integral)
}

gauss_laguerre_vec <- Vectorize(FUN = gauss_laguerre, vectorize.args = "y")

# Numerical integration using Monte Carlo method -----------------------
monte_carlo <- function(integrand, n_pts, y, mu, phi, power) {
    pts <- rtweedie(n_pts, mu = mu, phi = phi, power = power)
    norma <- dtweedie(pts, mu = mu, phi = phi, power = power)
    integral <- mean(integrand(pts, y = y, mu = mu, phi = phi, 
                               power = power)/norma)
    return(integral)
}

# Probability mass function Poisson-Tweedie ----------------------------
dptweedie_aux <- function(y, mu, phi, power, n_pts, method) {
    if(method == "laguerre" | y > 0) {
        pmf <- gauss_laguerre(integrand = integrand, n_pts = n_pts, 
                              y = y, mu = mu, phi = phi, power = power)
    }
    if(method == "laguerre" & y == 0) {
        v.y <- round(10*sqrt(mu + phi*mu^power),0)
        if(v.y > 1000) {v.y <- 1000}
        #print(v.y)
        y.grid <- 0:v.y
        temp <- gauss_laguerre_vec(integrand = integrand, n_pts = n_pts,
                                   y = y.grid, mu = mu, phi = phi,
                                   power = power)
        pmf <- 1-sum(temp)+temp[1]
    }
    if(method == "MC") {
        pmf <- monte_carlo(integrand = integrand, n_pts = n_pts, 
                           y = y, mu = mu, phi = phi, power = power)
    }
    return(pmf)
}
# Vectorize version
dptw <- Vectorize(FUN = dptweedie_aux, vectorize.args = c("y","mu"))

# Simulating from Poisson-Tweedie models -------------------------------
rptweedie <- function(n, mu, phi, power) {
    if(power != 3) {
        x <- rtweedie(n = n, mu = mu, phi = phi, power = power)
    }
    if(power == 3) {
        x <- rinvgauss(n = n, mean = mu, dispersion = phi)
    }
    y <- rpois(n = n, lambda = x)
    return(y)
}

# Plot function --------------------------------------------------------
plot_ptw <- function(mu, phi, power, title, n_sample = 100000, n_pts = 25, method = "laguerre") {
    require(statmod)
    grid_y <- 0:100
    obs <- rptweedie(n = n_sample, mu = mu, phi = phi, power = power)
    dens <- dptw(y = grid_y, mu = mu, phi = phi, power = power, n_pts = n_pts, method = method)
    plot(table(obs)/n_sample, ylab = "Mass function", xlab = "y", col = "gray", 
         main = title)
    nn <- length(grid_y)-1
    lines(c(0:nn), dens, lty = 2, lwd = 2, col = "black", type = "l")
}

# Expectation and variance ---------------------------------------------
moments_ptw <- function(mu, phi, power) {
    expec <- mu
    variance <- mu + phi*(mu^power)
    return(c("Expectation" = expec, "Variance" = variance))
}
moments_ptw <- Vectorize(moments_ptw, c("mu"))

# Dispersion index ---------------------------------------------------
disp_index_ptw <- function(mu, phi, power) {
    TEMP <- moments_ptw(mu = mu, phi = phi, power = power)
    DI <- TEMP[2]/TEMP[1]
    return(DI)
}
disp_index_ptw <- Vectorize(disp_index_ptw, c("mu"))

# Zero-inflation index -----------------------------------------------
zero_inflation_ptw <- function(mu, phi, power) {
    PX0 <- dptw(y = 0, mu = mu, phi = phi, power = power, n_pts = 30, 
                method = "laguerre")
    mu <- moments_ptw(mu = mu, phi = phi, power = power)[1]
    ZI <- 1 + log(PX0)/mu
    return(ZI)
}
zero_inflation_ptw <- Vectorize(zero_inflation_ptw, c("mu"))

# Heavy-tail index ---------------------------------------------------
heavy_tail_ptw <- function(x, mu, phi, power) {
    Px1 <- dptw(y = c(x, x+1), mu = mu, phi = phi, power = power, 
               n_pts = 100, method = "laguerre")
    LTI <- Px1[2]/Px1[1]
    return(LTI)
}
heavy_tail_ptw <- Vectorize(heavy_tail_ptw, c("x"))

# END-------------------------------------------------------------------