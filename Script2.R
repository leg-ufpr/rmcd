# Script 2: Poisson-Tweedie distribution -------------------------------
# Author: Wagner Hugo Bonat LEG/UFPR -----------------------------------
# Date: 27/02/2017 -----------------------------------------------------

# Probability mass function --------------------------------------------
dptw <- function(y, mu, phi, power, control_sample) {
    pts <- control_sample$pts
    norma <- control_sample$norma
    integral <- mean(integrand(pts, y = y, mu = mu, phi = phi, 
                               power = power)/norma)
    return(integral)
}
dptw <- Vectorize(dptw, vectorize.args = "y")

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
plot_ptw <- function(mu, phi, power, title, n_sample = 100000, 
                     grid_y = 0:100, mc_sample = 5000) {
    control_sample <- list()
    if(power <= 2) {
    control_sample$pts <- rtweedie(n = mc_sample, mu = mu, 
                                   phi = phi, power = power)
    control_sample$norma <- dtweedie(control_sample$pts, mu = mu, 
                                     phi = phi, power = power)
    }
    if(power > 2) {
        control_sample$pts <- rinvgauss(n = mc_sample, mean = mu, 
                                       dispersion = phi)
        control_sample$norma <- dinvgauss(control_sample$pts, 
                                          mean = mu, dispersion = phi)
    }
    obs <- rptweedie(n = n_sample, mu = mu, phi = phi, power = power)
    dens <- dptw(y = grid_y, mu = mu, phi = phi, power = power, 
                 control_sample = control_sample)
    plot(table(obs)/n_sample, ylab = "Mass function", xlab = "y", col = "gray", 
         main = title)
    nn <- length(grid_y)-1
    lines(c(0:nn), dens, lty = 2, lwd = 2, col = "black", type = "l")
}