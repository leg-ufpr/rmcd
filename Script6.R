# Script 6: Ortogonality study on count distributions ------------------
# Author: Eduardo Elias Ribeiro Junior LEG/UFPR ------------------------
# Date: 24/03/2017 -----------------------------------------------------

# Packages and functions -----------------------------------------------
library(MRDCr)
library(gridExtra)
library(parallel)

# Trellis layout
source("config/_setup.R")
ac <- list(pad1 = 0.5, pad2 = 0.5, tck = 0.5)
ps <- list(
    layout.widths = list(
        left.padding = 0.25,
        right.padding = -1,
        ylab.axis.padding = 0),
    layout.heights = list(
        bottom.padding = 0.25,
        top.padding = 0,
        axis.xlab.padding = 0,
        xlab.top = 0),
    axis.components = list(
        bottom = ac, top = ac,
        left = ac, right = ac)
)

# To simule counts following dispersion indexes
simule <- function(config, n, y, seed = NULL) {
    lapply(pars, function(x) {
        if (!is.null(seed)) set.seed(seed)
        if (any(names(x) %in% "power")) {
            sim <- do.call(x$fun, c(list(n), x[-1]))
        } else {
            prob <- do.call(x$fun, c(list(y = y), x[-1]))
            sim <- sample(y, size = n, replace = TRUE, prob = prob)
        }
        return(sim)
    })
}

# Get n characters
substrlast <- function(string, n = 2){
    substr(string, nchar(string) - n + 1, nchar(string))
}

# Simulate and organize data -------------------------------------------

# For PTW distribution

# Configuration
source("Script4.R")
pars <- list(
    GC_.5 = c(fun = dgcnt, lambda = aux[1, 4], alpha = aux[2, 4]),
    GC_02 = c(fun = dgcnt, lambda = aux[3, 4], alpha = aux[4, 4]),
    GC_05 = c(fun = dgcnt, lambda = aux[5, 4], alpha = aux[6, 4]),
    CP_.5 = c(fun = dcmp, sumto = 50,
              lambda = aux[7, 4], nu = aux[8, 4]),
    CP_02 = c(fun = dcmp, sumto = 80,
              lambda = aux[9, 4], nu = aux[10, 4]),
    CP_05 = c(fun = dcmp, sumto = 150,
              lambda = aux[11, 4], nu = aux[12, 4])
)

# Simulation process
samples <- simule(pars, n = 5000, y = 0:120, seed = 2017)

# Verify simulations
do.call("rbind",
    lapply(samples, function(x) {
        c("Mean" = mean(x), "Var" = var(x), "DI" = var(x)/mean(x))
    })
)

# info <- plyr::ldply(strsplit(rep(names(samples), each = n), "_"))
# info[, 2] <- as.numeric(info[, 2])
# colnames(info) <- c("model", "di")
# dasim <- cbind(info, y = unlist(samples))
# rownames(dasim) <- NULL

##----------------------------------------------------------------------
## Compute deviances models

##-------------------------------------------
## Gamma-count
indGC <- grep("GC", names(samples))
fitsGC <- mclapply(indGC, function(x) {
    ##-------------------------------------------
    ## Fit model
    data <- data.frame(y = samples[[x]])
    start <- log(rev(unlist(pars[[x]][-1])))
    names(start) <- c("lalpha", "lgamma")
    m0 <- gcnt(y ~ 1, data = data, start = start)
    ##-------------------------------------------
    ## Define grid
    interval <- confint(m0, level = 0.999, method = "quad")
    aux <- lapply(as.data.frame(t(interval)),
                  function(x) seq(x[1], x[2], length.out = 50))
    grid <- do.call("expand.grid", list(aux, KEEP.OUT.ATTRS = FALSE))
    ##-------------------------------------------
    ## Compute loglikelihood and deviance for grid
    grid$loglik <- apply(grid, 1, function(par) {
        -llgcnt(params = par, y = data$y, X = m0@data$X)
    })
    grid$deviance <- -2 * (grid$loglik - logLik(m0))
    ##-------------------------------------------
    return(list(model = m0, grid = grid))
}, mc.preschedule = FALSE, mc.cores = length(indGC))

##-------------------------------------------
## COM-Poisson
indCP <- grep("CP", names(samples))
fitsCP <- mclapply(indCP, function(x) {
    ##-------------------------------------------
    ## Fit model
    data <- data.frame(y = samples[[x]])
    start <- log(rev(unlist(pars[[x]][-(1:2)])))
    sumto <- pars[[x]]$sumto
    names(start) <- c("lnu", "llambda")
    m0 <- cmp(y ~ 1, data = data, start = start, sumto = sumto)
    ##-------------------------------------------
    ## Define grid
    interval <- confint(m0, level = 0.999, method = "quad")
    aux <- lapply(as.data.frame(t(interval)),
                  function(x) seq(x[1], x[2], length.out = 50))
    grid <- do.call("expand.grid", list(aux, KEEP.OUT.ATTRS = FALSE))
    ##-------------------------------------------
    ## Compute loglikelihood and deviance for grid
    grid$loglik <- apply(grid, 1, function(par) {
        -llcmp(params = par, y = data$y, X = m0@data$X, sumto = sumto)
    })
    grid$deviance <- -2 * (grid$loglik - logLik(m0))
    ##-------------------------------------------
    return(list(model = m0, grid = grid))
}, mc.preschedule = FALSE, mc.cores = length(indCP))

##======================================================================
## Describe the results

##-------------------------------------------
## Deviance contours
plotsGC <- lapply(seq(fitsGC), function(i) {
    ##-------------------------------------------
    main <- paste("Gamma-Count", "| DI =",
                  as.numeric(substrlast(names(pars)))[indGC][i])
    ##-------------------------------------------
    fit <- fitsGC[[i]]
    co <- exp(coef(fit$model))
    tco <- rev(unlist(pars[indGC][[i]][-1]))
    niveis <- c(0.9, 0.95, 0.99)
    cortes <- qchisq(niveis, df = 2)
    xy <- levelplot(
        deviance ~ exp(lalpha) + exp(lgamma),
        data = fit$grid, cuts = 30,
        xlab = expression(alpha),
        ylab = expression(gamma),
        main = list(main, cex = 0.8),
        scales = list(y = list(rot = 90)),
        colorkey = list(space = "top"),
        par.settings = ps,
        panel = function(x, y, z, at, region, ...){
            panel.levelplot(x, y, z, at = at, region = TRUE, ...)
            panel.contourplot(x, y, z, ..., at = cortes,
                              contour = TRUE, region = FALSE)
            panel.abline(v = co[1], h = co[2], lty = 2)
            panel.points(x = tco[1], y = tco[2], lty = 2,
                         col = 2, pch = 19)
        })
    return(xy)
})

## Graph
marrangeGrob(plotsGC, ncol = 3, nrow = 1, top = "")

plotsCP <- lapply(seq(fitsCP), function(i) {
    ##-------------------------------------------
    main <- paste("COM-Poisson", "| DI =",
                  as.numeric(substrlast(names(pars)))[indCP][i])
    ##-------------------------------------------
    fit <- fitsCP[[i]]
    co <- exp(coef(fit$model))
    tco <- rev(unlist(pars[indCP][[i]][-1]))
    niveis <- c(0.9, 0.95, 0.99)
    cortes <- qchisq(niveis, df = 2)
    xy <- levelplot(
        deviance ~ exp(lnu) + exp(llambda),
        data = fit$grid, cuts = 30,
        xlab = expression(nu),
        ylab = expression(lambda),
        main = list(main, cex = 0.8),
        scales = list(y = list(rot = 90)),
        colorkey = list(space = "top"),
        par.settings = ps,
        panel = function(x, y, z, at, region, ...){
            panel.levelplot(x, y, z, at = at, region = TRUE, ...)
            panel.contourplot(x, y, z, ..., at = cortes,
                              contour = TRUE, region = FALSE)
            panel.abline(v = co[1], h = co[2], lty = 2)
            panel.points(x = tco[1], y = tco[2], lty = 2,
                         col = 2, pch = 19)
        })
    return(xy)
})

## Graph
marrangeGrob(plotsCP, ncol = 3, nrow = 1, top = "")

## Save
save(fitsCP, fitsGC, indCP, indGC,
     pars, ps, samples, substrlast,
     file = "orthogonality.rda")
