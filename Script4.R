# Script 4 - Exploring count distributions -----------------------------
# Author: Wagner Hugo Bonat LEG/UFPR -----------------------------------
# Date: 02/03/2017 -----------------------------------------------------

# Loading extra functions ----------------------------------------------
source("Script1.R")
source("Script2.R")
source("Script3.R")

# Painel plot for lattice object used in Figure 2 ----------------------
require(lattice)
require(latticeExtra)
require(RColorBrewer)
ac <- list(pad1=0.5, pad2=0.5, tck=0.5)
mycol <- gray.colors(n=5)
ps <- list(box.rectangle=list(col=1, fill=c("gray70")),
           box.umbrella=list(col=1, lty=1),
           dot.symbol=list(col=1),
           dot.line=list(col=1, lty=3),
           plot.symbol=list(col="gray", cex=0.7, pch = 20),
           plot.line=list(col=1),
           plot.polygon=list(col="gray80"),
           superpose.line=list(col=mycol),
           superpose.symbol=list(col=mycol),
           superpose.polygon=list(col=mycol),
           strip.background=list(col=c("gray90","gray50")),
           layout.widths=list(
               left.padding=0.25,
               right.padding=-1,
               ylab.axis.padding=0),
           layout.heights=list(
               bottom.padding=0.25,
               top.padding=0,
               axis.xlab.padding=0,
               xlab.top=0),
           axis.components=list(bottom=ac, top=ac, left=ac, right=ac)
)

# Parameters
mu = 5:20
x = 40:100
# Gamma-Count ----------------------------------------------------------
# Case 1 - DI = 0.5
gc_mv_0.5 <- moments_gc(beta = mu, alpha = 2.1)
gc_di_0.5 <- disp_index_gc(beta = mu, alpha = 2.1)
gc_zi_0.5 <- zero_inflation_gc(beta = mu, alpha = 2.1)
gc_ht_0.5 <- heavy_tail_gc(x = x, beta = 10, alpha = 2.1)

# Case 2 - DI = 2
gc_mv_2 <- moments_gc(beta = mu, alpha = 0.465)
gc_di_2 <- disp_index_gc(beta = mu, alpha = 0.465)
gc_zi_2 <- zero_inflation_gc(beta = mu, alpha = 0.465)
gc_ht_2 <- heavy_tail_gc(x = x, beta = 10, alpha = 0.465)

# Case 3 - DI = 5
gc_mv_5 <- moments_gc(beta = 2:17, alpha = 0.15)
gc_di_5 <- disp_index_gc(beta = 2:17, alpha = 0.15)
gc_zi_5 <- zero_inflation_gc(beta = 2:17, alpha = 0.15)
gc_ht_5 <- heavy_tail_gc(x = x, beta = 10, alpha = 0.15)

# Case 4 - DI = 20
gc_mv_20 <- moments_gc(beta = seq(0.01, 5, l = 100), alpha = 0.0218)
gc_di_20 <- disp_index_gc(beta = seq(0.01, 5, l = 100), alpha = 0.0218)
gc_zi_20 <- zero_inflation_gc(beta = seq(0.01, 5, l = 100), alpha = 0.0218)
gc_ht_20 <- heavy_tail_gc(x = x, beta = 10, alpha = 0.0218)

# Poisson-Tweedie p = 1.1 ----------------------------------------------
ptw1.1_mv_0.5 <- NA
ptw1.1_di_0.5 <- NA
ptw1.1_zi_0.5 <- NA
ptw1.1_ht_0.5 <- NA

# Case 2 - DI = 2
ptw1.1_mv_2 <- moments_ptw(mu = mu, phi = 0.8, power = 1.1)
ptw1.1_di_2 <- disp_index_ptw(mu = mu, phi = 0.8, power = 1.1)
ptw1.1_zi_2 <- zero_inflation_ptw(mu = mu, phi = 0.8, power = 1.1)
ptw1.1_ht_2 <- heavy_tail_ptw(x = x, mu = 10, phi = 0.8, power = 1.1)

# Case 3 - DI = 5
ptw1.1_mv_5 <- moments_ptw(mu = mu, phi = 3.2, power = 1.1)
ptw1.1_di_5 <- disp_index_ptw(mu = mu, phi = 3.2, power = 1.1)
ptw1.1_zi_5 <- zero_inflation_ptw(mu = mu, phi = 3.2, power = 1.1)
ptw1.1_ht_5 <- heavy_tail_ptw(x = x, mu = 10, phi = 3.2, power = 1.1)

# Case 4 - DI = 5
ptw1.1_mv_20 <- moments_ptw(mu = mu, phi = 15, power = 1.1)
ptw1.1_di_20 <- disp_index_ptw(mu = mu, phi = 15, power = 1.1)
ptw1.1_zi_20 <- zero_inflation_ptw(mu = mu, phi = 15, power = 1.1)
ptw1.1_ht_20 <- heavy_tail_ptw(x = x, mu = 10, phi = 15, power = 1.1)

# Poisson-Tweedie p = 2 ------------------------------------------------
ptw2_mv_0.5 <- NA
ptw2_di_0.5 <- NA
ptw2_zi_0.5 <- NA
ptw2_ht_0.5 <- NA

# Case 2 - DI = 2
ptw2_mv_2 <- moments_ptw(mu = mu, phi = 0.1, power = 2)
ptw2_di_2 <- disp_index_ptw(mu = mu, phi = 0.1, power = 2)
ptw2_zi_2 <- zero_inflation_ptw(mu = mu, phi = 0.1, power = 2)
ptw2_ht_2 <- heavy_tail_ptw(x = x, mu = 10, phi = 0.1, power = 2)

# Case 3 - DI = 5
ptw2_mv_5 <- moments_ptw(mu = mu, phi = 0.4, power = 2)
ptw2_di_5 <- disp_index_ptw(mu = mu, phi = 0.4, power = 2)
ptw2_zi_5 <- zero_inflation_ptw(mu = mu, phi = 0.4, power = 2)
ptw2_ht_5 <- heavy_tail_ptw(x = x, mu = 10, phi = 0.4, power = 2)

# Case 4 - DI = 20
ptw2_mv_20 <- moments_ptw(mu = mu, phi = 1.9, power = 2)
ptw2_di_20 <- disp_index_ptw(mu = mu, phi = 1.9, power = 2)
ptw2_zi_20 <- zero_inflation_ptw(mu = mu, phi = 1.9, power = 2)
ptw2_ht_20 <- heavy_tail_ptw(x = x, mu = 10, phi = 1.9, power = 2)

# Poisson-Tweedie p = 2 ------------------------------------------------
# Case 1 - DI = 0.5
ptw3_mv_0.5 <- NA
ptw3_di_0.5 <- NA
ptw3_zi_0.5 <- NA
ptw3_ht_0.5 <- NA

# Case 2 - DI = 2
ptw3_mv_2 <- moments_ptw(mu = mu, phi = 0.01, power = 3)
ptw3_di_2 <- disp_index_ptw(mu = mu, phi = 0.01, power = 3)
ptw3_zi_2 <- zero_inflation_ptw(mu = mu, phi = 0.01, power = 3)
ptw3_ht_2 <- heavy_tail_ptw(x = x, mu = 10, phi = 0.01, power = 3)

# Case 3 - DI = 5
ptw3_mv_5 <- moments_ptw(mu = mu, phi = 0.04, power = 3)
ptw3_di_5 <- disp_index_ptw(mu = mu, phi = 0.04, power = 3)
ptw3_zi_5 <- zero_inflation_ptw(mu = mu, phi = 0.04, power = 3)
ptw3_ht_5 <- heavy_tail_ptw(x = x, mu = 10, phi = 0.04, power = 3)

# Case 4 - DI = 20
ptw3_mv_20 <- moments_ptw(mu = mu, phi = 0.19, power = 3)
ptw3_di_20 <- disp_index_ptw(mu = mu, phi = 0.19, power = 3)
ptw3_zi_20 <- zero_inflation_ptw(mu = mu, phi = 0.19, power = 3)
ptw3_ht_20 <- heavy_tail_ptw(x = x, mu = 10, phi = 0.19, power = 3)

# COM-Poisson ----------------------------------------------------------
# Case 1 - DI = 0.5
cp_mv_0.5 <- ap_moments_cp(lambda = 20:500, nu = 2.05)
cp_di_0.5 <- disp_index_cp(lambda = 20:500, nu = 2.05)
cp_zi_0.5 <- zero_inflation_cp(lambda = 20:500, nu = 2.05)
cp_ht_0.5 <- heavy_tail_cp(x = 40:100, lambda = 118.51, nu = 2.05)

# Case 2 - DI = 2
cp_mv_2 <- ap_moments_cp(lambda = seq(1.75,4, l = 100), nu = 0.47)
cp_di_2 <- disp_index_cp(lambda = seq(1.75,4, l = 100), nu = 0.47)
cp_zi_2 <- zero_inflation_cp(lambda = seq(1.75,4, l = 100), nu = 0.47)
cp_ht_2 <- heavy_tail_cp(x = 40:100, lambda = 2.88, nu = 0.47)

# Case 3 - DI = 5
cp_mv_5 <- ap_moments_cp(lambda = seq(1, 1.45, l = 100), nu = 0.13)
cp_di_5 <- disp_index_cp(lambda = seq(1, 1.45, l = 100), nu = 0.13)
cp_zi_5 <- zero_inflation_cp(lambda = seq(1, 1.45, l = 100), nu = 0.13)
cp_ht_5 <- heavy_tail_cp(x = 40:100, lambda = 1.30, nu = 0.13)

# Case 3 - DI = 20
cp_mv_20 <- NA
cp_di_20 <- NA
cp_zi_20 <- NA
cp_ht_20 <- NA

# Data set for plotting ------------------------------------------------

# Case 1 - DI = 0.5
# Mean and variance
mv_data0.5 <- data.frame("Mean" = c(gc_mv_0.5[1,], ptw1.1_mv_0.5, ptw2_mv_0.5,
                                 ptw3_mv_0.5, cp_mv_0.5[1,]),
                      "Variance" = c(gc_mv_0.5[2,], ptw1.1_mv_0.5, ptw2_mv_0.5,
                                     ptw3_mv_0.5, cp_mv_0.5[2,]),
                      "Configuration" = "DI = 0.5",
                      "Model" = c(rep("GC", length(gc_mv_0.5[1,])),
                                  rep("PTW1.1", length(ptw1.1_mv_0.5)),
                                  rep("PTW2", length(ptw2_mv_0.5)),
                                  rep("PTW3", length(ptw3_mv_0.5)),
                                  rep("CP", length(cp_mv_0.5[1,]))) )
# Dispersion index
di_data0.5 <- data.frame("Mean" = c(gc_mv_0.5[1,], ptw1.1_mv_0.5, ptw2_mv_0.5,
                                    ptw3_mv_0.5, cp_mv_0.5[1,]),
                         "DI" = c(gc_di_0.5, ptw1.1_di_0.5, ptw2_di_0.5,
                                        ptw3_di_0.5, cp_di_0.5),
                         "Configuration" = "DI = 0.5",
                         "Model" = c(rep("GC", length(gc_mv_0.5[1,])),
                                     rep("PTW1.1", length(ptw1.1_mv_0.5)),
                                     rep("PTW2", length(ptw2_mv_0.5)),
                                     rep("PTW3", length(ptw3_mv_0.5)),
                                     rep("CP", length(cp_mv_0.5[1,]))) )
# Zero-inflation
zi_data0.5 <- data.frame("Mean" = c(gc_mv_0.5[1,], ptw1.1_mv_0.5, ptw2_mv_0.5,
                                    ptw3_mv_0.5, cp_mv_0.5[1,]),
                         "ZI" = c(gc_zi_0.5, ptw1.1_zi_0.5, ptw2_zi_0.5,
                                  ptw3_zi_0.5, cp_zi_0.5),
                         "Configuration" = "DI = 0.5",
                         "Model" = c(rep("GC", length(gc_mv_0.5[1,])),
                                     rep("PTW1.1", length(ptw1.1_mv_0.5)),
                                     rep("PTW2", length(ptw2_mv_0.5)),
                                     rep("PTW3", length(ptw3_mv_0.5)),
                                     rep("CP", length(cp_mv_0.5[1,]))) )
# Heavy tail 
ht_data0.5 <- data.frame("x" = c(40:100, 1, 1,
                                    1, 40:100),
                         "HT" = c(gc_ht_0.5, ptw1.1_ht_0.5, ptw2_ht_0.5,
                                  ptw3_ht_0.5, cp_ht_0.5),
                         "Configuration" = "DI = 0.5",
                         "Model" = c(rep("GC", length(gc_ht_0.5)),
                                     rep("PTW1.1", length(ptw1.1_ht_0.5)),
                                     rep("PTW2", length(ptw2_ht_0.5)),
                                     rep("PTW3", length(ptw3_ht_0.5)),
                                     rep("CP", length(cp_ht_0.5))) )


# Case 2 - DI = 2 ------------------------------------------------------
# Mean and variance
mv_data2 <- data.frame("Mean" = c(gc_mv_2[1,], ptw1.1_mv_2[1,], ptw2_mv_2[1,],
                                    ptw3_mv_2[1,], cp_mv_2[1,]),
                         "Variance" = c(gc_mv_2[2,], ptw1.1_mv_2[2,], ptw2_mv_2[2,],
                                        ptw3_mv_2[2,], cp_mv_2[2,]),
                         "Configuration" = "DI = 2",
                         "Model" = c(rep("GC", length(gc_mv_2[1,])),
                                     rep("PTW1.1", length(ptw1.1_mv_2[1,])),
                                     rep("PTW2", length(ptw2_mv_2[1,])),
                                     rep("PTW3", length(ptw3_mv_2[1,])),
                                     rep("CP", length(cp_mv_2[1,]))) )

# Dispersion index
di_data2 <- data.frame("Mean" = c(gc_mv_2[1,], ptw1.1_mv_2[1,], ptw2_mv_2[1,],
                                    ptw3_mv_2[1,], cp_mv_2[1,]),
                         "DI" = c(gc_di_2, ptw1.1_di_2, ptw2_di_2,
                                  ptw3_di_2, cp_di_2),
                         "Configuration" = "DI = 2",
                         "Model" = c(rep("GC", length(gc_mv_2[1,])),
                                     rep("PTW1.1", length(ptw1.1_mv_2[1,])),
                                     rep("PTW2", length(ptw2_mv_2[1,])),
                                     rep("PTW3", length(ptw3_mv_2[1,])),
                                     rep("CP", length(cp_mv_2[1,]))) )
# Zero-inflation
zi_data2 <- data.frame("Mean" = c(gc_mv_2[1,], ptw1.1_mv_2[1,], ptw2_mv_2[1,],
                                    ptw3_mv_2[1,], cp_mv_2[1,]),
                         "ZI" = c(gc_zi_2, ptw1.1_zi_2, ptw2_zi_2,
                                  ptw3_zi_2, cp_zi_2),
                         "Configuration" = "DI = 2",
                         "Model" = c(rep("GC", length(gc_mv_2[1,])),
                                     rep("PTW1.1", length(ptw1.1_mv_2[1,])),
                                     rep("PTW2", length(ptw2_mv_2[1,])),
                                     rep("PTW3", length(ptw3_mv_2[1,])),
                                     rep("CP", length(cp_mv_2[1,]))) )

# Heavy-tail
ht_data2 <- data.frame("x" = c(40:100, 40:100 , 40:100,
                               40:100, 40:100),
                         "HT" = c(gc_ht_2, ptw1.1_ht_2, ptw2_ht_2,
                                  ptw3_ht_2, cp_ht_2),
                         "Configuration" = "DI = 2",
                         "Model" = c(rep("GC", length(gc_ht_2)),
                                     rep("PTW1.1", length(ptw1.1_ht_2)),
                                     rep("PTW2", length(ptw2_ht_2)),
                                     rep("PTW3", length(ptw3_ht_2)),
                                     rep("CP", length(cp_ht_2))) )

# Case 3 - DI = 5 ------------------------------------------------------
# Mean and variance
mv_data5 <- data.frame("Mean" = c(gc_mv_5[1,], ptw1.1_mv_5[1,], ptw2_mv_5[1,],
                                  ptw3_mv_5[1,], cp_mv_5[1,]),
                       "Variance" = c(gc_mv_5[2,], ptw1.1_mv_5[2,], ptw2_mv_5[2,],
                                      ptw3_mv_5[2,], cp_mv_5[2,]),
                       "Configuration" = "DI = 5",
                       "Model" = c(rep("GC", length(gc_mv_5[1,])),
                                   rep("PTW1.1", length(ptw1.1_mv_5[1,])),
                                   rep("PTW2", length(ptw2_mv_5[1,])),
                                   rep("PTW3", length(ptw3_mv_5[1,])),
                                   rep("CP", length(cp_mv_5[1,]))) )

# Dispersion index
di_data5 <- data.frame("Mean" = c(gc_mv_5[1,], ptw1.1_mv_5[1,], ptw2_mv_5[1,],
                                  ptw3_mv_5[1,], cp_mv_5[1,]),
                       "DI" = c(gc_di_5, ptw1.1_di_5, ptw2_di_5,
                                ptw3_di_5, cp_di_5),
                       "Configuration" = "DI = 5",
                       "Model" = c(rep("GC", length(gc_mv_5[1,])),
                                   rep("PTW1.1", length(ptw1.1_mv_5[1,])),
                                   rep("PTW2", length(ptw2_mv_5[1,])),
                                   rep("PTW3", length(ptw3_mv_5[1,])),
                                   rep("CP", length(cp_mv_5[1,]))) )

# Zero-inflation
zi_data5 <- data.frame("Mean" = c(gc_mv_5[1,], ptw1.1_mv_5[1,], ptw2_mv_5[1,],
                                  ptw3_mv_5[1,], cp_mv_5[1,]),
                       "ZI" = c(gc_zi_5, ptw1.1_zi_5, ptw2_zi_5,
                                ptw3_zi_5, cp_zi_5),
                       "Configuration" = "DI = 5",
                       "Model" = c(rep("GC", length(gc_mv_5[1,])),
                                   rep("PTW1.1", length(ptw1.1_mv_5[1,])),
                                   rep("PTW2", length(ptw2_mv_5[1,])),
                                   rep("PTW3", length(ptw3_mv_5[1,])),
                                   rep("CP", length(cp_mv_5[1,]))) )

# Heavy tail
ht_data5 <- data.frame("x" = c(40:100, 40:100 , 40:100,
                               40:100, 40:100),
                       "HT" = c(gc_ht_5, ptw1.1_ht_5, ptw2_ht_5,
                                ptw3_ht_5, cp_ht_5),
                       "Configuration" = "DI = 5",
                       "Model" = c(rep("GC", length(gc_ht_5)),
                                   rep("PTW1.1", length(ptw1.1_ht_5)),
                                   rep("PTW2", length(ptw2_ht_5)),
                                   rep("PTW3", length(ptw3_ht_5)),
                                   rep("CP", length(cp_ht_5))) )



# Case 4 - DI = 20 -----------------------------------------------------
# Mean and variance
mv_data20 <- data.frame("Mean" = c(gc_mv_20[1,], ptw1.1_mv_20[1,], ptw2_mv_20[1,],
                                  ptw3_mv_20[1,], cp_mv_20),
                       "Variance" = c(gc_mv_20[2,], ptw1.1_mv_20[2,], ptw2_mv_20[2,],
                                      ptw3_mv_20[2,], cp_mv_20),
                       "Configuration" = "DI = 20",
                       "Model" = c(rep("GC", length(gc_mv_20[1,])),
                                   rep("PTW1.1", length(ptw1.1_mv_20[1,])),
                                   rep("PTW2", length(ptw2_mv_20[1,])),
                                   rep("PTW3", length(ptw3_mv_20[1,])),
                                   rep("CP", length(cp_mv_20))) )
# Dispersion index
di_data20 <- data.frame("Mean" = c(gc_mv_20[1,], ptw1.1_mv_20[1,], ptw2_mv_20[1,],
                                  ptw3_mv_20[1,], cp_mv_20),
                       "DI" = c(gc_di_20, ptw1.1_di_20, ptw2_di_20,
                                ptw3_di_20, cp_di_20),
                       "Configuration" = "DI = 20",
                       "Model" = c(rep("GC", length(gc_mv_20[1,])),
                                   rep("PTW1.1", length(ptw1.1_mv_20[1,])),
                                   rep("PTW2", length(ptw2_mv_20[1,])),
                                   rep("PTW3", length(ptw3_mv_20[1,])),
                                   rep("CP", length(cp_mv_20))) )
# Zero-inflation
zi_data20 <- data.frame("Mean" = c(gc_mv_20[1,], ptw1.1_mv_20[1,], ptw2_mv_20[1,],
                                  ptw3_mv_20[1,], cp_mv_20),
                       "ZI" = c(gc_zi_20, ptw1.1_zi_20, ptw2_zi_20,
                                ptw3_zi_20, cp_zi_20),
                       "Configuration" = "DI = 20",
                       "Model" = c(rep("GC", length(gc_mv_20[1,])),
                                   rep("PTW1.1", length(ptw1.1_mv_20[1,])),
                                   rep("PTW2", length(ptw2_mv_20[1,])),
                                   rep("PTW3", length(ptw3_mv_20[1,])),
                                   rep("CP", length(cp_mv_20))) )

# Heavy-tail
ht_data20 <- data.frame("x" = c(40:100, 40:100 , 40:100,
                               40:100, 1),
                       "HT" = c(gc_ht_20, ptw1.1_ht_20, ptw2_ht_20,
                                ptw3_ht_20, cp_ht_20),
                       "Configuration" = "DI = 20",
                       "Model" = c(rep("GC", length(gc_ht_20)),
                                   rep("PTW1.1", length(ptw1.1_ht_20)),
                                   rep("PTW2", length(ptw2_ht_20)),
                                   rep("PTW3", length(ptw3_ht_20)),
                                   rep("CP", length(cp_ht_20))) )
mv_data <- rbind(mv_data0.5, mv_data2, mv_data5, mv_data20)
di_data <- rbind(di_data0.5, di_data2, di_data5, di_data20)
zi_data <- rbind(zi_data0.5, zi_data2, zi_data5, zi_data20)
ht_data <- rbind(ht_data0.5, ht_data2, ht_data5, ht_data20)
mv_data$Type <- "MV"
di_data$Type <- "DI"
zi_data$Type <- "ZI"
names(mv_data)[2] <- "Index"
names(di_data)[2] <- "Index"
names(zi_data)[2] <- "Index"
aux_data <- rbind(mv_data, di_data, zi_data)
# Plot -----------------------------------------------------------------
levels(aux_data$Configuration) <- c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
levels(mv_data$Configuration) <- c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
levels(di_data$Configuration) <- c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
levels(zi_data$Configuration) <- c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
levels(ht_data$Configuration) <- c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
key <- list(text=list(levels(di_data$Model)),
            title="Model", cex.title=1.1,
            lines=list(pch=c(1,2,15,16,18)), 
            col = c("black","gray10","gray30","gray50", "gray70"),
            lty = 1:5,
            divide=c(1), columns=5, type="l")

# Mean and variance relationship ---------------------------------------
aux_data$Type <- factor(aux_data$Type, levels = c("ZI","DI","MV"))
head(aux_data)
AUX.plot <- xyplot(Index ~ Mean | Configuration + Type, ylab = "Index",
                  xlab = expression(mu),
                  data = aux_data, group = Model, type="l",
                  par.settings = ps, scales = list(y = list(relation = "free")),
                  horizontal = TRUE, 
                  pch = c(1,2,15,16,18), 
                  lty = 1:5,
                  layout = c(4,3),
                  col = c("black","gray10","gray30","gray50","gray70"),
                  cex = 0.5, key = key)
AUX.plot <- useOuterStrips(AUX.plot)

HT.plot <- xyplot(HT ~ x | Configuration, ylab = "HT",
                  xlab = "y",
                  data = ht_data, group = Model, type="l",
                  par.settings = ps, scales = "free",
                  horizontal = TRUE, 
                  pch = c(1,2,15,16,18), 
                  lty = 1:5,
                  layout = c(4,1),
                  col = c("black","gray10","gray30","gray50","gray70"),
                  cex = 0.5, xlim = c(39,100))
# END ------------------------------------------------------------------
