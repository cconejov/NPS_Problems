## --------------------------------------------------
##Adapt the np_pred_CI function to include the argument type_boot, which can take either the value "naive" or "wild". 

## 1) If type_boot = "wild", then the function must perform the wild bootstrap algorithm described above, implemented from scratch following substeps i–iv. 

## 2) Compare and validate the correct behavior of the confidence intervals, for the two specifications of type_boot, in the model considered in Exercise 5.8 
## (without performing the full simulation study).
## --------------------------------------------------

rm(list =ls())

set.seed(12345)
B <- 500
n <- 100
eps <- rnorm(n, sd = 0.75)
m <- function(x) 0.25*x^2 - 0.75*x + 3
X <- rnorm(n, sd = 1.5)
Y <- m(X) + eps


plot(X,Y, main = "Simulated observations")

bw1 <- np::npregbw(formula = Y ~ X, regtype = "lc")
fit1 <- np::npreg(bw1)
summary(fit1)

#Asymptotic standard errors
fit1$merr

# Normal approximation confidence intervals + extraction of errors
npplot_std <- plot(fit1, 
                   plot.errors.method = "bootstrap",
                   plot.errors.type = "standard", 
                   plot.errors.boot.num = B,
                   plot.errors.style = "bar", 
                   plot.behavior = "plot-data",
                   lwd = 2,
                   main = "Bootstrap (npplot) \n Normal Approximation Confidence Intervals")
lines(npplot_std$r1$eval[, 1], npplot_std$r1$mean + npplot_std$r1$merr[, 1],
      col = 2, lty = 2)
lines(npplot_std$r1$eval[, 1], npplot_std$r1$mean + npplot_std$r1$merr[, 2],
      col = 2, lty = 2)
# These bootstrap standard errors are different from the asymptotic ones
head(npplot_std$r1$merr)


# Quantile confidence intervals + extraction of errors
npplot_qua <- plot(fit1, 
                   plot.errors.method = "bootstrap",
                   plot.errors.type = "quantiles", 
                   plot.errors.boot.num = B,
                   plot.errors.style = "bar", 
                   plot.behavior = "plot-data",
                   main = "Bootstrap (npplot) \n Quantile Confidence Intervals")
lines(npplot_qua$r1$eval[, 1], npplot_qua$r1$mean + npplot_qua$r1$merr[, 1],
      col = 2, lty = 2)
lines(npplot_qua$r1$eval[, 1], npplot_qua$r1$mean + npplot_qua$r1$merr[, 2],
      col = 2, lty = 2)

# There is no predict() method featuring bootstrap confidence intervals,
# it has to be coded manually!
# Function to predict and compute confidence intervals for m(x). Takes as main
# arguments a np::npreg object (npfit) and the values of the predictors where
# to carry out prediction (exdat). Requires that exdat is a data.frame of the
# same type than the one used for the predictors (e.g., there will be an error
# if one variable appears as a factor for computing npfit and then is passed
# as a numeric in exdat)
source("R_Code/np_pred_CI.R") 


# Obtain predictions and confidence intervals along a fine grid, using the
# same seed employed by np::npplot for proper comparison

ci1 <- np_pred_CI(npfit = fit1, 
                  exdat = seq(-5, 5, by = 0.1),
                  B = B, 
                  type_CI = "quantiles",
                  type_boot = "naive")
# Reconstruction of np::npplot’s figure -- the curves coincide perfectly
plot(fit1, 
     plot.errors.method = "bootstrap", 
     plot.errors.type = "quantiles",
     plot.errors.boot.num = B, 
     plot.errors.style = "bar", 
     lwd = 3,  
     main = "Confidence Intervals quantile \n Naive Bootstrap",
     xlim = c(-5,5))
lines(npplot_qua$r1$eval[, 1], npplot_qua$r1$mean + npplot_qua$r1$merr[, 1],
      col = 2, lwd = 3)
lines(npplot_qua$r1$eval[, 1], npplot_qua$r1$mean + npplot_qua$r1$merr[, 2],
      col = 2, lwd = 3)
lines(ci1$exdat, ci1$m_hat, col = 3)
lines(ci1$exdat, ci1$lwr, col = 4)
lines(ci1$exdat, ci1$upr, col = 4)
points(X,Y)

Y
   
ci1$lwr


## Wild bootstrap (Normal perturbation)


ci1 <- np_pred_CI(npfit = fit1, 
                  exdat = seq(-5, 5, by = 0.1),
                  B = B, 
                  type_CI = "quantiles",
                  type_boot = "wild")
# Reconstruction of np::npplot’s figure -- the curves does not coincide perfectly: Broad intervals
plot(fit1, 
     plot.errors.method = "bootstrap", 
     plot.errors.type = "quantiles",
     plot.errors.boot.num = B, 
     plot.errors.style = "bar", 
     lwd = 3,
     main = "Confidence Intervals quantile \n Wild Bootstrap - Normal perturbation",
     xlim = c(-5,5) )
lines(npplot_qua$r1$eval[, 1], npplot_qua$r1$mean + npplot_qua$r1$merr[, 1],
      col = 2, lwd = 3)
lines(npplot_qua$r1$eval[, 1], npplot_qua$r1$mean + npplot_qua$r1$merr[, 2],
      col = 2, lwd = 3)
lines(ci1$exdat, ci1$m_hat, col = 3)
lines(ci1$exdat, ci1$lwr, col = 4)
lines(ci1$exdat, ci1$upr, col = 4)

## Wild bootstrap (Golden ratio perturbation)
ci1 <- np_pred_CI(npfit = fit1, 
                  exdat = seq(-5, 5, by = 0.1),
                  B = B, 
                  type_CI = "quantiles",
                  type_boot = "wild",
                  perturbed_res = "golden")
# Reconstruction of np::npplot’s figure -- the curves does not coincide perfectly: Broad intervals
plot(fit1, 
     plot.errors.method = "bootstrap", 
     plot.errors.type = "quantiles",
     plot.errors.boot.num = B, 
     plot.errors.style = "bar", 
     lwd = 3,
     main = "Confidence Intervals quantile \n Wild Bootstrap - Golden ratio perturbation",
     xlim = c(-5,5) )
lines(npplot_qua$r1$eval[, 1], npplot_qua$r1$mean + npplot_qua$r1$merr[, 1],
      col = 2, lwd = 3)
lines(npplot_qua$r1$eval[, 1], npplot_qua$r1$mean + npplot_qua$r1$merr[, 2],
      col = 2, lwd = 3)
lines(ci1$exdat, ci1$m_hat, col = 3)
lines(ci1$exdat, ci1$lwr, col = 4)
lines(ci1$exdat, ci1$upr, col = 4)
