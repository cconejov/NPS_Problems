## --------------------------------------------------
##Adapt the np_pred_CI function to include the argument type_boot, which can take either the value "naive" or "wild". 

## 1) If type_boot = "wild", then the function must perform the wild bootstrap algorithm described above, implemented from scratch following substeps i–iv. 

## 2) Compare and validate the correct behavior of the confidence intervals, for the two specifications of type_boot, in the model considered in Exercise 5.8 
## (without performing the full simulation study).
## --------------------------------------------------

library(np)


## Generate some data
set.seed(12345)
n <- 100
eps <- rnorm(n, sd = 0.75)
m <- function(x) 0.25*x^2 - 0.75*x + 3
X <- rnorm(n, sd = 1.5)
Y <- m(X) + eps
x_grid <- seq(-5,5, by = 0.1)

round(seq(-5,5, length.out = 100),1)

## Local constant fit

bw0 <- np::npregbw(formula = Y ~ X, nmulti = 2)
kre0 <- np::npreg(bws = bw0)

plot(kre0, col = 2, type = "o")
points(X,Y)
rug(X, side = 1); rug(Y, side = 2)

## Local linear fit

bw1 <- np::npregbw(formula = Y ~ X, regtype = "ll", nmulti = 2)
kre1 <- np::npreg(bws = bw1)

plot(kre1, col = 2, type = "o")
points(X,Y)
rug(X, side = 1); rug(Y, side = 2)

## Both fits

x_grid <- seq(-5,5, by = 0.1)
kre0 <- np::npreg(bws = bw0, exdat = x_grid)
#kre1 <- np::npreg(bws = bw1, exdat = x_grid) # GOOD FOR VIS, BUT DOES NOT WORK with CI_np_pred
kre1 <- np::npreg(bws = bw1)

plot(X,Y)
lines(kre0$eval$x_grid, kre0$mean, col = 2)
lines(kre1$eval$x_grid, kre0$mean, col = 3)
rug(X, side = 1); rug(Y, side = 2)
legend("topright", legend = c("Nadaraya-Watson", "Local linear"),
       lwd = 2, col = 2:3)


# Asymptotic confidence bands for the marginal effects of each predictor on the response
plot(kre1, plot.errors.method = "asymptotic", common.scale = FALSE,
     plot.par.mfrow = FALSE, main = "CI Asympotic")

# Bootstrap confidence bands (using naive bootstrap, the default)
# They take more time to compute because a resampling + refitting takes place
B <- 500
plot(kre1, plot.errors.method = "bootstrap", common.scale = FALSE,
     plot.par.mfrow = FALSE, plot.errors.boot.num = B, random.seed = 42, main = "CI bootstrap")


# Extract the info
# Asymptotic (not bootstrap) standard errors
head(kre1$merr)

## [1] 7.585337 5.755449 4.446989 3.497504 2.798633 2.277110

# Normal approximation confidence intervals + extraction of errors
npplot_std <- plot(kre1, 
                   plot.errors.method = "bootstrap",
                   plot.errors.type = "standard", 
                   plot.errors.boot.num = B,
                   plot.errors.style = "bar", 
                   plot.behavior = "plot-data",
                   lwd = 2,
                   main = "Normal Confidence Interval")
lines(npplot_std$r1$eval[, 1], 
      npplot_std$r1$mean + npplot_std$r1$merr[, 1],
      col = 2, lty = 2)
lines(npplot_std$r1$eval[, 1], npplot_std$r1$mean + npplot_std$r1$merr[, 2],
      col = 2, lty = 2)

# These bootstrap standard errors are different from the asymptotic ones
head(npplot_std$r1$merr)

# Quantile confidence intervals + extraction of errors
npplot_qua <- plot(kre1, 
                   plot.errors.method = "bootstrap",
                   plot.errors.type = "quantiles", 
                   plot.errors.boot.num = B,
                   plot.errors.style = "bar", 
                   plot.behavior = "plot-data",
                   main = "Quantile Confidence Interval")
lines(npplot_qua$r1$eval[, 1], npplot_qua$r1$mean + npplot_qua$r1$merr[, 1],
      col = 2, lty = 2)
lines(npplot_qua$r1$eval[, 1], npplot_qua$r1$mean + npplot_qua$r1$merr[, 2],
      col = 2, lty = 2)
# These bootstrap standard errors are different from the asymptotic ones,
# and also from the previous bootstrap errors (different confidence
# interval method)
head(npplot_qua$r1$merr)

# There is no predict() method featuring bootstrap confidence intervals,
# it has to be coded manually!

# Function to predict and compute confidence intervals for m(x). Takes as main
# arguments a np::npreg object (npfit) and the values of the predictors where
# to carry out prediction (exdat). Requires that exdat is a data.frame of the
# same type than the one used for the predictors (e.g., there will be an error
# if one variable appears as a factor for computing npfit and then is passed
# as a numeric in exdat)
np_pred_CI <- function(npfit, 
                       exdat, 
                       B = 200, 
                       conf = 0.95,
                       type_CI = c("standard", "quantiles")[1]) {
  # Extract predictors
  xdat <- npfit$eval
  # Extract response, using a trick from np::npplot.rbandwidth
  
  tt <- terms(npfit$bws)
  tmf <- npfit$bws$call[c(1, match(c("formula", "data"),
                                   names(npfit$bws$call)))]
  tmf[[1]] <- as.name("model.frame")
  tmf[["formula"]] <- tt
  tmf <- eval(tmf, envir = environment(tt))
  ydat <- model.response(tmf)
  # Predictions
  m_hat <- np::npreg(txdat = xdat, tydat = ydat, exdat = exdat,
                     bws = npfit$bws)$mean
  # Function for performing step 3
  boot_function <- function(data, indices) {
    np::npreg(txdat = xdat[indices,], tydat = ydat[indices],
              exdat = exdat, bws = npfit$bws)$mean
  }
  # Carry out step 3
  m_hat_star <- boot::boot(data = data.frame(xdat), statistic = boot_function,
                           R = B)$t
  # Confidence intervals
  alpha <- 1 - conf
  if (type_CI == "standard") {
    z <- qnorm(p = 1 - alpha / 2)
    se <- apply(m_hat_star, 2, sd)
    lwr <- m_hat - z * se
    upr <- m_hat + z * se
  } else if (type_CI == "quantiles") {
    q <- apply(m_hat_star, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2))
    lwr <- q[1, ]
    upr <- q[2, ]
  } else {
    stop("Incorrect type_CI")
  }
  # Return evaluation points, estimates, and confidence intervals
  return(data.frame("exdat" = exdat, "m_hat" = m_hat, "lwr" = lwr, "upr" = upr))
}


B =50
# Obtain predictions and confidence intervals along a fine grid, using the
# same seed employed by np::npplot for proper comparison
set.seed(12345)
ci1 <- np_pred_CI(npfit = kre1, 
                  B = B, 
                  exdat = x_grid,
                  type_CI = "quantiles")

# Reconstruction of np::npplot’s figure -- the curves coincide perfectly
plot(fit1, plot.errors.method = "bootstrap", plot.errors.type = "quantiles",
     plot.errors.boot.num = B, plot.errors.style = "bar", lwd = 3)
lines(npplot_qua$r1$eval[, 1], npplot_qua$r1$mean + npplot_qua$r1$merr[, 1],
      col = 2, lwd = 3)
lines(npplot_qua$r1$eval[, 1], npplot_qua$r1$mean + npplot_qua$r1$merr[, 2],
      col = 2, lwd = 3)
lines(ci1$exdat, ci1$m_hat, col = 3)
lines(ci1$exdat, ci1$lwr, col = 4)
lines(ci1$exdat, ci1$upr, col = 4)
