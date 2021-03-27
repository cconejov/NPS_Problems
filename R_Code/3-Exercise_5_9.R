## --------------------------------------------------
##Adapt the np_pred_CI function to include the argument type_boot, which can take either the value "naive" or "wild". 

## 1) If type_boot = "wild", then the function must perform the wild bootstrap algorithm described above, implemented from scratch following substeps i–iv. 

## 2) Compare and validate the correct behavior of the confidence intervals, for the two specifications of type_boot, in the model considered in Exercise 5.8 
## (without performing the full simulation study).
## --------------------------------------------------

rm(list =ls())

set.seed(12345)
B <- 200
n <- 100
eps <- rnorm(n, sd = 0.75)
m <- function(x) 0.25*x^2 - 0.75*x + 3
X <- rnorm(n, sd = 1.5)
Y <- m(X) + eps

bw1 <- np::npregbw(formula = Y ~ X, regtype = "lc")
fit1 <- np::npreg(bw1)
summary(fit1)

#Asymptotic standard errors
fit1$merr


# Normal approximation confidence intervals + extraction of errors
npplot_std <- plot(fit1, plot.errors.method = "bootstrap",
                   plot.errors.type = "standard", plot.errors.boot.num = B,
                   plot.errors.style = "bar", plot.behavior = "plot-data",
                   lwd = 2)
lines(npplot_std$r1$eval[, 1], npplot_std$r1$mean + npplot_std$r1$merr[, 1],
      col = 2, lty = 2)
lines(npplot_std$r1$eval[, 1], npplot_std$r1$mean + npplot_std$r1$merr[, 2],
      col = 2, lty = 2)
# These bootstrap standard errors are different from the asymptotic ones
head(npplot_std$r1$merr)


# Quantile confidence intervals + extraction of errors
npplot_qua <- plot(fit1, plot.errors.method = "bootstrap",
                   plot.errors.type = "quantiles", plot.errors.boot.num = B,
                   plot.errors.style = "bar", plot.behavior = "plot-data")
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
np_pred_CI <- function(npfit, exdat, B = 200, conf = 0.95,
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
# Obtain predictions and confidence intervals along a fine grid, using the
# same seed employed by np::npplot for proper comparison
set.seed(42)
ci1 <- np_pred_CI(npfit = fit1, B = B, exdat = seq(-5, 5, by = 0.1),
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

