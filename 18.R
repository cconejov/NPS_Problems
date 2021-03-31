## Exercise 1 - Description  ----

## Assess S claims of battery failure from temperatures from a sample of previous batteries that experienced failure and from a sample
## in general 

## a. Perform a kernel density estimation for temps-7 and temps-other using what you consider is the most adequate bandwidth. Since the temperatures are positive,
## is it required to perform any transformation?.

## b. Is there any important difference on the results from considering the LSCV selector over the DPI selector?

## c. It seems that in temps-7 there is a secondary mode. Compute a kernel derivative estimation for temps-7 and temps-other using what you consider are the 
## most adequate bandwidths.

## d. Precisely determine the location of the extreme points

## e. Check with a kernel second derivative that the extreme points are actually modes.
## -------------------------------------------------------------

## 1.1 Preliminary Work ----


## Exercise 2 - Description------------------------------------

## Exercise 2: Compare the MISE and AMISE criteria in three densities in nor1mix of your choice.

## 1. Code (2.33) and the AMISE expression for the normal kernel, and compare the two error curves.

## 2. Compare them for n = 100, 200, 500, adding a vertical line to represent the $h_{MISE}$ and $h_{AMISE}$ bandwidths. Describe in detail the results and the major takeaways.
## -------------------------------------------------------------

## 2.1 Functions ----
## a) Omega_a

omega_a <- function(a, w, h, mu, sigma2){
  k <- length(w)
  M <- matrix(0, ncol = k, nrow = k)
  for(i in 1:k){
    for(j in 1:k){
      M[i,j] = dnorm(x = (mu[i] - mu[j]), 
                     mean = 0, 
                     sd = sqrt(a*h^2 + sigma2[i] + sigma2[j]))
      
    } #End for j
  } # End for i
  return(M)
}

## a) Omega

omega <- function(w, h, n, mu, sigma2){
  
  omega2 <- omega_a(a = 2, w, h, mu, sigma2)
  omega1 <- omega_a(a = 1, w, h, mu, sigma2)
  omega0 <- omega_a(a = 0, w, h, mu, sigma2)
  omega <- (1 - (1/n))*omega2 - 2*omega1 + omega0
  return(t(w) %*% omega %*% w)
}

## c) Explicit MISE 

MISE <- function(h, n, w, mu, sigma2){
  
  len_h   <- length(h)
  MISE_h1 <- numeric(len_h)
  MISE_h  <- numeric(len_h)
  for(i in 1:len_h){MISE_h1[i] <- omega(w, h[i], n, mu, sigma2)}
  MISE_h <- (2*sqrt(pi) *n*h)^(-1) + MISE_h1
  return(MISE_h)
  
}

## d) Explicit AMISE

## Inputs
## n: Number of sample observations
## h: Grid of bandwidth 
## R_K: Squared integral of kernel:
## R_f2: Second derivative of Squared integral of density f

AMISE <- function(n, h, R_K, R_f2){((n*h)^(-1))*R_K + 0.25*(h^4)*R_f2}

## e) Optimal AMISE bandwidth

## Inputs
## n: Number of sample observations.
## R_K: Squared integral of kernel.
## R_f2: Second derivative of Squared integral of density f

h_AMISE_n <- function(n, R_K, R_f2){(R_K/(n*R_f2))^(1/5)}

## f) Hermite Polynomial: Useful for the four derivative of Normal density

Hermite_4 <- function(x, mean, sd){
  
  H4 <- x^4 - 6*x^2 + 3
  Normal4 <-  dnorm(x = x, 
                    mean = mean, 
                    sd = sd)
  
  return(Normal4*H4)
  
}

## g) Comute the Second derivative of Squared integral of density f (Normal mixture)

## Inputs

## w: vector of weights
## mu: Vector of means
## sigma2: Vector of variance

R_f2 <- function(w, mu, sigma2){
  
  k <- length(w)
  sum_int <- 0
  for(i in 1:k){
    for(j in 1:k){
      sum_int <- sum_int + (w[i]*w[j]*Hermite_4(x = (mu[i] - mu[j]), 
                                                mean = 0, 
                                                sd = sqrt(sigma2[i] + sigma2[j])))
    }
  }
  return(sum_int)
}

## 2.2 General inputs ----
# Useful for all scenarios

# Available models
?nor1mix::MarronWand

# Density evaluation
x <- seq(-4, 4, length.out = 400)

# h_grid and Second derivative of the normal kernel
h <- seq(0.04,0.85, length.out = 100)
R_K_N01 <- (2*sqrt(pi))^(-1)


## 2.3 Study Cases ----

# Case 1: nor1mix::MW.nm1 #Gaussian N(0,1) =================================

# Simulating
set.seed(42)
x100 <- nor1mix::rnorMix(n = 100, obj = nor1mix::MW.nm1)
x200 <- nor1mix::rnorMix(n = 200, obj = nor1mix::MW.nm1)
x500 <- nor1mix::rnorMix(n = 500, obj = nor1mix::MW.nm1)


# Histogram plots
par(mfrow = c(1,3))
hist(x100, freq = FALSE, xlim = c(-4,4), main = "Histogram MW.nm1 n = 100")
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm1), col = 2)
hist(x200, freq = FALSE, xlim = c(-4,4), main = "Histogram MW.nm1 n = 200")
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm1), col = 2)
hist(x500, freq = FALSE, xlim = c(-4,4), main = "Histogram MW.nm1 n = 500")
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm1), col = 2)
par(mfrow = c(1,1))


# Second derivative of the density: In this case, it is equivalent to the zero-plug in.
R_f2_f01 <- (3/(8*sqrt(pi)))

## *********************************************
## CASE 1: MISE AND AMISE n = 100
## *********************************************

# AMISE n = 100

AMISE_100 <- AMISE(n = 100, h = h, R_K = R_K_N01, R_f2 = R_f2_f01)
h_AMISE_100 <- h_AMISE_n(n = 100, R_K = R_K_N01, R_f2 = R_f2_f01)


# MISE = 100

MISE_100 <- MISE(h, n = 100, w = 1, mu = 1, sigma2 = 1)
h_MISE_100 <- optimize(MISE, c(0, 1), tol = 0.0001, n = 100, w = 1, mu = 1, sigma2 = 1)

plot(h, MISE_100, type = "l", 
     main = "Comparison MISE(h) and AMISE(h) \n n = 100, nor1mix::MW.nm1",
     ylab = "MISE(h)")
lines(h, AMISE_100, col = "blue")
abline(v = h_AMISE_100, col = "blue")
abline(v = h_MISE_100)
text(0.47, 0.06, labels = paste("h_MISE = ", round(h_MISE_100$minimum,3) ), col = "black", cex = 0.5, font = 2, adj = c(0, 0))
text(0.25, 0.06, labels = paste("h_AMISE = ", round(h_AMISE_100,3) ), col = "blue", cex = 0.5, font = 2, adj = c(0, 0))
legend("topright", legend = c("MISE", "AMISE"), lwd = 1, col = c("black", "blue"))

## *********************************************
## CASE 2: MISE AND AMISE n = 200
## *********************************************

# AMISE n = 200

AMISE_200 <- AMISE(n = 200, h = h, R_K = R_K_N01, R_f2 = R_f2_f01)
h_AMISE_200 <- h_AMISE_n(n = 200, R_K = R_K_N01, R_f2 = R_f2_f01)

# MISE = 200

MISE_200 <- MISE(h, n = 200, w = 1, mu = 1, sigma2 = 1)
h_MISE_200 <- optimize(MISE, c(0, 1), tol = 0.0001, n = 200, w = 1, mu = 1, sigma2 = 1)

plot(h, MISE_200, type = "l", 
     main = "Comparison MISE(h) and AMISE(h) \n n = 200, nor1mix::MW.nm1",
     ylab = "MISE(h)")
lines(h, AMISE_200, col = "blue")
abline(v = h_AMISE_200, col = "blue")
abline(v = h_MISE_200)
text(0.4, 0.03, labels = paste("h_MISE = ", round(h_MISE_200$minimum,3) ), col = "black", cex = 0.5, font = 2, adj = c(0, 0))
text(0.2, 0.03, labels = paste("h_AMISE = ", round(h_AMISE_200,3) ), col = "blue", cex = 0.5, font = 2, adj = c(0, 0))
legend("topright", legend = c("MISE", "AMISE"), lwd = 1, col = c("black", "blue"))


## *********************************************
## CASE 3: MISE AND AMISE n = 500
## *********************************************

# AMISE n = 500

AMISE_500 <- AMISE(n = 500, h = h, R_K = R_K_N01, R_f2 = R_f2_f01)
h_AMISE_500 <- h_AMISE_n(n = 500, R_K = R_K_N01, R_f2 = R_f2_f01)

# MISE = 500

MISE_500 <- MISE(h, n = 500, w = 1, mu = 1, sigma2 = 1)
h_MISE_500 <- optimize(MISE, c(0, 1), tol = 0.0001, n = 500, w = 1, mu = 1, sigma2 = 1)

plot(h, MISE_500, type = "l", 
     main = "Comparison MISE(h) and AMISE(h) \n n = 500, nor1mix::MW.nm1",
     ylab = "MISE(h)")
lines(h, AMISE_500, col = "blue")
abline(v = h_AMISE_500, col = "blue")
abline(v = h_MISE_500)
text(0.4, 0.01, labels = paste("h_MISE = ", round(h_MISE_500$minimum,3) ), col = "black", cex = 0.5, font = 2, adj = c(0, 0))
text(0.1, 0.01, labels = paste("h_AMISE = ", round(h_AMISE_500,3) ), col = "blue", cex = 0.5, font = 2, adj = c(0, 0))
legend("topright", legend = c("MISE", "AMISE"), lwd = 1, col = c("black", "blue"))


# Case 2: nor1mix::MW.nm2 #Skewed =================================

# Simulating
set.seed(42)
x100 <- nor1mix::rnorMix(n = 100, obj = nor1mix::MW.nm2)
x200 <- nor1mix::rnorMix(n = 200, obj = nor1mix::MW.nm2)
x500 <- nor1mix::rnorMix(n = 500, obj = nor1mix::MW.nm2)

# Histogram plots
par(mfrow = c(1,3))
hist(x100, freq = FALSE, xlim = c(-4,4), main = "Histogram MW.nm2 n = 100")
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm2), col = 2)
hist(x200, freq = FALSE, xlim = c(-4,4), main = "Histogram MW.nm2 n = 200")
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm2), col = 2)
hist(x500, freq = FALSE, xlim = c(-4,4), main = "Histogram MW.nm2 n = 500")
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm2), col = 2)
par(mfrow = c(1,1))

# Parameters: Weights, means. sigma and R(f^{2]})

mu_skewed <-  c(-.3, .3, 1)
sigma2_skewed <- c(1.44, .64, 4/9)
w_skewed <- c(.2, .2, .6)
R_f2_skewed <- R_f2(w = w_skewed, mu = mu_skewed, sigma2 = sigma2_skewed)

## *********************************************
## CASE 1: MISE AND AMISE n = 100
## *********************************************

# AMISE n = 100

AMISE_100 <- AMISE(n = 100, h = h, R_K = R_K_N01, R_f2 = R_f2_skewed)
h_AMISE_100 <- h_AMISE_n(n = 100, R_K = R_K_N01, R_f2 = R_f2_skewed)

# MISE = 100

MISE_100 <- MISE(h, n = 100, w = w_skewed, mu = mu_skewed, sigma2 = sigma2_skewed)
h_MISE_100 <- optimize(MISE, c(0, 1), tol = 0.0001, n = 100, w = w_skewed, mu = mu_skewed, sigma2 = sigma2_skewed)

plot(h, MISE_100, type = "l", 
     main = "Comparison MISE(h) and AMISE(h) \n n = 100, nor1mix::MW.nm2",
     ylab = "MISE(h)")
lines(h, AMISE_100, col = "blue")
abline(v = h_AMISE_100, col = "blue")
abline(v = h_MISE_100)
text(0.4, 0.06, labels = paste("h_MISE = ", round(h_MISE_100$minimum,3) ), col = "black", cex = 0.5, font = 2, adj = c(0, 0))
text(0.18, 0.06, labels = paste("h_AMISE = ", round(h_AMISE_100,3) ), col = "blue", cex = 0.5, font = 2, adj = c(0, 0))
legend("topright", legend = c("MISE", "AMISE"), lwd = 1, col = c("black", "blue"))

## *********************************************
## CASE 2: MISE AND AMISE n = 200
## *********************************************

# AMISE n = 200

AMISE_200 <- AMISE(n = 200, h = h, R_K = R_K_N01, R_f2 = R_f2_skewed)
h_AMISE_200 <- h_AMISE_n(n = 200, R_K = R_K_N01, R_f2 = R_f2_skewed)

# MISE = 200

MISE_200 <- MISE(h, n = 200, w = w_skewed, mu = mu_skewed, sigma2 = sigma2_skewed)
h_MISE_200 <- optimize(MISE, c(0, 1), tol = 0.0001, n = 200, w = w_skewed, mu = mu_skewed, sigma2 = sigma2_skewed)

plot(h, MISE_200, type = "l", 
     main = "Comparison MISE(h) and AMISE(h) \n n = 200, nor1mix::MW.nm2",
     ylab = "MISE(h)")
lines(h, AMISE_200, col = "blue")
abline(v = h_AMISE_200, col = "blue")
abline(v = h_MISE_200)
text(0.4, 0.03, labels = paste("h_MISE = ", round(h_MISE_200$minimum,3) ), col = "black", cex = 0.5, font = 2, adj = c(0, 0))
text(0.15, 0.03, labels = paste("h_AMISE = ", round(h_AMISE_200,3) ), col = "blue", cex = 0.5, font = 2, adj = c(0, 0))
legend("topright", legend = c("MISE", "AMISE"), lwd = 1, col = c("black", "blue"))


## *********************************************
## CASE 3: MISE AND AMISE n = 500
## *********************************************

# AMISE n = 500

AMISE_500 <- AMISE(n = 500, h = h, R_K = R_K_N01, R_f2 = R_f2_skewed)
h_AMISE_500 <- h_AMISE_n(n = 500, R_K = R_K_N01, R_f2 = R_f2_skewed)

# MISE = 500

MISE_500 <- MISE(h, n = 500, w = w_skewed, mu = mu_skewed, sigma2 = sigma2_skewed)
h_MISE_500 <- optimize(MISE, c(0, 1), tol = 0.0001, n = 500, w = w_skewed, mu = mu_skewed, sigma2 = sigma2_skewed)

plot(h, MISE_500, type = "l", 
     main = "Comparison MISE(h) and AMISE(h) \n n = 500, nor1mix::MW.nm2",
     ylab = "MISE(h)")
lines(h, AMISE_500, col = "blue")
abline(v = h_AMISE_500, col = "blue")
abline(v = h_MISE_500)
text(0.1, 0.02, labels = paste("h_MISE = ", round(h_MISE_500$minimum,3) ), col = "black", cex = 0.5, font = 2, adj = c(0, 0))
text(0.3, 0.02, labels = paste("h_AMISE = ", round(h_AMISE_500,3) ), col = "blue", cex = 0.5, font = 2, adj = c(0, 0))
legend("topright", legend = c("MISE", "AMISE"), lwd = 1, col = c("black", "blue"))

# Case 3: nor1mix::MW.nm6 # Bimodal =================================

# Simulating
set.seed(42)
x100 <- nor1mix::rnorMix(n = 100, obj = nor1mix::MW.nm6)
x200 <- nor1mix::rnorMix(n = 200, obj = nor1mix::MW.nm6)
x500 <- nor1mix::rnorMix(n = 500, obj = nor1mix::MW.nm6)

# Histogram plots
par(mfrow = c(1,3))
hist(x100, freq = FALSE, xlim = c(-4,4), main = "Histogram MW.nm6 n = 100")
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm6), col = 2)
hist(x200, freq = FALSE, xlim = c(-4,4), main = "Histogram MW.nm6 n = 200")
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm6), col = 2)
hist(x500, freq = FALSE, xlim = c(-4,4), main = "Histogram MW.nm6 n = 500")
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm6), col = 2)
par(mfrow = c(1,1))

# Parameters: weights, means, variance and R(f^{2})

mu_Bimodal <-  c(-1, 1)
sigma2_Bimodal <- c(4/9, 4/9)
w_Bimodal <- c(.5, .5)
R_f2_Bimodal <- R_f2(w = w_Bimodal, mu = mu_Bimodal, sigma2 = sigma2_Bimodal)

## *********************************************
## CASE 1: MISE AND AMISE n = 100
## *********************************************

# AMISE n = 100

AMISE_100 <- AMISE(n = 100, h = h, R_K = R_K_N01, R_f2 = R_f2_Bimodal)
h_AMISE_100 <- h_AMISE_n(n = 100, R_K = R_K_N01, R_f2 = R_f2_Bimodal)

# MISE = 100

MISE_100 <- MISE(h, n = 100, w = w_Bimodal, mu = mu_Bimodal, sigma2 = sigma2_Bimodal)
h_MISE_100 <- optimize(MISE, c(0, 1), tol = 0.0001, n = 100, w = w_Bimodal, mu = mu_Bimodal, sigma2 = sigma2_Bimodal)

plot(h, MISE_100, type = "l", 
     main = "Comparison MISE(h) and AMISE(h) \n n = 100, nor1mix::MW.nm6",
     ylab = "MISE(h)")
lines(h, AMISE_100, col = "blue")
abline(v = h_AMISE_100, col = "blue")
abline(v = h_MISE_100)
text(0.4, 0.06, labels = paste("h_MISE = ", round(h_MISE_100$minimum,3) ), col = "black", cex = 0.5, font = 2, adj = c(0, 0))
text(0.1, 0.06, labels = paste("h_AMISE = ", round(h_AMISE_100,3) ), col = "blue", cex = 0.5, font = 2, adj = c(0, 0))
legend("topright", legend = c("MISE", "AMISE"), lwd = 1, col = c("black", "blue"))

## *********************************************
## CASE 2: MISE AND AMISE n = 200
## *********************************************

# AMISE n = 200

AMISE_200 <- AMISE(n = 200, h = h, R_K = R_K_N01, R_f2 = R_f2_Bimodal)
h_AMISE_200 <- h_AMISE_n(n = 200, R_K = R_K_N01, R_f2 = R_f2_Bimodal)

# MISE = 200

MISE_200 <- MISE(h, n = 200, w = w_Bimodal, mu = mu_Bimodal, sigma2 = sigma2_Bimodal)
h_MISE_200 <- optimize(MISE, c(0, 1), tol = 0.0001, n = 200, w = w_Bimodal, mu = mu_Bimodal, sigma2 = sigma2_Bimodal)

plot(h, MISE_200, type = "l", 
     main = "Comparison MISE(h) and AMISE(h) \n n = 200, nor1mix::MW.nm6",
     ylab = "MISE(h)")
lines(h, AMISE_200, col = "blue")
abline(v = h_AMISE_200, col = "blue")
abline(v = h_MISE_200)
text(0.4, 0.03, labels = paste("h_MISE = ", round(h_MISE_200$minimum,3) ), col = "black", cex = 0.5, font = 2, adj = c(0, 0))
text(0.1, 0.03, labels = paste("h_AMISE = ", round(h_AMISE_200,3) ), col = "blue", cex = 0.5, font = 2, adj = c(0, 0))
legend("topright", legend = c("MISE", "AMISE"), lwd = 1, col = c("black", "blue"))


## *********************************************
## CASE 3: MISE AND AMISE n = 500
## *********************************************

# AMISE n = 500

AMISE_500 <- AMISE(n = 500, h = h, R_K = R_K_N01, R_f2 = R_f2_Bimodal)
h_AMISE_500 <- h_AMISE_n(n = 500, R_K = R_K_N01, R_f2 = R_f2_Bimodal)

# MISE = 500

MISE_500 <- MISE(h, n = 500, w = w_Bimodal, mu = mu_Bimodal, sigma2 = sigma2_Bimodal)
h_MISE_500 <- optimize(MISE, c(0, 1), tol = 0.0001, n = 500, w = w_Bimodal, mu = mu_Bimodal, sigma2 = sigma2_Bimodal)

plot(h, MISE_500, type = "l", 
     main = "Comparison MISE(h) and AMISE(h) \n n = 500, nor1mix::MW.nm2",
     ylab = "MISE(h)")
lines(h, AMISE_500, col = "blue")
abline(v = h_AMISE_500, col = "blue")
abline(v = h_MISE_500)
text(0.3, 0.01, labels = paste("h_MISE = ", round(h_MISE_500$minimum,3) ), col = "black", cex = 0.5, font = 2, adj = c(0, 0))
text(0.05, 0.01, labels = paste("h_AMISE = ", round(h_AMISE_500,3) ), col = "blue", cex = 0.5, font = 2, adj = c(0, 0))
legend("topright", legend = c("MISE", "AMISE"), lwd = 1, col = c("black", "blue"))



## Exercise 3 - Description  ------------------------------------
##Adapt the np_pred_CI function to include the argument type_boot, which can take either the value "naive" or "wild". 

## 1) If type_boot = "wild", then the function must perform the wild bootstrap algorithm described above, implemented from scratch following substeps i–iv. 

## 2) Compare and validate the correct behavior of the confidence intervals, for the two specifications of type_boot, in the model considered in Exercise 5.8 
## (without performing the full simulation study).
## -------------------------------------------------------------

rm(list = ls())

## 3.1 Functions ----

# Function to predict and compute confidence intervals for m(x). 

# 1) Inputs

## 1.1) npfit: A np::npreg object (npfit) 
## 1.2) exdat: Values of the predictors where to carry out prediction (exdat)
## 1.3) B:     Number of bootstrap iterations.
## 1.4) conf:  Range of confidence interval
## 1.5) type_CI: Type of confidence interval. (Based on normal standard or quantiles)
## 1.6) type_boot: Type of bootstrap procedure. Options Naive and Wild bootstrap
## 1.7) perturbed_res: Valid only for Wild Bootstrap. Type of perturbation on the residuals. Options are "normal"or "golden"
## 1.8) seed:

# 2) Outputs

## 2.1) exdat: Values of the predictors where to carry out prediction
## 2.2) m_hat: Predicted regression
## 2.3) lwr:   Lower confidence interval
## 2.4) upr:   Upper confidence interval

np_pred_CI <- function(npfit, 
                       exdat, 
                       B = 200, 
                       conf = 0.95,
                       type_CI = c("standard", "quantiles")[1],
                       type_boot = c("naive", "wild")[1],
                       perturbed_res = c("normal", "golden")[1],
                       seed = 42) {
  
  # Fix seed
  set.seed(seed)
  
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
  
  # Predictions m_hat from the original sample
  m_hat <- np::npreg(txdat = xdat, 
                     tydat = ydat, 
                     exdat = exdat,
                     bws   = npfit$bws)$mean
  
  if (type_boot == "naive") {
    
    # Function for performing naive bootstrap
    boot_function_naive <- function(data, indices) {
      np::npreg(txdat = xdat[indices,], 
                tydat = ydat[indices],
                exdat = exdat, 
                bws   = npfit$bws)$mean
    }
    
    # Carry out the bootstrap estimator
    m_hat_star <- boot::boot(data = data.frame(xdat), 
                             statistic = boot_function_naive,
                             R = B)$t
    
  } else if (type_boot == "wild") {
    
    # Sample size of the predictors
    n <- length(xdat)
    
    # Y fitted
    Y_hat <- npfit$mean
    
    # Ordinary residuals
    residuals_O <- Y_hat - ydat
    
    # Type of perturbation
    if(perturbed_res == "normal"){
      
      # Function for performing wild bootstrap
      boot_function_wild <- function(data, indices) {
        
        # Step i: Simulate V_{i} copies of V (Mean 0 and variance 1)
        V_n <- rnorm(n)
        
        # Step iii. Obtain the bootstrap sample
        ydat_bt <- Y_hat + data[indices]*V_n
        
        np::npreg(txdat = xdat, 
                  tydat = ydat_bt,
                  exdat = exdat, 
                  bws = npfit$bws)$mean
      }
      
      # Step iv. Carry out the wild bootstrap estimator
      m_hat_star <- boot::boot(data = residuals_O, 
                               statistic = boot_function_wild,
                               R = B)$t
      
      
    } else if(perturbed_res == "golden"){
      
      # Function for performing wild bootstrap
      boot_function_wild <- function(data, indices) {
        
        # Step i: Simulate V_{i} copies of V (Mean 0 and variance 1)
        phi <- (1 + sqrt(5))/2
        prob <- (phi + 2)/5  
        
        
        golden <- sample(x = c(1-phi,phi), size = n, prob = c(prob, 1 - prob), replace=T)
        
        # Step iii. Obtain the bootstrap sample
        ydat_bt <- Y_hat + data[indices]*golden
        
        np::npreg(txdat = xdat, 
                  tydat = ydat_bt,
                  exdat = exdat, 
                  bws = npfit$bws)$mean
      }
      
      # Step iv. Carry out the wild bootstrap estimator
      m_hat_star <- boot::boot(data = residuals_O, 
                               statistic = boot_function_wild,
                               R = B)$t
      
    }
    
    else{stop("Incorrect type of peturbation")}
    
  }else{stop("Incorrect type_boot")}
  
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

## 3.2 Simulation ----

set.seed(12345)
B <- 500
n <- 100
eps <- rnorm(n, sd = 0.75)
m <- function(x) 0.25*x^2 - 0.75*x + 3
X <- rnorm(n, sd = 1.5)
Y <- m(X) + eps

plot(X,Y, main = "Simulated observations")


## 3.3 Fitted model npregbw ----

bw1 <- np::npregbw(formula = Y ~ X, regtype = "lc")
fit1 <- np::npreg(bw1)
summary(fit1)

#Asymptotic standard errors
fit1$merr

## 3.4 npplot: Normal approximation ----

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


## 3.5 npplot: Quantile approximation ----

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

## 3.6 Naive Bootstrap ----

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

## 3.7 Wild Bootstrap: Normal Perturbation ----
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
points(X,Y)


## 3.8 Wild Bootstrap: Golden Perturbation ----

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
points(X,Y)