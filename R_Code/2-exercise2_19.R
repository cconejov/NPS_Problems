## --------------------------------------------------
## Exercise 2: Compare the MISE and AMISE criteria in three densities in nor1mix of your choice.

## 1. Code (2.33) and the AMISE expression for the normal kernel, and compare the two error curves.

## 2. Compare them for n = 100, 200, 500, adding a vertical line to represent the $h_{MISE}$ and $h_{AMISE}$ bandwidths. Describe in detail the results and the major takeaways.
## --------------------------------------------------

# Available models
?nor1mix::MarronWand

# Functions =================================

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

omega <- function(w, h, n, mu, sigma2){
  
  omega2 <- omega_a(a = 2, w, h, mu, sigma2)
  omega1 <- omega_a(a = 1, w, h, mu, sigma2)
  omega0 <- omega_a(a = 0, w, h, mu, sigma2)
  omega <- (1 - (1/n))*omega2 - 2*omega1 + omega0
  return(t(w) %*% omega %*% w)
}

MISE <- function(h, n, w, mu, sigma2){
  
  len_h   <- length(h)
  MISE_h1 <- numeric(len_h)
  MISE_h  <- numeric(len_h)
  for(i in 1:len_h){MISE_h1[i] <- omega(w, h[i], n, mu, sigma2)}
  MISE_h <- (2*sqrt(pi) *n*h)^(-1) + MISE_h1
  return(MISE_h)
  
}

# Case 1: nor1mix::MW.nm1 N(0,1) =================================

# Simulating
set.seed(42)
x100 <- nor1mix::rnorMix(n = 100, obj = nor1mix::MW.nm1)
x200 <- nor1mix::rnorMix(n = 200, obj = nor1mix::MW.nm1)
x500 <- nor1mix::rnorMix(n = 500, obj = nor1mix::MW.nm1)

# Density evaluation
x <- seq(-4, 4, length.out = 400)

# Histogram plots
par(mfrow = c(1,3))
hist(x100, freq = FALSE, xlim = c(-4,4), main = "Histogram MW.nm1 n = 100")
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm1), col = 2)
hist(x200, freq = FALSE, xlim = c(-4,4), main = "Histogram MW.nm1 n = 200")
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm1), col = 2)
hist(x500, freq = FALSE, xlim = c(-4,4), main = "Histogram MW.nm1 n = 500")
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm1), col = 2)
par(mfrow = c(1,1))


# Computing MISE AND AMISE

h <- seq(0.04,0.85, length.out = 100)
R_f2 <- (3/(8*sqrt(pi)))
R_K <- (2*sqrt(pi))^(-1)

AMISE <- function(n, h){((n*h)^(-1))*R_K + 0.25*(h^4)*R_f2}
h_AMISE_n <- function(n){(R_K/(n*R_f2))^(1/5)}


## *********************************************
## CASE 1: MISE AND AMISE n = 100
## *********************************************

# AMISE n = 100

h_AMISE_100 <- h_AMISE_n(n = 100)
AMISE_100 <- AMISE(n = 100, h = h)

# MISE = 100

MISE_100 <- MISE(h, n = 100, w = 1, mu = 1, sigma2 = 1)
h_MISE_100 <- optimize(MISE, c(0, 1), tol = 0.0001, n = 100, w = 1, mu = 1, sigma2 = 1)

plot(h, MISE_100, type = "l", 
     main = "Comparison MISE(h) and AMISE(h) \n n = 100, nor1mix::MW.nm1",
     ylab = "MISE(h)")
lines(h, AMISE_100, col = "3")
abline(v = h_AMISE_100, col = "3")
abline(v = h_MISE_100)
legend("topright", legend = c("MISE", "AMISE"), lwd = 1, col = c(1,3))


## *********************************************
## CASE 2: MISE AND AMISE n = 200
## *********************************************

# AMISE n = 200

h_AMISE_200 <- h_AMISE_n(n = 200)
AMISE_200 <- AMISE(n = 200, h = h)

# MISE = 200

MISE_200 <- MISE(h, n = 200, w = 1, mu = 1, sigma2 = 1)
h_MISE_200 <- optimize(MISE, c(0, 1), tol = 0.0001, n = 200, w = 1, mu = 1, sigma2 = 1)

plot(h, MISE_200, type = "l", 
     main = "Comparison MISE(h) and AMISE(h) \n n = 200, nor1mix::MW.nm1",
     ylab = "MISE(h)")
lines(h, AMISE_200, col = "3")
abline(v = h_AMISE_200, col = "3")
abline(v = h_MISE_200)
legend("topright", legend = c("MISE", "AMISE"), lwd = 1, col = c(1,3))

## *********************************************
## CASE 3: MISE AND AMISE n = 500
## *********************************************

# AMISE n = 500

h_AMISE_500 <- h_AMISE_n(n = 500)
AMISE_500 <- AMISE(n = 500, h = h)

# MISE = 500

MISE_500 <- MISE(h, n = 500, w = 1, mu = 1, sigma2 = 1)
h_MISE_500 <- optimize(MISE, c(0, 1), tol = 0.0001, n = 500, w = 1, mu = 1, sigma2 = 1)

plot(h, MISE_500, type = "l", 
     main = "Comparison MISE(h) and AMISE(h) \n n = 500, nor1mix::MW.nm1",
     ylab = "MISE(h)")
lines(h, AMISE_500, col = "3")
abline(v = h_AMISE_500, col = "3")
abline(v = h_MISE_500)
legend("topright", legend = c("MISE", "AMISE"), lwd = 1, col = c(1,3))