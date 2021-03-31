## --------------------------------------------------
## Exercise 2: Compare the MISE and AMISE criteria in three densities in nor1mix of your choice.

## 1. Code (2.33) and the AMISE expression for the normal kernel, and compare the two error curves.

## 2. Compare them for n = 100, 200, 500, adding a vertical line to represent the $h_{MISE}$ and $h_{AMISE}$ bandwidths. Describe in detail the results and the major takeaways.
## --------------------------------------------------

# Available models
?nor1mix::MarronWand

# Functions =================================

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


# Inputs: =================================
# Useful for all scenarios

# Density evaluation
x <- seq(-4, 4, length.out = 400)

# h_grid and Second derivative of the normal kernel
h <- seq(0.04,0.85, length.out = 100)
R_K_N01 <- (2*sqrt(pi))^(-1)

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


# Second derivative of the density: In this case, it is equivalent to the zero-plug in.


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


# Second derivative of the density: In this case, it is equivalent to the zero-plug in.


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

