## --------------------------------------------------
## Exercise 2: Compare the MISE and AMISE criteria in three densities in nor1mix of your choice.

## 1. Code (2.33) and the AMISE expression for the normal kernel, and compare the two error curves.

## 2. Compare them for n = 100, 200, 500, adding a vertical line to represent the $h_{MISE}$ and $h_{AMISE}$ bandwidths. Describe in detail the results and the major takeaways.
## --------------------------------------------------

# Available models
?nor1mix::MarronWand

# Simulating
set.seed(42)
samp <- nor1mix::rnorMix(n = 500, obj = nor1mix::MW.nm9) # Trimodal
# MW object in the second argument
hist(samp, freq = FALSE)
# Density evaluation
x <- seq(-4, 4, length.out = 400)
lines(x, nor1mix::dnorMix(x = x, obj = nor1mix::MW.nm9), col = 2)
# Plot a MW object directly
# A normal with the same mean and variance is plotted in dashed lines --
# you can remove it with argument p.norm = FALSE
par(mfrow = c(2, 2))
plot(nor1mix::MW.nm5)
plot(nor1mix::MW.nm7)
plot(nor1mix::MW.nm10)
plot(nor1mix::MW.nm12)
lines(nor1mix::MW.nm7, col = 2) # Also possible
par(mfrow = c(1, 1))


## f1

# Kernel: K ~ N(0,1)
# f density f ~ N(0,1)

plot(nor1mix::MW.nm1)


n <- 100
h <- seq(0.01, 1, by = 0.01 )
w <- c(1)

sigma2 <- 1

MISE1 <- (1/(n*h)) + (1 - (1/n))*((sigma2 + h^2)^(-1/2)) + (1/sqrt(sigma2)) - 2^(3/2)*((2*sigma2 + h^2)^(-1/2))
MISE <- (2*sqrt(pi))^(-1)*MISE1

plot(log(h), MISE1)


# nn



## ----------------------
rm(list = ls())


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


n <- 100
h <- seq(0.01, 1, by = 0.01 )

w <- c(3/4,1/4)
mu <- c(0, 3/2)
sigma2 <- c(1, (1/3)^2) 
MISE1 <- vector(length = length(h))

for(i in 1:length(h)){MISE1[i] <- omega(w, h[i], n, mu, sigma2)}

MISE <- (2*sqrt(pi) *n*h)^(-1) + MISE1

plot(h, MISE)
plot(log(h), MISE,type = "l")












