## Exercise 1 - Description  ----

## Exercise 1: Assess S claims of battery failure from temperatures from a sample of previous batteries that experienced failure and from a sample
## in general 

## a. Perform a kernel density estimation for temps-7 and temps-other using what you consider is the most adequate bandwidth. Since the temperatures are positive,
## is it required to perform any transformation?.

## b. Is there any important difference on the results from considering the LSCV selector over the DPI selector?

## c. It seems that in temps-7 there is a secondary mode. Compute a kernel derivative estimation for temps-7 and temps-other using what you consider are the 
## most adequate bandwidths.

## d. Precisely determine the location of the extreme points

## e. Check with a kernel second derivative that the extreme points are actually modes.

## Preliminary Work ----

# Load data

problemPhones <- read.table("temps-7.txt",     header = TRUE)$x
pastPhones    <- read.table("temps-other.txt", header = TRUE)$x

# Summaries of data

summary(problemPhones)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.85   12.51   15.68   17.34   19.52   62.84 

summary(pastPhones)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.39   12.69   17.22   17.76   22.01   46.86

# Histograms - default bandwidths

par(mfrow=c(1,2))
p1 <- hist(problemPhones, xlim=c(0,65), probability = T, col=rgb(1,0,0,1/4), main = "Problematic Phones")
p2 <- hist(pastPhones, xlim=c(0,65), probability = T, col=rgb(0,0,1,1/4), main = "Past Phones")
par(mfrow=c(1,1))


# Density plots - default bandwidths

par(mfrow = c(1,2))
plot(density(x = problemPhones), xlim = c(-5,65), main = 'Problematic battery')
plot(density(x = pastPhones), xlim = c(-5,65), main = 'Past battery')
par(mfrow = c(1,1))


## Part a. ----
## Perform a kernel density estimation for temps-7 and temps-other using what you consider is the most adequate bandwidth. 
## Since the temperatures are positive, is it required to perform any transformation?.

# Default bandwidth bw.nrd0:

bw.nrd0(x = problemPhones)
bw.nrd0(x = pastPhones)


# TRANSFORMATIONS

# !Problematic temperatures =================================

## kde with log-transformed data
kde_problemPhones <- density(log(problemPhones))
plot(kde_problemPhones, main = "Problematic battery: transformed data")
range(kde_problemPhones$x)

# Untransform kde$x so the grid is in (0, infty)
kde_transf_problemPhones <- kde_problemPhones
kde_transf_problemPhones$x <- exp(kde_transf_problemPhones$x)
# Transform the density using the chain rule
kde_transf_problemPhones$y <- kde_transf_problemPhones$y * 1 / kde_transf_problemPhones$x

par(mfrow = c(1,2))
plot(density(x = problemPhones), xlim = c(-5,65), main = 'Problematic battery')
# Transformed kde
plot(kde_transf_problemPhones, main = "Problematic battery: Transformed", xlim = c(-5, 65))
par(mfrow = c(1,1))

# !Past temperatures =================================

## kde with log-transformed data
kde_pastPhones <- density(log(pastPhones))
plot(kde_pastPhones, main = "Kde of transformed data")
range(kde_pastPhones$x)

# Untransform kde$x so the grid is in (0, infty)
kde_transf_pastPhones <- kde_pastPhones
kde_transf_pastPhones$x <- exp(kde_transf_pastPhones$x)

# Transform the density using the chain rule
kde_transf_pastPhones$y <- kde_transf_pastPhones$y * 1 / kde_transf_pastPhones$x

par(mfrow = c(1,2))
plot(density(x = pastPhones), xlim = c(-5,60), main = 'Past battery')
# Transformed kde
plot(kde_transf_pastPhones, main = "Past battery: Transformed", xlim = c(-5, 60))
par(mfrow = c(1,1))


# Transformed bandwidths:
kde_transf_problemPhones$bw
kde_transf_pastPhones$bw



# BANDWIDTH SELECTION

# The four methods are as below:

# RT
# DPI
# LSCV
# BCV


# RT


#Past temperatures
bw.nrd(x = pastPhones) 
#  1.151129
#similar to:
lengthpastPhones <- length(pastPhones)
iqr <- diff(quantile(pastPhones, c(0.25, 0.75))) / diff(qnorm(c(0.25, 0.75))) 
1.06 * lengthpastPhones^(-1/5) * min(sd(pastPhones), iqr)
#  1.151129


#Problematic temperatures
bw.nrd(x = problemPhones) 
#  0.873937
#similar to:
lengthproblemPhones <- length(problemPhones) 
iqr <- diff(quantile(problemPhones, c(0.25, 0.75))) / diff(qnorm(c(0.25, 0.75))) 
1.06 * lengthproblemPhones^(-1/5) * min(sd(problemPhones), iqr)
#  0.8797933



# DPI


#Past temperatures
bw.SJ(x = pastPhones, method = "dpi")
#  0.841367
# Similar but faster:
ks::hpi(x = pastPhones) 
#  0.8413892


#Problematic temperatures
bw.SJ(x = problemPhones, method = "dpi")
#  0.6303672
# Similar but faster:
ks::hpi(x = problemPhones) 
#  0.6307256


# LSCV 


# Past temperatures 
# UCV R-default function
bw.ucv(x = pastPhones) #does not provide a warning
# 0.6420771
# UCV R-default function with extended grid
bw.ucv(x = pastPhones, lower = 0.01, upper = 1) #also is fine
#  0.64334

#Class function

bw.ucv.mod <- function(x, nb = 1000L,
                       h_grid = 10^seq(-3, log10(1.2 * sd(x) *
                                                   length(x)^(-1/5)), l = 200),
                       plot_cv = FALSE) {
  if ((n <- length(x)) < 2L)
    stop("need at least 2 data points")
  n <- as.integer(n)
  if (is.na(n))
    stop("invalid length(x)")
  if (!is.numeric(x))
    stop("invalid 'x'")
  nb <- as.integer(nb)
  if (is.na(nb) || nb <= 0L)
    stop("invalid 'nb'")
  storage.mode(x) <- "double"
  hmax <- 1.144 * sqrt(var(x)) * n^(-1/5)
  Z <- .Call(stats:::C_bw_den, nb, x)
  d <- Z[[1L]]
  cnt <- Z[[2L]]
  fucv <- function(h) .Call(stats:::C_bw_ucv, n, d, cnt, h)
  ## Original
  # h <- optimize(fucv, c(lower, upper), tol = tol)$minimum
  # if (h < lower + tol | h > upper - tol)
  #   warning("minimum occurred at one end of the range")
  ## Modification
  obj <- sapply(h_grid, function(h) fucv(h))
  h <- h_grid[which.min(obj)]
  if (h %in% range(h_grid)) 
    warning("minimum occurred at one end of h_grid")
  if (plot_cv) {
    plot(h_grid, obj, type = "o")
    rug(h_grid)
    abline(v = h, col = 2, lwd = 2)
  }
  h
}

bw.ucv.mod(x = pastPhones, h_grid = 10^seq(-1.5, 0.5, l = 200))





# Problem temperatures 
# UCV R-default function
bw.ucv(x = problemPhones) #does not provide a warning
# 0.6524331
# UCV R-default function with extended grid
bw.ucv(x = problemPhones, lower = 0.01, upper = 1) #also is fine
# 0.6514426
# Class function
bw.ucv.mod(x = problemPhones, h_grid = 10^seq(-1.5, 0.5, l = 200))
# minimum occurred at end of the grid - value: 0.001


# BCV


# Past temperatures
# BCV R-default function
bw.bcv(x = pastPhones)
## 0.8934792
# BCV R-default function extended search interval
bw.bcv(x = x, lower = 0.01, upper = 1)
## 0.8922003

# Class function

bw.bcv.mod <- function(x, nb = 1000L,
                       h.grid = 10^seq(-3, log10(1.2 * sd(x) *
                                                   length(x)^(-1/5)), l = 200),
                       plot.cv = FALSE) {
  if ((n <- length(x)) < 2L)
    stop("need at least 2 data points")
  n <- as.integer(n)
  if (is.na(n))
    stop("invalid length(x)")
  if (!is.numeric(x))
    stop("invalid 'x'")
  nb <- as.integer(nb)
  if (is.na(nb) || nb <= 0L)
    stop("invalid 'nb'")
  storage.mode(x) <- "double"
  hmax <- 1.144 * sqrt(var(x)) * n^(-1/5)
  Z <- .Call(stats:::C_bw_den, nb, x)
  d <- Z[[1L]]
  cnt <- Z[[2L]]
  fbcv <- function(h) .Call(stats:::C_bw_bcv, n, d, cnt, h)
  # h <- optimize(fbcv, c(lower, upper), tol = tol)$minimum
  # if (h < lower + tol | h > upper - tol)
  #   warning("minimum occurred at one end of the range")
  obj <- sapply(h.grid, function(h) fbcv(h))
  h <- h.grid[which.min(obj)]
  if (h %in% range(h.grid)) 
    warning("minimum occurred at one end of h.grid")
  if (plot.cv) {
    plot(h.grid, obj, type = "o")
    rug(h.grid)
    abline(v = h, col = 2, lwd = 2)
  }
  h
}
bw.bcv.mod(x = pastPhones  , nb = 1000L)
# 0.8766243


# Problematic temperatures
# BCV R-default function
bw.bcv(x = problemPhones)
## 0.6142516
# BCV R-default function extened search interval
bw.bcv(x = problemPhones, lower = 0.01, upper = 1)
##  0.6138301
# Class function
bw.bcv.mod(x = problemPhones  , nb = 1000L)
# 0.6149622

# OPTIMAL BANDWIDTH PLOTS


# Past temperature bandwidths:
# - RT: 1.15
# - DPI: 0.84
# - LSCV: 0.64
# - BCV: 0.88


# Problematic temperature bandwidths:
# - RT: 0.87 
# - DPI: 0.63
# - LSCV: 0.66
# - BCV:  0.61




## Part b. ----
## Is there any important difference on the results from considering the LSCV selector over the DPI selector?

# !Problematic temperatures =================================

#LSCV

bw.lscv.problemPhones <- bw.ucv(x = problemPhones)

#DPI

bw_dpi_problemPhones <- bw.SJ(x = problemPhones, method = "dpi")

# !Past temperatures =================================

#LSCV

bw.lscv.pastPhones <- bw.ucv(x = pastPhones)

#DPI

bw_dpi_pastPhones <- bw.SJ(x = pastPhones, method = "dpi")


#PLOTS

# !Problematic temperatures =================================

par(mfrow=c(1,2))
plot(density(x = problemPhones, bw = bw.lscv.problemPhones), xlim = c(-5,65), main = "Density LSCV: Problematic")
plot(density(x = problemPhones, bw = bw_dpi_problemPhones), xlim = c(-5,65), main = "Density DPI: Problematic")
par(mfrow=c(1,1))


# !Past temperatures =================================

par(mfrow=c(1,2))
plot(density(x = pastPhones, bw = bw.lscv.pastPhones), xlim = c(-5,65), main = "Density LSCV: Past") 
plot(density(x = pastPhones, bw = bw_dpi_pastPhones), xlim = c(-5,65), main = "Density DPI: Past")
par(mfrow=c(1,1))



## Part c. ----
## It seems that in temps-7 there is a secondary mode. Compute a kernel derivative estimation for temps-7 and temps-other using what you consider are the most adequate bandwidths.

# !Problematic temperatures =================================

## Detect Global mode
posMax1 <- which.max(density(x = problemPhones)$y)

## Detect minima
YY <- density(x = problemPhones)$y[density(x = problemPhones)$x > 20 & density(x = problemPhones)$x < 40]
minY <- min(YY)
posMin <- which(density(x = problemPhones)$y == minY)

## Detect Secondary mode
YY <- density(x = problemPhones)$y[density(x = problemPhones)$x > 33 & density(x = problemPhones)$x < 60]
maxY <- max(YY)
posMax2 <- which(density(x = problemPhones)$y == maxY)


kdde_0_problemPhones <- ks::kdde(x = problemPhones, deriv.order = 0)
kdde_1_problemPhones <- ks::kdde(x = problemPhones, deriv.order = 1)

par(mfrow = c(1,2))
plot(kdde_0_problemPhones, xlab = "x", main = "Problematic battery: Density estimation", xlim = c(-5,65))
abline(v = density(x = problemPhones)$x[posMax1], col = "3")
abline(v = density(x = problemPhones)$x[posMax2], col = "3")
abline(v = density(x = problemPhones)$x[posMin], col = "2")
plot(kdde_1_problemPhones, xlab = "x", main = "Problematic battery: Density derivative estimation")
abline(v = density(x = problemPhones)$x[posMax1], col = "3")
abline(v = density(x = problemPhones)$x[posMax2], col = "3")
abline(v = density(x = problemPhones)$x[posMin], col = "2")
abline(h = 0)
par(mfrow = c(1,1))

# !Past temperatures =================================

#plot(density(x = pastPhones), xlim = c(-5,65), main = 'Past battery')

## Detect Global mode
pos_Max_PastPhone <- which.max(density(x = pastPhones)$y)


kdde_0_pastPhones <- ks::kdde(x = pastPhones, deriv.order = 0)
kdde_1_pastPhones <- ks::kdde(x = pastPhones, deriv.order = 1)

par(mfrow = c(1,2))
plot(kdde_0_pastPhones, xlab = "x", main = "Past battery: Density estimation", xlim = c(-5,65))
abline(v = density(x = pastPhones)$x[pos_Max_PastPhone], col = "3")
plot(kdde_1_pastPhones, xlab = "x", main = "Past battery: Density derivative estimation")
abline(v = density(x = pastPhones)$x[pos_Max_PastPhone], col = "3")
abline(h = 0)
par(mfrow = c(1,1))


## Part d. ---- 
## Precisely determine the location of the extreme points.

# !Problematic temperatures =================================

density(x = problemPhones)$x[posMax1]  
density(x = problemPhones)$x[posMax2]
density(x = problemPhones)$x[posMin]

density(x = problemPhones)$y[posMax1]
density(x = problemPhones)$y[posMax2] 
density(x = problemPhones)$y[posMin]


# !Problematic temperatures =================================

density(x = pastPhones)$x[posMax2]
density(x = pastPhones)$y[posMax2]  


## Part e. ---- 
## Check with a kernel second derivative that the extreme points are actually modes.

# !Problematic temperatures =================================

kdde_2_problemPhones <- ks::kdde(x = problemPhones, deriv.order = 2)
par(mfrow = c(1,3))
plot(kdde_0_problemPhones, xlab = "x", main = "Problematic series: \n Density estimation", xlim = c(-5,65))
abline(v = density(x = problemPhones)$x[posMax1], col = "3")
abline(v = density(x = problemPhones)$x[posMax2], col = "3")
abline(v = density(x = problemPhones)$x[posMin], col = "2")
plot(kdde_1_problemPhones, xlab = "x", main = "Problematic series: \n Density derivative estimation")
abline(v = density(x = problemPhones)$x[posMax1], col = "3")
abline(v = density(x = problemPhones)$x[posMax2], col = "3")
abline(v = density(x = problemPhones)$x[posMin], col = "2")
abline(h = 0)
plot(kdde_2_problemPhones, xlab = "x", main = "Problematic series: \n  Density second derivative estimation")
abline(v = density(x = problemPhones)$x[posMax1], col = "3")
abline(v = density(x = problemPhones)$x[posMax2], col = "3")
abline(v = density(x = problemPhones)$x[posMin], col = "2")
abline(h = 0)
par(mfrow = c(1,1))



# !Past temperatures =================================

kdde_2_pastPhones <- ks::kdde(x = pastPhones, deriv.order = 2)

par(mfrow = c(1,3))
plot(kdde_0_pastPhones, xlab = "x", main = "Past series \n Density estimation", xlim = c(-5,65))
abline(v = density(x = pastPhones)$x[pos_Max_PastPhone], col = "3")
plot(kdde_1_pastPhones, xlab = "x", main = "Past battery \n Density derivative estimation")
abline(v = density(x = pastPhones)$x[pos_Max_PastPhone], col = "3")
abline(h = 0)
plot(kdde_2_pastPhones, xlab = "x", main = "Past series \n Density Second derivative estimation")
abline(v = density(x = pastPhones)$x[pos_Max_PastPhone], col = "3")
abline(h = 0)
par(mfrow = c(1,1))


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