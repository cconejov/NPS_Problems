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

