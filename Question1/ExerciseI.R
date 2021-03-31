# Read data
rm(list = ls())

temps_7 <- read.table(file = "temps-7.txt", header = TRUE)
temps_other <- read.table(file = "temps-other.txt", header = TRUE)

temps_7 <- temps_7$x
temps_other <- temps_other$x

sum


##a. Perform a kernel density estimation for temps-7 and temps-other using what you consider is
## the most adequate bandwidth. Since the temperatures are positive, is it required to perform any
## transformation?

summary(temps_7)
summary(temps_other)

par(mfrow = c(1,2))
hist(temps_7, main = 'Problematic battery')
hist(temps_other, main = 'Past battery')
par(mfrow = c(1,1))

plot(density(x = temps_7))

par(mfrow = c(1,2))
plot(density(x = temps_7), xlim = c(-5,65), main = 'Problematic battery')
plot(density(x = temps_other), xlim = c(-5,65), main = 'Past battery')
par(mfrow = c(1,1))


# Deal with the problem of Select a particular bandwidth
par(mfrow = c(1,2))
plot(density(x = temps_7))
plot(density(x = temps_7, bw = 0.5))
par(mfrow = c(1,1))


# Transformed the data: No necessary. The min observation is away from zero! (No assigment of probability)

## kde with log-transformed data
kde <- density(log(temps_7))
plot(kde, main = "Kde of transformed data")
range(kde$x)

# Untransform kde$x so the grid is in (0, infty)
kde_transf <- kde
kde_transf$x <- exp(kde_transf$x)
# Transform the density using the chain rule
kde_transf$y <- kde_transf$y * 1 / kde_transf$x

par(mfrow = c(1,2))
plot(density(x = temps_7), xlim = c(-5,65), main = 'Problematic battery')
# Transformed kde
plot(kde_transf, main = "Problematic battery: Transformed kde", xlim = c(-5, 65))
par(mfrow = c(1,1))


## SAME WITH temps_other

## kde with log-transformed data
kde <- density(log(temps_other))
plot(kde, main = "Kde of transformed data")
range(kde$x)

# Untransform kde$x so the grid is in (0, infty)
kde_transf <- kde
kde_transf$x <- exp(kde_transf$x)

# Transform the density using the chain rule
kde_transf$y <- kde_transf$y * 1 / kde_transf$x

par(mfrow = c(1,2))
plot(density(x = temps_other), xlim = c(-5,60), main = 'Past battery')

# Transformed kde
plot(kde_transf, main = "Past battery: Transformed kde", xlim = c(-5, 60))
par(mfrow = c(1,1))


## b. Is there any important difference on the results from considering the LSCV selector over the DPI selector?

summary(temps_7)
par(mfrow = c(1,2))
plot(density(x = temps_7))
plot(density(x = temps_7), xlim = c(-5,65))
par(mfrow = c(1,1))

bw.nrd(x = temps_7)  # Rule Thumb rue
bw.nrd0(x = temps_7) # Default option!!!


bw.ucv(x = temps_7)

# Eduardo function

source("functions/bw_ucv_mod.R")
bw_ucv <- bw.ucv.mod(x = temps_7$x, plot_cv = TRUE, h_grid = 10^seq(-1.5, 0.5, l = 200))



bw_dpi <- bw.SJ(x = temps_7, method = "dpi")
ks::hpi(temps_7)
  
plot(density(x = temps_7), xlim = c(-5,65))
plot(density(x = temps_7, bw = bw_ucv), xlim = c(-5,65))
plot(density(x = temps_7, bw = bw_dpi), xlim = c(-5,65))


##c. It seems that in temps-7 there is a secondary mode. Compute a kernel derivative estimation for
## temps-7 and temps-other using what you consider are the most adequate bandwidths.

# OFE

plot(density(x = temps_7), xlim = c(-5,65))
pos_Max <- which.max(density(x = temps_7)$y)
pos_Max
density(x = temps_7)$x[pos_Max]

plot(density(x = temps_7), xlim = c(-5,65))
abline(v = density(x = temps_7)$x[pos_Max] )


## Find the minimun (Then, from the minumn, look up the max in the rest of the dist)

# MIN
YY <- density(x = temps_7)$y[density(x = temps_7)$x > 20 & density(x = temps_7)$x < 40]
minY <- min(YY)
posMin <- which(density(x = temps_7)$y == minY)
density(x = temps_7)$x[posMin]
## [1] 33.83911

# MAX
YY <- density(x = temps_7)$y[density(x = temps_7)$x > 33 & density(x = temps_7)$x < 60]
maxY <- max(YY)
posMax <- which(density(x = temps_7)$y == maxY)
density(x = temps_7)$x[posMax]



## CREATE A FUNCTION:

densMode <- function(x){
  td <- density(x)
  maxDens <- which.max(td$y)
  list(x=td$x[maxDens], y=td$y[maxDens])
}
densMode(temps_7)

plot(density(x = temps_7), xlim = c(-5,65), main = "temp_7: Density Estimation")
abline(v = density(x = temps_7)$x[pos_Max], col = "3")
abline(v = density(x = temps_7)$x[posMax], col = "3")


kdde_0 <- ks::kdde(x = temps_7, deriv.order = 0)
plot(kdde_0, xlab = "x", main = "Density estimation")
abline(v = density(x = temps_7)$x[pos_Max], col = "3")
abline(v = density(x = temps_7)$x[posMax], col = "3")

# Density derivative estimation (automatically chosen bandwidth, but different
# from kdde_0!)

kdde_1 <- ks::kdde(x = temps_7, deriv.order = 1)
plot(kdde_1, xlab = "x", main = "Density derivative estimation")
abline(v = density(x = temps_7)$x[pos_Max], col = "3")
abline(v = density(x = temps_7)$x[posMax], col = "3")
abline(h = 0)

kdde_2 <- ks::kdde(x = temps_7, deriv.order = 2)
plot(kdde_2, xlab = "x", main = "Density second derivative estimation")
abline(v = density(x = temps_7)$x[pos_Max], col = "3")
abline(v = density(x = temps_7)$x[posMax], col = "3")
abline(h = 0)


# Found the max and minima
plot(kde$x, kde$y)
density(x = temps_7)$y
x <- temps_7
dens <- function(x){density(x = x)$y}
nlm(f = dens, p = 15)


minus_dens <- function(x){-dens(x)}
extrema <- c(nlm(f = minus_dens, p = 15)$estimate,
             nlm(f = minus_dens, p = 40)$estimate)

extrema

## d. Precisely determine the location of the extreme points.
##e. Check with a kernel second derivative that the extreme points are actually modes.
