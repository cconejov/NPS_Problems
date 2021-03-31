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

problemphones <- read.table("temps-7.txt", header = TRUE)
pastphones <- read.delim("temps-other.txt", header = TRUE)


# Rename variables

names(problemphones)[names(problemphones) == "x"] <- "problemtemps"
names(pastphones)[names(pastphones) == "x"] <- "pasttemps"

# Attaching data for ease of use

attach(problemphones)
attach(pastphones)

# Summaries of data

summary(problemtemps)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.85   12.51   15.68   17.34   19.52   62.84 

summary(pastphones)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.39   12.69   17.22   17.76   22.01   46.86


# Histograms - default bandwidths

p1 <- hist(problemtemps)
p2 <- hist(pasttemps)

plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,65))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,65), add=T)  

# Density plots - default bandwidths

par(mfrow = c(1,2))
plot(density(x = temps_7), xlim = c(-5,65), main = 'Problematic battery')
plot(density(x = temps_other), xlim = c(-5,65), main = 'Past battery')
par(mfrow = c(1,1))

## what is it using is the default bandwidth? I believe that it's the DPI.


## Part a. ----
## Perform a kernel density estimation for temps-7 and temps-other using what you consider is the most adequate bandwidth. Since the temperatures are positive,
## is it required to perform any transformation?.


##    BANDWIDTH SELECTION


# We will use several techniques to select the optimal bandwidth we discussed over the term. These process are processes are the "rule of thumb" bandwidth selection 
# derived from the combining from minimizing the AMISE with a normal parametric assumption for the unknown $R(f'')$, part of the AMISE (asympotic mean integrated
# squared error)  derivation: $h_{AMISE} = [\frac{R(K)}{\mu^2_2(K)R(f'')n}]^\frac{1}{5}$. We then look at the Direct Plug In (DPI) selector, which instead of using a normal 
# parameterization in the zero stage, the DPI buries the parametric assumption, usually, two stages between bias (lower with number of stages) and variance
#(higher with the number of stages). We then look at two cross-validation techniques, where we attempt to minimize the MISE (mean integrated squared error),
# through means of cross validation, or leave-one-out techniques. There are two types: the Least Squared Cross-Validation (or Unbiased Cross Validation) and the
# Biased Cross Validation. Numerical optimization is required for the LSCV and therefore we can get trapped in spurious solutions, necessitating the limiting
# of the bandwidth grid for the $h_{LSCV}$. The BCV adapts a hybrid strategy estimating a modification of $R(f'')$ with leave-out-diagonals. It is important that
# the plot or the grid is reviewed in this selector as well as we are looking for the local minimizer not the global as h, the bandwidth length, goes to infinity.
# The $\hat{h}_{LSCV}$ has more much more variance than the $\hat{h}_{BCV}$, this reduction in variance tends to come at the cost of the increased bias, making
# $\hat{h}_{BCV}$ usually larger than $\hat{h}_{LSCV}$. 

# The four methods are as below:

# RT
# DPI
# LSCV
# BCV



# RT


#Past temperatures

bw.nrd(x = pasttemps) 

# [1] 1.151129

#similar to:

lengthpasttemps <- length(pasttemps)

iqr <- diff(quantile(pasttemps, c(0.25, 0.75))) / diff(qnorm(c(0.25, 0.75))) 
1.06 * lengthpasttemps^(-1/5) * min(sd(pasttemps), iqr)

# [1] 1.151129


#Problematic temperatures

bw.nrd(x = problemtemps) 

# [1] 0.873937

#similar to:

lengthproblemtemps <- length(problemtemps)

iqr <- diff(quantile(problemtemps, c(0.25, 0.75))) / diff(qnorm(c(0.25, 0.75))) 
1.06 * lengthproblemtemps^(-1/5) * min(sd(problemtemps), iqr)

# [1] 0.8797933



# DPI


#Past temperatures

bw.SJ(x = pasttemps, method = "dpi")

# [1] 0.841367

# Similar but faster:

ks::hpi(x = pasttemps) 

# [1] 0.8413892


#Problematic temperatures

bw.SJ(x = problemtemps, method = "dpi")

# [1] 0.6303672

# Similar but faster:

ks::hpi(x = problemtemps) 

# [1] 0.6307256



# LSCV 


# Past temperatures 

# UCV R-default function

bw.ucv(x = pasttemps) #does not provide a warning

# [1] 0.6420771

# UCV R-default function with extended grid

bw.ucv(x = pasttemps, lower = 0.01, upper = 1) #also is fine

# [1] 0.64334

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

bw.ucv.mod(x = pasttemps, nb = 5000)

# minimum occurred at end of the grid - value: 0.001



# Problem temperatures 

# UCV R-default function

bw.ucv(x = problemtemps) #does not provide a warning

# [1] 0.6524331

# UCV R-default function with extended grid

bw.ucv(x = problemtemps, lower = 0.01, upper = 1) #also is fine

# [1] 0.6514426

# Class function

bw.ucv.mod(x = pasttemps, nb = 5000)

# minimum occurred at end of the grid - value: 0.001




# BCV


# Past temperatures

# BCV R-default function

bw.bcv(x = pasttemps)

## [1] 0.8934792

# BCV R-default function extended search interval

bw.bcv(x = x, lower = 0.01, upper = 1)

## [1] 0.8922003

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

bw.bcv.mod(x = pasttemps  , nb = 1000L)

# [1] 0.8766243


# Problematic temperatures


# BCV R-default function

bw.bcv(x = problemtemps)

## [1] 0.6142516

# BCV R-default function extened search interval

bw.bcv(x = problemtemps, lower = 0.01, upper = 1)

## [1] 0.6138301

# Class function

bw.bcv.mod(x = problemtemps  , nb = 1000L)

# [1] 0.6149622


# Reviewing our results, we see that that the RT gives us bandwidths that are much quite large in comparison with the other bandwidth selectors
#, which is common for non-normal data like our own. From the literature, we know from  the DPI selector has a convergence rate that is much faster 
# than the cross validation selectors, LSCV being one of them. However, we do have quite large sample sizes. Therefore, with less
# sample size it is likely to perform better than the LSCV. Cross-validtory selectors can be much better suited to highly non-normal or rough selectors, where a DPI 
# selectors can over smooth.  We can see that comparing the plots and the bandwidths below..


# TRANSFORMATIONS


# Transformations must be performed in order to avoid boundary bias, which is to assign probability, with the function mapping, when there is no real density to support it.
# Hence, we are looking at temperatures that theoretically can take values below zero or close to zero to  a very  high values that can be considered
# postive infinity. In reality, in our sample, does not have any values close to zero, therefore boundary bias is a risk. From this, we can deduct that a log
# transfomration could be possible but we will check if this is necesaary in our sample. 

# Past temperatures

## kde with log-transformed data
kde <- density(log(pasttemps))
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


# Problematic temperatures

## kde with log-transformed data
kde <- density(log(problemtemps))
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


# SHOULD WE TAKE THE PROBIT?

# We see that transformation is not necessary as our estimation does not place probability where there is no density. A transformed curve almost replicates that 
# of a non-transformed curve.


## Part b. ----
## Is there any important difference on the results from considering the LSCV selector over the DPI selector?

#LSCV

#Past temperatures

bw.lscv.pasttemps <- bw.ucv(x = pasttemps)
                                              
#Problem Temperatures

bw.lscv.problemtemps <- bw.ucv(x = problemtemps)


#DPI

#Past temperatures

bw_dpi_pasttemps <- bw.SJ(x = pasttemps, method = "dpi")
bw_dpi_pasttemps

#Problem Temperatures

bw_dpi_problemtemps <- bw.SJ(x = problemtemps, method = "dpi")
bw_dpi_problemtemps


#PLOTS

# Problem Temperatures

par(mfrow=c(1,2))
plot(density(x = problemtemps, bw = bw.lscv.problemtemps), xlim = c(-5,65)) 
plot(density(x = problemtemps, bw = bw_dpi_problemtemps), xlim = c(-5,65))


# Past Temperatures

par(mfrow=c(1,2))
plot(density(x = pasttemps, bw = bw.lscv.pasttemps), xlim = c(-5,65)) #this is incorrect.
plot(density(x = pasttemps, bw = bw_dpi_pasttemps), xlim = c(-5,65))




# We can see that the bandwidths for the past temperatures differ significantly more than those for the problematic temperatures. The difference between the bandwidths
# is less than 0.02 for the problematic temperatures. It is almost 0.2 for the past temperatures. In the literature, the DPI has much higher convergence rate
# than cross-validation methods, like LSCV, however, with very volatile or non-normal data the DPI may over smooth. We may be seeing this in practice
# as the LSCV does have a smaller bandwidth that is capturing the disruption in the past temperatures series. The sample size is also considering larger in the 
# problematic temperature series so this may allow the two approaches to more closely approach one another.


## Part c. ----
## It seems that in temps-7 there is a secondary mode. Compute a kernel derivative estimation for
## temps-7 and temps-other using what you consider are the most adequate bandwidths.


# CHANGE THE BANDWIDTHS

# We can create a plot for the two modes of the problematic series below and see the two modes, which are identified by the two lines in green in the plot.
# We then take the first derivative and can see that the derivaive estimator crosses the the x-axis at thoes times - indicating it is truly a peak in the series.

densMode <- function(x){
  td <- density(x)
  maxDens <- which.max(td$y)
  list(x=td$x[maxDens], y=td$y[maxDens])
}
densMode(problemtemps)

par(mfrow=c(1,2))

kdde_0 <- ks::kdde(x = problemtemps, deriv.order = 0)
plot(kdde_0, xlab = "x", main = "Density estimation")
abline(v = density(x = problemtemps)$x[pos_Max], col = "3")
abline(v = density(x = problemtemps)$x[posMax], col = "3")


kdde_1 <- ks::kdde(x = problemtemps, deriv.order = 1)
plot(kdde_1, xlab = "x", main = "Density derivative estimation")
abline(v = density(x = problemtemps)$x[pos_Max], col = "3")
abline(v = density(x = problemtemps)$x[posMax], col = "3")
abline(h = 0)


# We can then review the results of the problematic series, where because of the distortion towards, the


densMode <- function(x){
  td <- density(x)
  maxDens <- which.max(td$y)
  list(x=td$x[maxDens], y=td$y[maxDens])
}
densMode(problemtemps)

kdde_0 <- ks::kdde(x = problemtemps, deriv.order = 0)
plot(kdde_0, xlab = "x", main = "Density estimation")
abline(v = density(x = problemtemps)$x[pos_Max], col = "3")
abline(v = density(x = problemtemps)$x[posMax], col = "3")


kdde_1 <- ks::kdde(x = problemtemps, deriv.order = 1)
plot(kdde_1, xlab = "x", main = "Density derivative estimation")
abline(v = density(x = problemtemps)$x[pos_Max], col = "3")
abline(v = density(x = problemtemps)$x[posMax], col = "3")
abline(h = 0)



# We can then see that in the past telephone series that the modes are more difficult to identify as the first derivative crosses the axis many times going to zero.
# This reflects what we touched upon earlier that there is some distortion in this series particularly towards the high point of the series.


densMode <- function(x){
  td <- density(x)
  maxDens <- which.max(td$y)
  list(x=td$x[maxDens], y=td$y[maxDens])
}
densMode(pasttemps)

par(mfrow=c(1,2))

kdde_0 <- ks::kdde(x = pasttemps, deriv.order = 0)
plot(kdde_0, xlab = "x", main = "Density estimation")
abline(v = density(x = pasttemps)$x[pos_Max], col = "3")
abline(v = density(x = pasttemps)$x[posMax], col = "3")


kdde_1 <- ks::kdde(x = pasttemps, deriv.order = 1)
plot(kdde_1, xlab = "x", main = "Density derivative estimation")
abline(v = density(x = pasttemps)$x[pos_Max], col = "3")
abline(v = density(x = pasttemps)$x[posMax], col = "3")
abline(h = 0)


## Part d. ---- 
## Precisely determine the location of the extreme points.

kde <


plot(kde$x, kde$y)
density(x = temps_7)$y
x <- temps_7
dens <- function(x){density(x = x)$y}
nlm(f = dens, p = 15)


minus_dens <- function(x){-dens(x)}
minus_dens(problemtemps)
extrema <- c(nlm(f = minus_dens, p = 15)$estimate,
             nlm(f = minus_dens, p = 40)$estimate)

extrema



## Part e. ---- 
## Check with a kernel second derivative that the extreme points are actually modes.


kdde_2 <- ks::kdde(x = problemphones, deriv.order = 2)
plot(kdde_2, xlab = "x", main = "Density second derivative estimation")
abline(v = density(x = problemtemps)$x[pos_Max], col = "3")
abline(v = density(x = problemtemps)$x[posMax], col = "3")
abline(h = 0)






.





