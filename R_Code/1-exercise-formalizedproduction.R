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

par(mfrow = c(1,2))
p1 <- hist(problemPhones, xlim=c(0,65))
p2 <- hist(pastPhones, xlim=c(0,65))
par(mfrow = c(1,1))

plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,65), main = "Comparison Histograms")  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,65), add=T)  
rm(p1,p2)

# Density plots - default bandwidths

par(mfrow = c(1,2))
plot(density(x = problemPhones), xlim = c(-5,65), main = 'Problematic battery')
plot(density(x = pastPhones), xlim = c(-5,65), main = 'Past battery')
par(mfrow = c(1,1))


## Part a. ----
## Perform a kernel density estimation for temps-7 and temps-other using what you consider is the most adequate bandwidth. 
## Since the temperatures are positive, is it required to perform any transformation?.

# In Preliminary Work we see how the default bandwidth bw.nrd0 discover the properties of both sample distributions:

bw.nrd0(x = problemPhones)
bw.nrd0(x = pastPhones)

# TRANSFORMATIONS

# Transformations must be performed in order to avoid boundary bias, which is to assign probability, 
# with the function mapping, when there is no real density to support it.
# Hence, we are looking at temperatures that theoretically can take values below zero or close to zero to  a very  high values that can be considered
# positive infinity. In reality, in our sample, does not have any values close to zero, therefore boundary bias is a risk. From this, we can deduct that a log
# transformation could be possible but we will check if this is necesaary in our sample. 

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
plot(kde_transf_problemPhones, main = "Problematic battery: Transformed kde", xlim = c(-5, 65))
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
plot(kde_transf_pastPhones, main = "Past battery: Transformed kde", xlim = c(-5, 60))
par(mfrow = c(1,1))

head(kde_transf_pastPhones)


# Although the KDE does not assign it in the negative part of the support, the transform data
# smooths the KDE. Also, we notice how the size of the bandwidth is smaller in the transformed data.

kde_transf_problemPhones$bw
kde_transf_pastPhones$bw

## Part b. ----
## Is there any important difference on the results from considering the LSCV selector over the DPI selector?


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

## Part c. ----
## It seems that in temps-7 there is a secondary mode. Compute a kernel derivative estimation for temps-7 and temps-other using what you consider are the most adequate bandwidths.

# !Problematic temperatures =================================

## Detect Global mode
pos_Max <- which.max(density(x = problemPhones)$y)

## Detect minima
YY <- density(x = problemPhones)$y[density(x = problemPhones)$x > 20 & density(x = problemPhones)$x < 40]
minY <- min(YY)
posMin <- which(density(x = problemPhones)$y == minY)

## Detect Secondary mode
YY <- density(x = problemPhones)$y[density(x = problemPhones)$x > 33 & density(x = problemPhones)$x < 60]
maxY <- max(YY)
posMax <- which(density(x = problemPhones)$y == maxY)


kdde_0_problemPhones <- ks::kdde(x = problemPhones, deriv.order = 0)
kdde_1_problemPhones <- ks::kdde(x = problemPhones, deriv.order = 1)

par(mfrow = c(1,2))
plot(kdde_0_problemPhones, xlab = "x", main = "Problematic battery: Density estimation", xlim = c(-5,65))
abline(v = density(x = problemPhones)$x[pos_Max], col = "3")
abline(v = density(x = problemPhones)$x[posMax], col = "3")
abline(v = density(x = problemPhones)$x[posMin], col = "2")
plot(kdde_1_problemPhones, xlab = "x", main = "Problematic battery: Density derivative estimation")
abline(v = density(x = problemPhones)$x[pos_Max], col = "3")
abline(v = density(x = problemPhones)$x[posMax], col = "3")
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



## Part e. ---- 
## Check with a kernel second derivative that the extreme points are actually modes.

# !Problematic temperatures =================================
kdde_2_problemPhones <- ks::kdde(x = problemPhones, deriv.order = 2)

par(mfrow = c(1,3))
plot(kdde_0_problemPhones, xlab = "x", main = "Problematic battery: \n Density estimation", xlim = c(-5,65))
abline(v = density(x = problemPhones)$x[pos_Max], col = "3")
abline(v = density(x = problemPhones)$x[posMax], col = "3")
abline(v = density(x = problemPhones)$x[posMin], col = "2")
plot(kdde_1_problemPhones, xlab = "x", main = "Problematic battery: \n Density derivative estimation")
abline(v = density(x = problemPhones)$x[pos_Max], col = "3")
abline(v = density(x = problemPhones)$x[posMax], col = "3")
abline(v = density(x = problemPhones)$x[posMin], col = "2")
abline(h = 0)
plot(kdde_2_problemPhones, xlab = "x", main = "Problematic battery: \n  Density second derivative estimation")
abline(v = density(x = problemPhones)$x[pos_Max], col = "3")
abline(v = density(x = problemPhones)$x[posMax], col = "3")
abline(v = density(x = problemPhones)$x[posMin], col = "2")
abline(h = 0)
par(mfrow = c(1,1))



# !Past temperatures =================================

kdde_2_pastPhones <- ks::kdde(x = pastPhones, deriv.order = 2)

par(mfrow = c(1,3))
plot(kdde_0_pastPhones, xlab = "x", main = "Past battery \n Density estimation", xlim = c(-5,65))
abline(v = density(x = pastPhones)$x[pos_Max_PastPhone], col = "3")
plot(kdde_1_pastPhones, xlab = "x", main = "Past battery \n Density derivative estimation")
abline(v = density(x = pastPhones)$x[pos_Max_PastPhone], col = "3")
abline(h = 0)
plot(kdde_2_pastPhones, xlab = "x", main = "Past battery \n Density Second derivative estimation")
abline(v = density(x = pastPhones)$x[pos_Max_PastPhone], col = "3")
abline(h = 0)
par(mfrow = c(1,1))

