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

p1 <- hist(problemPhones, xlim=c(0,65), probability = T, col=rgb(1,0,0,1/4), main = "Problematic Phones")
p2 <- hist(pastPhones, xlim=c(0,65), probability = T, col=rgb(0,0,1,1/4), main = "Past Phones")

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

# !Problematic temperatures =================================

#LSCV

bw.lscv.problemtemps <- bw.ucv(x = problemPhones)

#DPI

bw_dpi_problemtemps <- bw.SJ(x = problemPhones, method = "dpi")

# !Past temperatures =================================

#LSCV

bw.lscv.pasttemps <- bw.ucv(x = pastPhones)

#DPI

bw_dpi_pasttemps <- bw.SJ(x = pastPhones, method = "dpi")


#PLOTS

# Problem Temperatures

par(mfrow=c(1,2))
plot(density(x = problemtemps, bw = bw.lscv.problemtemps), xlim = c(-5,65)) 
plot(density(x = problemtemps, bw = bw_dpi_problemtemps), xlim = c(-5,65))
par(mfrow=c(1,1))


# Past Temperatures

par(mfrow=c(1,2))
plot(density(x = pasttemps, bw = bw.lscv.pasttemps), xlim = c(-5,65))
plot(density(x = pasttemps, bw = bw_dpi_pasttemps), xlim = c(-5,65))
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

# !Problematic temperatures =================================

density(x = problemPhones)$x[posMax1]  
density(x = problemPhones)$x[posMax2]

density(x = problemPhones)$y[posMax1]
density(x = problemPhones)$y[posMax2]  

# !Problematic temperatures =================================

density(x = pastPhones)$y[posMax1]
density(x = pastPhones)$y[posMax2]  

density(x = pastPhones)$x[posMax1]  
density(x = pastPhones)$x[posMax2]


## Part e. ---- 
## Check with a kernel second derivative that the extreme points are actually modes.

# !Problematic temperatures =================================
kdde_2_problemPhones <- ks::kdde(x = problemPhones, deriv.order = 2)

par(mfrow = c(1,3))
plot(kdde_0_problemPhones, xlab = "x", main = "Problematic battery: \n Density estimation", xlim = c(-5,65))
abline(v = density(x = problemPhones)$x[posMax1], col = "3")
abline(v = density(x = problemPhones)$x[posMax2], col = "3")
abline(v = density(x = problemPhones)$x[posMin], col = "2")
plot(kdde_1_problemPhones, xlab = "x", main = "Problematic battery: \n Density derivative estimation")
abline(v = density(x = problemPhones)$x[posMax1], col = "3")
abline(v = density(x = problemPhones)$x[posMax2], col = "3")
abline(v = density(x = problemPhones)$x[posMin], col = "2")
abline(h = 0)
plot(kdde_2_problemPhones, xlab = "x", main = "Problematic battery: \n  Density second derivative estimation")
abline(v = density(x = problemPhones)$x[posMax1], col = "3")
abline(v = density(x = problemPhones)$x[posMax2], col = "3")
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

