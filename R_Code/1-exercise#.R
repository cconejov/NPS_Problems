#libraries

library(tidyverse)

# load data

problemphones <- read.table("temps-7.txt", header = TRUE)
pastphones <- read.delim("temps-other.txt", header = TRUE)


#rename variables

names(problemphones)[names(problemphones) == "x"] <- "problemtemps"
names(pastphones)[names(pastphones) == "x"] <- "pasttemps"

## maintain this order of temperatures of problematic phones and the temperatires of past phones 


#attaching for ease of use

attach(problemphones)
attach(pastphones)

#sumamries of data

summary(problemtemps)
summary(pastphones)


#histograms

p1 <- hist(problemtemps)
p2 <- hist(pasttemps)

plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,65))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,65), add=T)  

#pastemps in red, so lower, makes since they didn't burn


#probability histograms


p1p <- hist(problemtemps, probability = T)
p2p <- hist(pasttemps, probability = T)

plot( p1p, col=rgb(0,0,1,1/4), xlim=c(0,65))  # first histogram
plot( p2p, col=rgb(1,0,0,1/4), xlim=c(0,65), add=T)  

#looks very similar
#could adjust bandwidth but not necessary i believe

#KDE and Density

#doesn't work very well

plot(density(problemtemps), ylim = c(0, 0.8))
curve(dlnorm(x), from = -2, to = 10, n = 500, col = 2, add = TRUE)
rug(problemtemps)

#log transform - much much better - but why?

kde <- density(log(problemtemps))
plot(kde, main = "Kde of transformed data")
rug(log(problemtemps))

# I think we literally need to use the KS library

## KS only considers normal kernels so I need to use other stuff





# BANDWIDTH SELECTION


#OPTIONS: AMISE



# DPI (?) - Rule of thumb

# Rule-of-thumb - not useful if i remember correctly
bw.nrd(x = problemtemps)
## [1]0.8797933
# bwd.nrd employs 1.34 as an approximation for diff(qnorm(c(0.25, 0.75))) - (?)

# Same as
iqr <- diff(quantile(problemtemps, c(0.25, 0.75))) / diff(qnorm(c(0.25, 0.75)))
1.06 * 9947^(-1/5) * min(sd(problemtemps), iqr)
## [1] 0.873937


# DPI


# DPI selector
bw.SJ(x = problemtemps, method = "dpi")
## [1] 0.6303672 ,but why so different?

# Similar to
ks::hpi(x) # Default is two-stage
## [1] 0.4999456 - gives most different result


### LSCV ### 

# UCV gives a warning
bw.ucv(x = problemtemps)
## [1] 0.6524331

# Extend search interval
bw.ucv(x = problemtemps, lower = 0.01, upper = 1)
## [1] 0.6514426


## DOESN'T WORK

# bw.ucv.mod replaces the optimization routine of bw.ucv with an exhaustive
# search on "h_grid" (chosen adaptatively from the sample) and optionally
# plots the LSCV curve with "plot_cv"

bw.ucv.mod <- function(problemphones, nb = 1000L,
                       h_grid = 10^seq(-3, log10(1.2 * sd(problemphones) *
                                                   length(problemphones)^(-1/5)), l = 200),
                       plot_cv = FALSE) {
  if ((n <- length(problemphones)) < 2L)
    stop("need at least 2 data points")
  n <- as.integer(n)
  if (is.na(n))
    stop("invalid length(x)")
  if (!is.numeric(problemphones))
    stop("invalid 'x'")
  nb <- as.integer(nb)
  if (is.na(nb) || nb <= 0L)
    stop("invalid 'nb'")
  storage.mode(problemphones) <- "double"
  hmax <- 1.144 * sqrt(var(problemphones)) * n^(-1/5)
  Z <- .Call(stats:::C_bw_den, nb, problemphones)
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

# Compute the bandwidth and plot the LSCV curve
bw.ucv.mod(x = problemphones, plot_cv = TRUE, h_grid = 10^seq(-1.25, 0.5, l = 200))
## [1] 0.5431561

# We can compare with the default bw.ucv output
abline(v = bw.ucv(problemphones), col = 3)


## DOESN'T WORK


# USING KS PACKAGE - NORMAL KERNELS


bw <- 0.6514426
plot(kde <- ks::kde(x = problemphones, h = bw), lwd = 3) # ?ks::plot.kde for options
lines(density(x = problemphones, bw = bw), col = 2)






# b) theory I think LSCV has similar variance but less bias if I can remember correctly

# c) don't get how this changes that much - need to verify this

# d) 

# e)







