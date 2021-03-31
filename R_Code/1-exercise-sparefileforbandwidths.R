library(shiny)
library(nor1mix)
library(ks)


bw.ucv.mod <- function(x, nb = 1000L,
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
  fucv <- function(h) .Call(stats:::C_bw_ucv, n, d, cnt, h)
  # h <- optimize(fucv, c(lower, upper), tol = tol)$minimum
  # if (h < lower + tol | h > upper - tol)
  #   warning("minimum occurred at one end of the range")
  obj <- sapply(h.grid, function(h) fucv(h))
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

bw.bcv.mod(x = pasttemps, nb = 1000L, h.grid = 10^seq(-3, log10(1.2 * sd(x) * length(x)^(-1/5)), l = 200))


bw.ucv.mod <- function(x, nb = 1000L,
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
  fucv <- function(h) .Call(stats:::C_bw_ucv, n, d, cnt, h)
  # h <- optimize(fucv, c(lower, upper), tol = tol)$minimum
  # if (h < lower + tol | h > upper - tol)
  #   warning("minimum occurred at one end of the range")
  obj <- sapply(h.grid, function(h) fucv(h))
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

bw.ucv.mod(x = pasttemps, nb = 512)


# four types of bandwidth selectors

bw.nrd(x = pasttemps) #unrealistic
hpi(x = pasttemps, binned = TRUE, bgridsize = 512)
bw.ucv.mod(x = pasttemps, nb = 512) #problem
bw.bcv.mod(x = pasttemps  , nb = 512)



