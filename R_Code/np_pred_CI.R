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