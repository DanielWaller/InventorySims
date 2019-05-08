# Case 2 - Promotions #

#### (a) Specific DGP (we call it model 2) ####

DGP_2 <- function(baseline, sigma, Length, price.cut, promoprop, elasticity){
  
  errors <- exp(rnorm(mean = 0, sd = sigma, n = Length))
  beta <- elasticity ; P <- price.cut
  promoinds <- numeric(Length)
  for(j in 1:(floor(Length*promoprop))){
    promoinds[(j*20)] <- 1
  }
  dataseries <- baseline * ((P^beta)^promoinds) * errors
  
  return(list(2,dataseries,promoinds,P,beta))
  
}

#### (b) Estimating baseline and elasticity ####

get.parameter.estimates <- function(dataseries, promoinds, P){
  N <- length(dataseries) ; log.data <- log(dataseries) ; regressors <- promoinds * log(P)
  linear.model <- lm(log.data ~ regressors)
  alpha.est <- unname(linear.model$coefficients[1]) ; beta.est <- unname(linear.model$coefficients[2]) 
  return(list(alpha.est, beta.est))
}

#### (c) Forecast functions based on parameter estimates ####

#### (i) Base forecast function ####

get.forecast.base <- function(alpha.est, beta.est, promoind, price.cut){
  P <- price.cut
  forecast <- alpha.est * ((P^beta.est)^promoind)
  return(forecast)
}

#### (d) Calculate SS based on forecast errors ####

calculate.SS.CSL <- function(forecast.errors, CSL){
  sorted.errors <- sort(forecast.errors, decreasing = FALSE)
  N <- length(forecast.errors) ; quant <- floor(CSL * N)
  SS.quant <- sorted.errors[quant]
  return(SS.quant)
}

#### (e) Initialise the simulations ####

initialise.inventory.sim <- function(datalist,review.period,lead.time,CSL,History){
  R <- review.period ; L <- lead.time : H <- History
  dataseries <- datalist[[2]] ; promoinds <- datalist[[3]] ; P <- datalist[[4]]
  promoprop <- datalist[[5]]
  
  forecasts.base <- numeric()
  forecasts.base.errors <- numeric()
  
  # First H forecasts
  
  for(j in 1:H){
    datasample <- dataseries[j:(H-1+j)] ; promoindssample <- promoinds[j:(H-1+j)]
    param.ests <- get.parameter.estimates(dataseries = datasample, promoinds = promoindssample, P)
    alpha.est <- paramests[[1]] ; beta.est <- paramests[[2]] ; promoind <- promoinds[(j+H-1+L)]
    forecasts.base[j] <- get.forecast.base(alpha.est = alpha.est, beta.est = beta.est,
                                           promoind = promoind, price.cut = P)
    forecasts.base.errors[(j+1+L)] <- dataseries[(j+1+L)] - forecasts.base[j]
  }
  
  # Forecasts produced whilst waiting for first H forecast errors
  
  for(j in (H+1):(H+1+L)){
    datasample <- dataseries[j:(H-1+j)] ; promoindssample <- promoinds[j:(H-1+j)]
    param.ests <- get.parameter.estimates(dataseries = datasample, promoinds = promoindssample, P)
    alpha.est <- paramests[[1]] ; beta.est <- paramests[[2]] ; promoind <- promoinds[(j+H-1+L)]
    forecasts.base[j] <- get.forecast.base(alpha.est = alpha.est, beta.est = beta.est,
                                           promoind = promoind, price.cut = P)
  }
  
  # Use the log of the first H forecast errors (they start indexed from 2+L) to calculate SS
  
  log.forecast.sample <- log(forecasts.base.errors[(2+L):(2+L+H-1)])
  log.safety.stock <- calculate.SS.CSL(forecast.errors = log.forecast.sample,CSL = CSL)
  safety.stock <- exp(log.safety.stock)
  
  # Initialise S, IP, OHS, ordering system
  
  S = safety.stock + forecasts.base[(H+1)]
  
  orders <- forecasts.base[(H+2):(H+1+L)] ; num.orders <- L ; periods.next <- 1# Assumption: R = 1
  Initial.IP <- S[1] ; Initial.OHS <- S[1] - (order.size * num.orders)
  Time.since.R <- 0 ; k <- S[4]
  order.times <- c(1:num.orders)
 
  
  
}

#### (f) Simulation burn-in ####




####

datalist <- DGP_2(baseline = 100, sigma = 0.4, Length = 20, price.cut = 0.5, promoprop = 0.05, elasticity = -4)
dataseries <- datalist[[2]] ; promoinds <- datalist[[3]] ; P <- datalist[[4]]

