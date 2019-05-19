# Case 2 - Promotions #

#### (a) Specific DGP (we call it model 2) ####

DGP_2 <- function(baseline, sigma, Length, price.cut, promoprop, elasticity){
  
  errors <- exp(rnorm(mean = 0, sd = sigma, n = Length))
  beta <- elasticity ; P <- price.cut
  promoinds <- numeric(Length)
  for(j in 1:(floor(Length*promoprop))){
    promoinds[(j*(1/promoprop))] <- 1
  }
  dataseries <- baseline * ((P^beta)^promoinds) * errors
  
  return(list(2,dataseries,promoinds,P,promoprop))
  
}

#### (b) Estimating baseline and elasticity ####

get.parameter.estimates <- function(dataseries, promoinds, P){
  N <- length(dataseries) ; log.data <- log(dataseries) ; regressors <- promoinds * log(P)
  linear.model <- lm(log.data ~ regressors)
  alpha.est <- exp(unname(linear.model$coefficients[1])) ; beta.est <- unname(linear.model$coefficients[2]) 
  return(list(alpha.est, beta.est))
}

#### (c) Forecast functions based on parameter estimates ####

####     (c.i) Base forecast function ####

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
  R <- review.period ; L <- lead.time ; H <- History
  dataseries <- datalist[[2]] ; promoinds <- datalist[[3]] ; P <- datalist[[4]]
  promoprop <- datalist[[5]]
  
  forecasts.base <- numeric()
  forecasts.base.errors <- numeric()
  
  # First H forecasts produced at end of Timesteps (H through (H + (H - 1)) )
  # The forecasts are recorded as forecasts for Period (Timestep + L + 1)
  
  T <- H
  
  for(j in 1:H){
    datasample <- dataseries[(T-H+1):T] ; promoindssample <- promoinds[(T-H+1):(T)]
    param.ests <- get.parameter.estimates(dataseries = datasample, promoinds = promoindssample, P)
    alpha.est <- param.ests[[1]] ; beta.est <- param.ests[[2]] ; promoind <- promoinds[(T+L+1)]
    forecasts.base[(T+L+1)] <- get.forecast.base(alpha.est = alpha.est, beta.est = beta.est,
                                           promoind = promoind, price.cut = P)
    forecasts.base.errors[(T+L+1)] <- log(dataseries[(T+L+1)]) - log(forecasts.base[(T+L+1)])
    T <- T + 1
  }
  
  # Forecasts produced whilst waiting for first H forecast errors to be completed
  # These are produced during periods (H + H) through (H + H + L)
  # (H + H) is the period where the first forecast is produced not needed for the initial forecast error sample
  # (H + H + L) is the period where we observe the data forecasted in (H + (H - 1))
  # The forecasts are recorded as forecasts for Period (T + L + 1), as usual.
  
  T <- H+H
  
  for(j in 1:(1+L)){
    datasample <- dataseries[(T-H+1):T] ; promoindssample <- promoinds[(T-H+1):T]
    param.ests <- get.parameter.estimates(dataseries = datasample, promoinds = promoindssample, P)
    alpha.est <- param.ests[[1]] ; beta.est <- param.ests[[2]] ; promoind <- promoinds[(T+L+1)]
    forecasts.base[(T+L+1)] <- get.forecast.base(alpha.est = alpha.est, beta.est = beta.est,
                                           promoind = promoind, price.cut = P)
    T <- T + 1
  }
  
  # Use the first H forecast errors (they start indexed from 2+L) to calculate SS
  
  log.forecast.sample <- forecasts.base.errors[(2+L):(2+L+H-1)]
  log.safety.stock <- calculate.SS.CSL(forecast.errors = log.forecast.sample,CSL = CSL)
  # Initialise S, IP, OHS, ordering system
  
  S = exp(log.safety.stock + log(forecasts.base[(H+1)])) ; safety.stock <- S - forecasts.base[(H+1)]
  
  orders <- forecasts.base[(H+2):(H+1+L)] ; num.orders <- L ; periods.next <- 1# Assumption: R = 1
  order.times <- c(1:num.orders)
  Initial.IP <- S ; Initial.OHS <- S - sum(orders)
  Time.since.R <- 0 
  
  return(list(2, dataseries, promoinds, P, promoprop, R, L, H, forecasts.base, 
              forecasts.base.errors,orders, num.orders, periods.next, order.times,
              Initial.IP, Initial.OHS, Time.since.R))
  
}