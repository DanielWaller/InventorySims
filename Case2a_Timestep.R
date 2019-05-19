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
  
  # Use the first H forecast errors (they start indexed from H+L+1) to calculate the first SS
  # We are in period H + H + L now
  
  T <- H + H + L
  log.forecast.error.sample <- forecasts.base.errors[(T-H+1):T]
  log.safety.stock <- calculate.SS.CSL(forecast.errors = log.forecast.error.sample,CSL = CSL)
  # Initialise S, IP, OHS, ordering system
  
  S = exp(log.safety.stock + log(forecasts.base[(T+2)])) + forecasts.base[(T+1)] ; safety.stock <- S - forecasts.base[(T+1)] - forecasts.base[(T+2)]
  
  orders <- forecasts.base[(T+1):(T+L)] ; num.orders <- L ; periods.next <- 1# Assumption: R = 1
  order.times <- c(1:num.orders) ; Initial.OHS <- S - sum(orders) ; Initial.IP <- S
  Time.since.R <- 0 
  
  init.T <- T
  
  return(list(2, dataseries, promoinds, P, promoprop, R, L, H, forecasts.base, 
              forecasts.base.errors,orders, num.orders, periods.next, order.times,
              safety.stock, Initial.OHS, Initial.IP, Time.since.R,CSL,init.T))
  
}

#### (f) Simulation burn-in ####

burn.in.simulation <- function(datachunk,burnin.length){
  dataseries <- datachunk[[2]] ; promoinds <- datachunk[[3]] ; P <- datachunk[[4]] ; promoprop <- datachunk[[5]]
  R <- datachunk[[6]] ; L <- datachunk[[7]] ; H <- datachunk[[8]] ; forecasts.base <- datachunk[[9]]
  forecasts.base.errors <- datachunk[[10]] ; orders <- datachunk[[11]] ; num.orders <- datachunk[[12]]
  periods.next <- datachunk[[13]] ; order.times <- datachunk[[14]] ; safety.stock <- datachunk[[15]] ; OHS <- datachunk[[16]]
  Initial.IP <- datachunk[[17]] ; Time.since.R <- datachunk[[18]] ; CSL <- datachunk[[19]] ;
  T <- datachunk[[20]]
  
  N <- length(dataseries)
  B <- burnin.length
  
  lostdemand <- numeric(N)
  alphas <- numeric(N)
  betas <- numeric(N)
  IPs <- numeric(N)
  OHS.As <- numeric(N)
  OHS.Bs <- numeric(N)
  Ss <- numeric(N)
  ordersizes <- numeric(N)
  safetystocks <- numeric(N)
  
  Ss[T] <- Initial.IP
  IPs[T] <- Initial.IP
  OHS.Bs[T] <- OHS
  safetystocks[T] <- safety.stock
  
  # The burnin period is B periods long, starting from period H + H + L + 1, since:
  # H + H + L was the period we gained enough forecast errors for the sample
  # H + H + L + B is the burnin end. We will record this and pass it as a parameter
  
  B.end <- H + H + L + B
  
  T <- H + H + L + 1
  
  for(i in 1:B){
    
    # Observe demand
    data.fulfilled <- min(OHS, dataseries[T])
    lostdemand[T] <- min(OHS - dataseries[T],0) ; OHS <- max(OHS - dataseries[T], 0)
    
    OHS.As[T] <- OHS
    
    # Mark the forecast error
    
    forecasts.base.errors[T] <- log(dataseries[T]) - log(forecasts.base[T])
    
    # Receive any orders
    order.times <- order.times - 1
    if(order.times[1] == 0){
      OHS <- OHS + orders[1] ; orders <- orders[-1] ; order.times <- order.times[-1]
      OHS.Bs[T] <- OHS
    }
    
    # Consider a replenishment
    Time.since.R <- Time.since.R + 1
    if(R == Time.since.R){
      new.param.ests <- get.parameter.estimates(dataseries = dataseries[(T-H+1):T],promoinds = promoinds[(T-H+1):T],P = P)
      alpha <- new.param.ests[[1]] ; beta <- new.param.ests[[2]] ; alphas[T] <- alpha ; betas[T] <- beta
      
      # Make the new forecast
      
      p.forecast <- get.forecast.base(alpha.est = alpha, beta.est = beta, promoind = promoinds[(T+L+1)], price.cut = P)
      forecasts.base[(T+L+1)] <- p.forecast
      
      # Calculate the safety stock
      
      log.safety.stock <- calculate.SS.CSL(forecast.errors = forecasts.base.errors[(T-H+1):T],CSL = CSL)
      
      # Calculate S
      
      S = exp(log.safety.stock + log(forecasts.base[(T+2)])) + forecasts.base[(T+1)] ; safety.stock <- S - forecasts.base[(T+1)] - forecasts.base[(T+2)]
      
      IP <- OHS + sum(orders) 
      order.size <- max(S - IP,0) ; ordersizes[T] <- order.size
      t <- length(order.times) ; orders[[(t+1)]] <- order.size ; order.times[(t+1)] <- L
      Time.since.R <- 0
    
      
      # Update vectors
      
      Ss[T] <- S ;  IPs[T] <- IP ; safetystocks[T] <- safety.stock
      
    } 
    
    # Advance time
      
    T <- T+1
    
  }
  
  datamass <- list(2, dataseries, promoinds, P, promoprop, R, L, H, forecasts.base, 
                   forecasts.base.errors,orders, num.orders, periods.next, order.times,
                   Time.since.R,CSL,lostdemand,alphas,betas,IPs,OHS.As,
                   OHS.Bs,Ss,ordersizes,safetystocks,B.end)
  
  return(datamass)
  
}










####

datalist <- DGP_2(baseline = 100, sigma = 0.1, Length = 150, price.cut = 0.5, promoprop = 0.1, elasticity = -4)
datachunk <- initialise.inventory.sim(datalist, 1,1,0.95,20)
datamass <- burn.in.simulation(datachunk, burnin.length = 59)