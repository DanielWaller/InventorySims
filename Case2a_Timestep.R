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
  sdalpha <- summary(linear.model)$coefficients[1,2]
  s <- summary(linear.model)[[6]]
  return(list(alpha.est, beta.est,sdalpha,s))
}

#### (c) Forecast functions based on parameter estimates ####

####     (c.i) Base forecast function ####

get.forecast.base <- function(alpha.est, beta.est, promoind, price.cut){
  P <- price.cut
  forecast <- alpha.est * ((P^beta.est)^promoind)
  return(forecast)
}

####     (c.ii) Miller's approximation ####

get.forecast.miller <- function(alpha.est, beta.est, promoind, price.cut, s){
  P <- price.cut 
  forecast <- alpha.est * ((P^beta.est)^promoind) * exp((s^2)/2)
  return(forecast)
}

####     (c.iii) New approximation ####

get.forecast.approx <- function(alpha.est, beta.est, promoind, price.cut, s, sdalpha,promoinds.hist){
  P <- price.cut 
  P.vector <- P^promoinds.hist ; sumPvec.sq <- sum(P.vector * P.vector)
  forecast <- alpha.est * ((P^beta.est)^promoind) * exp(((s^2)/2) * (1 - ((log(P^promoind) * log(P^promoind))/sumPvec.sq) )) * exp(-1*(sdalpha^2)/2)
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

initialise.inventory.sim <- function(datalist,review.period,lead.time,CSL,History,fcastmethod){
  R <- review.period ; L <- lead.time ; H <- History
  dataseries <- datalist[[2]] ; promoinds <- datalist[[3]] ; P <- datalist[[4]]
  promoprop <- datalist[[5]] ; fcastmethod <- fcastmethod
  
  forecasts.base <- numeric()
  forecasts.base.errors <- numeric()
  
  # First H forecasts produced at end of Timesteps (H through (H + (H - 1)) )
  # The forecasts are recorded as forecasts for Period (Timestep + L + 1)
  
  T <- H
  
  for(j in 1:H){
    datasample <- dataseries[(T-H+1):T] ; promoindssample <- promoinds[(T-H+1):(T)]
    param.ests <- get.parameter.estimates(dataseries = datasample, promoinds = promoindssample, P)
    alpha.est <- param.ests[[1]] ; beta.est <- param.ests[[2]] ; sdalpha <- param.ests[[3]]
    s <- param.ests[[4]] ; promoind <- promoinds[(T+L+1)]
    if(fcastmethod == 1){
      forecasts.base[(T+L+1)] <- get.forecast.base(alpha.est = alpha.est, beta.est = beta.est,
                                                   promoind = promoind, price.cut = P)
    }
    if(fcastmethod == 2){
      forecasts.base[(T+L+1)] <- get.forecast.miller(alpha.est = alpha.est, beta.est = beta.est,
                                                     promoind = promoind, price.cut = P, s = s)
    }
    if(fcastmethod == 3){
      forecasts.base[(T+L+1)] <- get.forecast.approx(alpha.est = alpha.est, beta.est = beta.est,
                                                     promoind = promoind, price.cut = P, s = s,
                                                     sdalpha = sdalpha, promoinds.hist = promoinds[(T-H+1):(T)])
    }
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
    alpha.est <- param.ests[[1]] ; beta.est <- param.ests[[2]] ; sdalpha <- param.ests[[3]]
    s <- param.ests[[4]] ; promoind <- promoinds[(T+L+1)]
    if(fcastmethod == 1){
      forecasts.base[(T+L+1)] <- get.forecast.base(alpha.est = alpha.est, beta.est = beta.est,
                                                   promoind = promoind, price.cut = P)
    }
    if(fcastmethod == 2){
      forecasts.base[(T+L+1)] <- get.forecast.miller(alpha.est = alpha.est, beta.est = beta.est,
                                                     promoind = promoind, price.cut = P, s = s)
    }
    if(fcastmethod == 3){
      forecasts.base[(T+L+1)] <- get.forecast.approx(alpha.est = alpha.est, beta.est = beta.est,
                                                     promoind = promoind, price.cut = P, s = s, 
                                                     sdalpha = sdalpha, promoinds.hist = promoinds[(T-H+1):(T)])
    }
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
              safety.stock, Initial.OHS, Initial.IP, Time.since.R,CSL,init.T,fcastmethod))
  
}

#### (f) Simulation burn-in ####

burn.in.simulation <- function(datachunk,burnin.length){
  dataseries <- datachunk[[2]] ; promoinds <- datachunk[[3]] ; P <- datachunk[[4]] ; promoprop <- datachunk[[5]]
  R <- datachunk[[6]] ; L <- datachunk[[7]] ; H <- datachunk[[8]] ; forecasts.base <- datachunk[[9]]
  forecasts.base.errors <- datachunk[[10]] ; orders <- datachunk[[11]] ; num.orders <- datachunk[[12]]
  periods.next <- datachunk[[13]] ; order.times <- datachunk[[14]] ; safety.stock <- datachunk[[15]] ; OHS <- datachunk[[16]]
  Initial.IP <- datachunk[[17]] ; Time.since.R <- datachunk[[18]] ; CSL <- datachunk[[19]] ;
  T <- datachunk[[20]] ; fcastmethod <- datachunk[[21]]
  
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
      alpha <- new.param.ests[[1]] ; beta <- new.param.ests[[2]] ; sdalpha <- new.param.ests[[3]]
      s <- new.param.ests[[4]] ; alphas[T] <- alpha ; betas[T] <- beta
      
      # Make the new forecast
      
      if(fcastmethod == 1){
        p.forecast <- get.forecast.base(alpha.est = alpha, beta.est = beta, promoind = promoinds[(T+L+1)], price.cut = P)
      }
      if(fcastmethod == 2){
        forecasts.base[(T+L+1)] <- get.forecast.miller(alpha.est = alpha.est, beta.est = beta.est,
                                                       promoind = promoind, price.cut = P, s = s)
      }
      if(fcastmethod == 3){
        forecasts.base[(T+L+1)] <- get.forecast.approx(alpha.est = alpha.est, beta.est = beta.est,
                                                       promoind = promoind, price.cut = P, s = s,
                                                       sdalpha = sdalpha, promoinds.hist = promoinds[(T-H+1):(T)])
      } 
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
                   OHS.Bs,Ss,ordersizes,safetystocks,B.end,fcastmethod)
  
  return(datamass)
  
}

#### (g) Test period of the simulation ####

simulation.test.period <- function(datamass){
  
  dataseries <- datamass[[2]] ; promoinds <- datamass[[3]] ; P <- datamass[[4]] ; promoprop <- datamass[[5]]
  R <- datamass[[6]] ; L <- datamass[[7]] ; H <- datamass[[8]] ; forecasts.base <- datamass[[9]]
  forecasts.base.errors <- datamass[[10]] ; orders <- datamass[[11]] ; num.orders <- datamass[[12]]
  periods.next <- datamass[[13]] ; order.times <- datamass[[14]] ; Time.since.R <- datamass[[15]]
  CSL <- datamass[[16]] ; lostdemand <- datamass[[17]] ; alphas <- datamass[[18]]
  betas <- datamass[[19]] ; IPs <- datamass[[20]] ; OHS.As <- datamass[[21]] ; OHS.Bs <- datamass[[22]]
  Ss <- datamass[[23]] ; ordersizes <- datamass[[24]] ; safetystocks <- datamass[[25]] ; B.end <- datamass[[26]]
  fcastmethod <- datamass[[27]]
  
  # We are now at B.end + 1 - calculate distance until end
  
  T <- B.end + 1 ; Test.length <- length(dataseries) - B.end
  
  for(i in 1:Test.length){
    
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
      alpha <- new.param.ests[[1]] ; beta <- new.param.ests[[2]] ; sdalpha <- new.param.ests[[3]]
      s <- new.param.ests[[4]] ; alphas[T] <- alpha ; betas[T] <- beta
      
      # Make the new forecast
      
      if(fcastmethod == 1){
        p.forecast <- get.forecast.base(alpha.est = alpha, beta.est = beta, promoind = promoinds[(T+L+1)], price.cut = P)
      }
      if(fcastmethod == 2){
        forecasts.base[(T+L+1)] <- get.forecast.miller(alpha.est = alpha.est, beta.est = beta.est,
                                                       promoind = promoind, price.cut = P,s = s)
      }
      if(fcastmethod == 3){
        forecasts.base[(T+L+1)] <- get.forecast.approx(alpha.est = alpha.est, beta.est = beta.est,
                                                       promoind = promoind, price.cut = P, s = s,
                                                       sdalpha = sdalpha, promoinds.hist = promoinds[(T-H+1):(T)])
      }
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
  
  # Things to return
  
  replication.output <- list(dataseries, promoinds, forecasts.base, forecasts.base.errors,
                             lostdemand, alphas, betas, IPs, OHS.As, OHS.Bs, Ss,
                             ordersizes, safetystocks,T)
  
  
}

#### 2. Run simulations ####

set.seed(5006)  
list.outputs <- list()
coverages <- numeric(500)
  
for(z in 1:500){
  datalist <- DGP_2(baseline = 100, sigma = 1, Length = 202, price.cut = 0.5, 
                    promoprop = 0.1, elasticity = -4)
  datachunk <- initialise.inventory.sim(datalist, lead.time = 1,review.period = 1,CSL = 0.75,
                                        History = 20,fcastmethod = 2)
  datamass <- burn.in.simulation(datachunk, burnin.length = 59)
  output <- simulation.test.period(datamass)
  list.outputs[[z]] <- output
  coverages[z] <- 1 - (sum(output[[5]][101:200] < 0)/100)
  if(z%%25 == 0){
    print(z)
  }
}

#### Data dictionary ####

# Ouptuts:
# 1. dataseries
# 2. promoinds
# 3. forecasts.base
# 4. forecasts.base.errors
# 5. lostdemand
# 6. alphas
# 7. betas
# 8. IPs
# 9. OHS.As
# 10. OHS.Bs
# 11. Ss
# 12. ordersizes
# 13. safetystocks
# 14. T

#### Plots and data storing ####
#### Process plot with last outputs: ####
output.curr <- list.outputs[[500]]
output.curr10amend <- output.curr[[10]]
output.curr10amend[is.na(output.curr10amend)] <- 0
OHS.matrix <- matrix(c(output.curr[[9]],output.curr10amend),ncol = 202, nrow = 2, byrow = TRUE)
topline <- as.vector(OHS.matrix)
seq1 <- rep(1:200, each = 2)

plot(x = seq1, y = topline[1:400],pch = 16, 
     ylim = c(-400, 2200),
     xlab = "Time period", ylab = "On-hand-stock",
     main = "Inventory process - Multiplicative promo. case")
abline(h = 0, col = "red")
for(i in 1:200){
  lines(x = seq1[(2*i):(2*i + 1)], y = topline[(2*i):(2*i + 1)],lty = 2)
}
for(i in 1:200){
  lines(x = seq1[(2*i - 1):(2*i)], y = topline[(2*i - 1):(2*i)],lty = 1)
}
points(x = which(output.curr[[5]] < 0), y = output.curr[[5]][output.curr[[5]] < 0 ], pch = 18, col = "blue")

pts1 <- which(output.curr[[5]] < 0) ; pts2 <- output.curr[[5]][output.curr[[5]] < 0]

for(i in 1:(length(pts1))){
  lines(x = rep(pts1[i],2), y = c(0, pts2[i]),col = "blue")
}

legend("topleft",legend = c("OHS","Lost sales"),pch = c(16,18), col = c(1,4))

abline(v = c(41,100),lty = 2, col = "red")

text(x = 20, y = -200, labels = "Initialisation",col = 3, cex = 1.25)
text(x = 70, y = -200, labels = "Burn-in",col = 3, cex = 1.25)
text(x = 155, y = -200, labels = "Test period",col = 3, cex = 1.25)

#### Coverage plots for 75%, 85%, 95% CSLs with increasing sigma = 0.1,0.3,0.5 ####

sigma1 <- numeric(3) ; sigma2 <- numeric(3) ; sigma3 <- numeric(3) ; sigma4 <- numeric(3)

sigma4[3] <- mean(coverages)

# Plot

frseq <- c(0.75,0.85,0.95)
plot(x = frseq, y = sigma1, ylim = c(0.7,1), xlim = c(0.75, 1), pch = 15, cex = 1.5, lwd = 1.5, type = "o",col = "red",
     xlab = "CSL (theoretical)", ylab = "CSL (empirical)", main = "Base forecast CSL (Emp. vs. Theor.) - varying sigma")
abline(a = 0, b = 1, col = "black",lty = 2)
lines(x = frseq, y = sigma2, pch = 16, cex=  1.25, lwd = 2, col = "green",type= "o")
lines(x = frseq, y = sigma3, pch = 17, cex=  1.25, lwd = 2, col = "blue",type= "o")
lines(x = frseq, y = sigma4, pch = 18, cex=  1.25, lwd = 2, col = 5,type= "o")

legend("bottomright", legend = c("sigma = 0.1","sigma = 0.3","sigma = 0.5","sigma = 1"),pch = c(15,16,17,18),
       col = c(2,3,4,5),cex = 1.5)

#### Forecast error distributions ####

base.errors.03 <- matrix(NA, ncol = 100, nrow = 500)
for(i in 1:500){
  base.errors.03[i,] <- list.outputs[[i]][[4]][101:200]
}

hist(base.errors.03, main = "Log forecast errors (sigma = 0.3)")

promoinds.matrix <- matrix(NA, ncol = 100, nrow = 500)
for(i in 1:500){
  promoinds.matrix[i,] <- list.outputs[[i]][[2]][101:200]
}

base.errors.03.nonpromo <- base.errors.03[promoinds.matrix == 0]
base.errors.03.promo <- base.errors.03[promoinds.matrix == 1]

hist(base.errors.03.nonpromo)
hist(base.errors.03.promo)

sd(base.errors.03) ; sd(base.errors.03.promo) ; sd(base.errors.03.nonpromo)
mean(base.errors.03) ; mean(base.errors.03.promo) ; mean(base.errors.03.nonpromo)

#### Lost demand distributions ####

lostdemand.03 <- matrix(NA, ncol = 100, nrow = 500)
for(i in 1:500){
  lostdemand.03[i,] <- list.outputs[[i]][[5]][101:200]
}

promoinds.matrix <- matrix(NA, ncol = 100, nrow = 500)
for(i in 1:500){
  promoinds.matrix[i,] <- list.outputs[[i]][[2]][101:200]
}

relativepromoinds.matrix <- matrix(NA, ncol = 100, nrow = 500)
for(i in 1:500){
  relativepromoinds.matrix[i,] <- rep(c(1:9,0),10)
}

lostdemand.03.nonpromo <- lostdemand.03[promoinds.matrix == 0]
lostdemand.03.promo <- lostdemand.03[promoinds.matrix == 1]

coverages.03 <- 1 - mean(lostdemand.03 < 0)
coverages.03.promo <- 1 - mean(lostdemand.03.promo < 0)
coverages.03.nonpromo <- 1 - mean(lostdemand.03.nonpromo < 0)

coverages.03 ; coverages.03.promo ; coverages.03.nonpromo

# 

lostdemand.03.0 <- lostdemand.03[relativepromoinds.matrix == 0]
lostdemand.03.1 <- lostdemand.03[relativepromoinds.matrix == 1]
lostdemand.03.2 <- lostdemand.03[relativepromoinds.matrix == 2]
lostdemand.03.3 <- lostdemand.03[relativepromoinds.matrix == 3]
lostdemand.03.4 <- lostdemand.03[relativepromoinds.matrix == 4]
lostdemand.03.5 <- lostdemand.03[relativepromoinds.matrix == 5]
lostdemand.03.6 <- lostdemand.03[relativepromoinds.matrix == 6]
lostdemand.03.7 <- lostdemand.03[relativepromoinds.matrix == 7]
lostdemand.03.8 <- lostdemand.03[relativepromoinds.matrix == 8]
lostdemand.03.9 <- lostdemand.03[relativepromoinds.matrix == 9]

ld.03.0 <- 1 - mean(lostdemand.03.0 < 0)
ld.03.1 <- 1 - mean(lostdemand.03.1 < 0)
ld.03.2 <- 1 - mean(lostdemand.03.2 < 0)
ld.03.3 <- 1 - mean(lostdemand.03.3 < 0)
ld.03.4 <- 1 - mean(lostdemand.03.4 < 0)
ld.03.5 <- 1 - mean(lostdemand.03.5 < 0)
ld.03.6 <- 1 - mean(lostdemand.03.6 < 0)
ld.03.7 <- 1 - mean(lostdemand.03.7 < 0)
ld.03.8 <- 1 - mean(lostdemand.03.8 < 0)
ld.03.9 <- 1 - mean(lostdemand.03.9 < 0)

ld.03.0 ; ld.03.1 ; ld.03.2 ; ld.03.3 ; ld.03.4 ;ld.03.5 ; ld.03.6 ; ld.03.7
ld.03.8 ; ld.03.9

#barplot(height = c(ld.03.0,ld.03.1,ld.03.2,ld.03.3,ld.03.4,ld.03.5,ld.03.6,ld.03.7,
                   #ld.03.8,ld.03.9), names = c("0","1","2","3","4","5","6","7","8","9"),
                  #xlab = "Periods after most recent promo.",ylab = "CSL (achieved)",
                  #main = "CSL achieved by period relative to recent promo.")

l.1 <- c(ld.03.0,ld.03.1,ld.03.2,ld.03.3,ld.03.4,ld.03.5,ld.03.6,ld.03.7,
  ld.03.8,ld.03.9)

xseq <- 0:9

plot(x = xseq, y = l.01,type="o",ylim = c(0.55,0.95),pch = 16, col = 1,
     xlab = "Days since recent promotion",ylab = "CSL (achieved)",
     main = "Achieved CSL by period relative to recent promo.")
lines(x = xseq, y = l.03, type = "o", pch = 17, col = 2)
lines(x = xseq, y = l.05, type = "o", pch = 18, col = 3)
lines(x = xseq, y = l.1, type = "o", pch = 21, col = 4)
abline(h = 0.75, col = 1, lty = 2)

legend("bottomright",legend = c("sigma = 0.1","sigma = 0.3","sigma = 0.5","sigma = 1"),pch = c(16,17,18,21), col = 1:4)

#### Forecast function plots ####
