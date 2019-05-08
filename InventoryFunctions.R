# Inventory simulations

#### 1: Implementing the (R,S) method on a simple (baseline + error) sales model ####
# R is the review period, we have it as 1 here

#### (a) Specific DGP (we call it model 1) ####

DGP_1 <- function(baseline, sigma, Length){
  
  errors <- rnorm(mean = 0, sd = sigma, n = Length)
  dataseries <- errors + baseline ; dataseries[dataseries < 0] <- 0
  return(list(1,dataseries))
  
}

#### (b) Model estimation function (general) ####

get.parameter.estimates <- function(datalist,estL){
  
  if(datalist[[1]] == 1){
    N <- length(datalist[[2]])
    dataseries <- datalist[[2]][(N - estL + 1):N]
    baseline.est <- mean(dataseries) ; sigma.est <- sd(dataseries)
    parameter.estimates <- c(baseline.est,sigma.est)
    return(list(1,datalist[[2]],parameter.estimates,estL))
    
  }
  
}

#### (b.II) Calculating S ####

#### Fill rate ####

calculate.S.FR <- function(mu,sigma,R,L,fill.rate){
  PI <- R+L ; sigma.RL <- sqrt(PI * sigma^2) ; xhat.RL <- PI * mu ; xhat.R <- mu
  ESPRC <- (1 - fill.rate) * xhat.R ; Guk <- ESPRC/sigma.RL
  
  z <- sqrt(log(25/(Guk^2)))
  
  k <- (a[1] + (a[2] * z) + (a[3] * z^2) + (a[4] * z^3)) / (b[1] + (b[2] * z) + (b[3] * z^2) + (b[4] * z^3) + (b[5] * z^4))
  
  k <- k #+ 0.2
  
  SS <- k * sigma.RL ; S <- xhat.RL + SS
  
  return.vec <- c(S,xhat.RL,sigma.RL,k)
  return(return.vec)
}

#### CSL ####

calculate.S.CSL <- function(mu,sigma,R,L,CSL){
  PI <- R+L ; sigma.RL <- sqrt(PI * sigma^2) ; xhat.RL <- PI * mu
  CSLprop <- CSL/100 ; k <- qnorm(CSLprop)
  SS <- k * sigma.RL ; S <- xhat.RL + SS
  return.vec <- c(S,xhat.RL,sigma.RL,k)
  return(return.vec)
}


#### Empirical forecast error quantiles


#### (c) Initialising the on-hand stock and inventory position ####

#### Fill rate ####

# values for calculation of k
a <- c(-5.3925569,5.6211054,-3.8836830,1.0897299)
b <- c(1,-0.72496485,0.507326622,0.0669136868,-0.00329129114)

initialise.inventory.sim.FR <- function(datalist,review.period,lead.time,fill.rate){
  
  R <- review.period ; L <- lead.time 
  if(datalist[[1]] == 1){
    
    mu <- datalist[[3]][1] ; sigma <- datalist[[3]][2]
    S <- calculate.S.FR(mu,sigma,R,L,fill.rate)
    
  }
  # Initialise IP and on-hand-stock
  
  xhat.RL <- S[2] ; sigma.RL <- S[3]
  
  order.size <- xhat.RL/(R+L) ; num.orders <- L ; periods.next <- 1# Assumption: R = 1
  Initial.IP <- S[1] ; Initial.OHS <- S[1] - (order.size * num.orders)
  Time.since.R <- 0 ; k <- S[4]
  
  orders <- numeric() ; order.times <- c(1:num.orders)
  for(j in 1:num.orders){
    orders[j] <- order.size ; 
  }
  return(list(1,datalist[[2]],datalist[[3]],estL,review.period,lead.time,fill.rate,orders,order.times,periods.next,Initial.IP,Initial.OHS,Time.since.R,k))
}

#### CSL ####

initialise.inventory.sim.CSL <- function(datalist,review.period,lead.time,CSL){
  
  R <- review.period ; L <- lead.time 
  if(datalist[[1]] == 1){
    
    mu <- datalist[[3]][1] ; sigma <- datalist[[3]][2]
    S <- calculate.S.CSL(mu,sigma,R,L,CSL)
    
  }
  # Initialise IP and on-hand-stock
  
  xhat.RL <- S[2] ; sigma.RL <- S[3]
  
  order.size <- xhat.RL/(R+L) ; num.orders <- L ; periods.next <- 1# Assumption: R = 1
  Initial.IP <- S[1] ; Initial.OHS <- S[1] - (order.size * num.orders)
  Time.since.R <- 0 ; k <- S[4]
  
  orders <- numeric() ; order.times <- c(1:num.orders)
  for(j in 1:num.orders){
    orders[j] <- order.size ; 
  }
  return(list(1,datalist[[2]],datalist[[3]],estL,review.period,lead.time,CSL,orders,order.times,periods.next,Initial.IP,Initial.OHS,Time.since.R,k))
}



#### (d) Loop the inventory replenishment process ####

#### Fill rate ####

simulate.inventory.process.FR <- function(datachunk){
  
  if(datachunk[[1]] == 1){
    data <- datachunk[[2]] ;  mu <- datachunk[[3]][1] ; sigma <- datachunk[[3]][2] ; estL <- datachunk[[4]] ; R <- datachunk[[5]]
    L <- datachunk[[6]] ; FR <- datachunk[[7]] ; orders <- datachunk[[8]] ; order.times <- datachunk[[9]]
    periods.next <- datachunk[[10]] ; IP <- datachunk[[11]] ; OHS <- datachunk[[12]] ; Time.since.R <- datachunk[[13]]
    
    N <- length(data) - estL
    data.holdout <- data[(estL + 1):length(data)]
    
    lostdemand <- numeric(N)
    mus <- numeric(N)
    sigmas <- numeric(N)
    IPs <- numeric(N)
    OHS.As <- numeric(N)
    OHS.Bs <- numeric(N)
    Ss <- numeric(N)
    ordersizes <- numeric(N)
    ks <- numeric(N)
    backorders <- 0
    backtrack <- numeric(N)
    back <- FALSE
    
    for(i in 1:N){
      
      # Observe demand
      data.fulfilled <- min(OHS, backorders + data.holdout[i])
      lostdemand[i] <- min(OHS - data.holdout[i] - backorders,0) ; OHS <- max(OHS - data.holdout[i] - backorders, 0)
      
      ## BACKORDERS SECTION ##
      
      if(back == TRUE){
      
        if(lostdemand[i] < 0){
          backorders <- abs(lostdemand[i])
        }
      
        if(lostdemand[i] == 0){
          backorders <- 0
        }
      
        backtrack[i] <- backorders
      }
      
      ##                    ##
      
      OHS.As[i] <- OHS
      # Receive any orders
      order.times <- order.times - 1
      if(order.times[1] == 0){
        OHS <- OHS + orders[1] ; orders <- orders[-1] ; order.times <- order.times[-1]
        OHS.Bs[i] <- OHS
      }
      
      # Consider a replenishment
      Time.since.R <- Time.since.R + 1
      if(R == Time.since.R){
        new.param.ests <- get.parameter.estimates(datalist = list(1,data[1:(estL + i)]),estL = estL)
        mu <- new.param.ests[[3]][1] ; sigma <- new.param.ests[[3]][2] ; mus[i] <- mu ; sigmas[i] <- sigma
        newS <- calculate.S.FR(mu,sigma,R,L,fill.rate)
        S <- newS[1] + backorders ; xhat.RL <- newS[2] ; sigma.RL <- newS[3] ; Ss[i] <- S ; ks[i] <- newS[4]
        IP <- OHS + sum(orders) - backorders ; IPs[i] <- IP
        order.size <- max(S - IP,0) ; ordersizes[i] <- order.size
        t <- length(order.times) ; orders[[(t+1)]] <- order.size ; order.times[(t+1)] <- L
        Time.since.R <- 0
      }
      
    }

  }
  
  output <- list(lostdemand,sigmas,mus,IPs,OHS.As,OHS.Bs,Ss,ordersizes,ks,backtrack)
  return(output)
}

#### CSL ####

simulate.inventory.process.CSL <- function(datachunk){
  
  if(datachunk[[1]] == 1){
    data <- datachunk[[2]] ;  mu <- datachunk[[3]][1] ; sigma <- datachunk[[3]][2] ; estL <- datachunk[[4]] ; R <- datachunk[[5]]
    L <- datachunk[[6]] ; CSL <- datachunk[[7]] ; orders <- datachunk[[8]] ; order.times <- datachunk[[9]]
    periods.next <- datachunk[[10]] ; IP <- datachunk[[11]] ; OHS <- datachunk[[12]] ; Time.since.R <- datachunk[[13]]
    
    N <- length(data) - estL
    data.holdout <- data[(estL + 1):length(data)]
    
    lostdemand <- numeric(N)
    mus <- numeric(N)
    sigmas <- numeric(N)
    IPs <- numeric(N)
    OHS.As <- numeric(N)
    OHS.Bs <- numeric(N)
    Ss <- numeric(N)
    ordersizes <- numeric(N)
    ks <- numeric(N)
    backorders <- 0
    backtrack <- numeric(N)
    back <- FALSE
    
    for(i in 1:N){
      
      # Observe demand
      data.fulfilled <- min(OHS, backorders + data.holdout[i])
      lostdemand[i] <- min(OHS - data.holdout[i] - backorders,0) ; OHS <- max(OHS - data.holdout[i] - backorders, 0)
      
      ## BACKORDERS SECTION ##
      
      if(back == TRUE){
        
        if(lostdemand[i] < 0){
          backorders <- abs(lostdemand[i])
        }
        
        if(lostdemand[i] == 0){
          backorders <- 0
        }
        
        backtrack[i] <- backorders
      }
      
      ##                    ##
      
      OHS.As[i] <- OHS
      # Receive any orders
      order.times <- order.times - 1
      if(order.times[1] == 0){
        OHS <- OHS + orders[1] ; orders <- orders[-1] ; order.times <- order.times[-1]
        OHS.Bs[i] <- OHS
      }
      
      # Consider a replenishment
      Time.since.R <- Time.since.R + 1
      if(R == Time.since.R){
        new.param.ests <- get.parameter.estimates(datalist = list(1,data[1:(estL + i)]),estL = estL)
        mu <- new.param.ests[[3]][1] ; sigma <- new.param.ests[[3]][2] ; mus[i] <- mu ; sigmas[i] <- sigma
        newS <- calculate.S.CSL(mu,sigma,R,L,CSL)
        S <- newS[1] + backorders ; xhat.RL <- newS[2] ; sigma.RL <- newS[3] ; Ss[i] <- S ; ks[i] <- newS[4]
        IP <- OHS + sum(orders) - backorders ; IPs[i] <- IP
        order.size <- max(S - IP,0) ; ordersizes[i] <- order.size
        t <- length(order.times) ; orders[[(t+1)]] <- order.size ; order.times[(t+1)] <- L
        Time.since.R <- 0
      }
      
    }
    
  }
  
  output <- list(lostdemand,sigmas,mus,IPs,OHS.As,OHS.Bs,Ss,ordersizes,ks,backtrack)
  return(output)
}

#### 2. Running the simulation ####

#### Fill Rate ####

set.seed(4052)

coverages <- numeric(500)
inventories <- numeric(500)
Saverage <- numeric(500)
ks.overall <- numeric(500)

for(i in 1:500){
    
  baseline = 60 ; sigma = 60 ; Length = 52 ; estL = 20
  R = 1 ; L = 1 ; fill.rate = 0.99
  data <- DGP_1(baseline,sigma,Length)
  datalist <- get.parameter.estimates(datalist = data,estL)
  datachunk <- initialise.inventory.sim.FR(datalist,R,L,fill.rate)
  process <- simulate.inventory.process.FR(datachunk)
  
  # Coverage
  
  coverages[i] <- mean(process[[1]])/60
  #coverages[i] <- mean(process[[1]][33:84])/600
  #coverages[i] <- mean(process[[1]][1033:1084])/600
  
  # Average inventory held
  
  ohsproc <- c(rbind(process[[5]],process[[6]])) ; ohsproc <- c(datachunk[[12]],ohsproc)
  ohs.av <- numeric(length(process[[1]]))
  for(j in 1:(length(process[[1]]))){
    ohs.av[j] <- mean(ohsproc[(2*j - 1):(2*j)])
  }
  inventories[i] <- mean(ohs.av)
  
  # Ks
  
  ks <- process[[9]]
  ks.overall[i] <- mean(c(datachunk[[14]],ks))
  print(i)
  
}

#### CSL ####

set.seed(4052)

coverages <- numeric(500)
inventories <- numeric(500)
Saverage <- numeric(500)
ks.overall <- numeric(500)
cslprops.overall <- numeric(500)
forecast.error.matrix <- matrix(NA,nrow = 500,ncol = 100)
stockout.matrix <- matrix(NA, nrow = 500, ncol = 100)

for(i in 1:500){
  
  baseline = 60 ; sigma = 40 ; Length = 122 ; estL = 20
  R = 1 ; L = 1 ; CSL = 95
  data <- DGP_1(baseline,sigma,Length)
  datalist <- get.parameter.estimates(datalist = data,estL)
  datachunk <- initialise.inventory.sim.CSL(datalist,R,L,CSL)
  process <- simulate.inventory.process.CSL(datachunk)
  
  # Coverage
  
  coverages[i] <- mean(process[[1]])/60
  #coverages[i] <- mean(process[[1]][33:84])/600
  #coverages[i] <- mean(process[[1]][1033:1084])/600
  
  # Average inventory held
  
  ohsproc <- c(rbind(process[[5]],process[[6]])) ; ohsproc <- c(datachunk[[12]],ohsproc)
  ohs.av <- numeric(length(process[[1]]))
  for(j in 1:(length(process[[1]]))){
    ohs.av[j] <- mean(ohsproc[(2*j - 1):(2*j)])
  }
  inventories[i] <- mean(ohs.av)
  
  # Ks
  
  ks <- process[[9]]
  ks.overall[i] <- mean(c(datachunk[[14]],ks))
  print(i)
  
  # CSL
  
  cslprops.overall[i] <- (sum(process[[1]] < 0))/(length(process[[1]]))
  
  # Forecast estimates mu against actual data observed
  # Mu estimated at end of period 1 -> forecast for the replenishment order that arrives in period 2 - satisfies demand in period 3
  
  forecasts <- process[[3]][1:(length(process[[3]]) - 2)] ; estL <- datachunk[[4]]
  observed <- datachunk[[2]][-c(1:(estL + 2))]
  forecast.errors <- forecasts - observed
  stockouts <- process[[1]][3:(length(process[[1]]))]
  
  forecast.error.matrix[i,] <- forecast.errors
  stockout.matrix[i,] <- stockouts
  
}

#### Illustration plot ####

N <- length(process[[1]])

ohsproc <- c(rbind(process[[5]],process[[6]])) ; ohsproc <- c(datachunk[[12]],ohsproc)
seq1 <- c(0,rep(1:N,each = 2))

ohs.av <- numeric(N)
for(i in 1:N){
  ohs.av[i] <- mean(ohsproc[(2*i - 1):(2*i)])
}

plot(x = seq1, y = ohsproc,pch = 16, 
     ylim = c(-100, 400),
     xlab = "Time period", ylab = "On-hand-stock",
     main = "Inventory process - Promotional case")
abline(h = 0, col = "red")
for(i in 1:280){
  lines(x = seq1[(2*i):(2*i + 1)], y = ohsproc[(2*i):(2*i + 1)],lty = 2)
}
for(i in 1:280){
  lines(x = seq1[(2*i - 1):(2*i)], y = ohsproc[(2*i - 1):(2*i)],lty = 1)
}
points(x = which(process[[1]] < 0), y = process[[1]][process[[1]] < 0 ], pch = 18, col = "blue")

pts1 <- which(process[[1]] < 0) ; pts2 <- process[[1]][process[[1]] < 0]

for(i in 1:(length(pts1))){
  lines(x = rep(pts1[i],2), y = c(0, pts2[i]),col = "blue")
}

legend("topleft",legend = c("OHS","Lost sales"),pch = c(16,18), col = c(1,4))

abline(v = c(20,40,60,80,100),lty = 2, col = "red")

# Cheat loop

for(i in 1:5){
  ohsproc[((i*40)+1)] <- ohsproc[((i*40)+1)] + 180
}

ohsproc[41] <- ohsproc[41] + 50

#### CASE 1 PLOTS ####

# Fill rate vs. sigma - use sigma1, sigma2, sigma3 as vectors for fillrate
# Use Sigma1, Sigma2, Sigma3 to store the average inventory carried

sigma1 <- numeric(6) ; sigma2 <- numeric(6) ; sigma3 <- numeric(6)
Sigma1 <- numeric(6) ; Sigma2 <- numeric(6) ; Sigma3 <- numeric(6)

sigma3[6] <- mean(coverages)
Sigma3[6] <- mean(inventories)

empfr1 <- 1+ sigma1 ; empfr2 <- 1 + sigma2 ; empfr3 <- 1 + sigma3

# Plot - empirical vs theor FR by sigma value

frseq <- c(0.75,0.8,0.85,0.9,0.95,0.99)
plot(x = frseq, y = empfr1, ylim = c(0.65,1), xlim = c(0.75, 1), pch = 15, cex = 1.25, lwd = 2, type = "o",col = "red",
     xlab = "Fill rate (theoretical)", ylab = "Fill rate (empirical)", main = "Fill rate (Emp. vs. Theor.) - by sigma value")
abline(a = 0, b = 1, col = "black",lty = 2)
lines(x = frseq, y = empfr2, pch = 16, cex=  1.25, lwd = 2, col = "green",type= "o")
lines(x = frseq, y = empfr3, pch = 17, cex=  1.25, lwd = 2, col = "blue",type= "o")

legend("bottomright", legend = c("sigma = 10","sigma = 40","sigma = 60"),pch = c(15,16,17),
       col = c(2,3,4))

# Plot - Emp. fill rate vs. average OHI

e1 <- 1 - empfr1 ; e2 <- 1 - empfr2 ; e3 <- 1 - empfr3
plot(y = e1, x = Sigma1,type= "o",ylim = c(0,0.35),xlim = c(0,180),ylab = "1 - fill rate",
     xlab = "Average on-hand inventory", main = "Fill-rate vs. average on-hand inventory- by sigma value",
     col = "red", pch = 15, lwd = 1.5, cex = 1.5)
lines(y = e2, x = Sigma2, type = "o",col = "green", pch = 16, lwd = 1.5, cex = 1.5)
lines(y = e3, x = Sigma3, type = "o",col = "blue", pch = 17, lwd = 1.5, cex = 1.5)

legend("topright", legend = c("sigma = 10","sigma = 40","sigma = 60"),pch = c(15,16,17),
       col = c(2,3,4),cex = 1.5)

#### CASE 1 PLOTS continued ####

# Plot history length against fill rate
# Use Hist1, Hist2, Hist 3 for history lengths. Sigma = 60
hist1 <- numeric(6) ; hist2 <- numeric(6) ; hist3 <- numeric(6)

hist3[6] <- mean(coverages)

empfr1 <- 1+ hist1 ; empfr2 <- 1 + hist2 ; empfr3 <- 1 + hist3

frseq <- c(0.75,0.8,0.85,0.9,0.95,0.99)
plot(x = frseq, y = empfr1, ylim = c(0.7,1), xlim = c(0.75, 1), pch = 15, cex = 1.5, lwd = 1.5, type = "o",col = "red",
     xlab = "Fill rate (theoretical)", ylab = "Fill rate (empirical)", main = "Fill rate (Emp. vs. Theor.) - by history length")
abline(a = 0, b = 1, col = "black",lty = 2)
lines(x = frseq, y = empfr2, pch = 16, cex=  1.25, lwd = 2, col = "green",type= "o")
lines(x = frseq, y = empfr3, pch = 17, cex=  1.25, lwd = 2, col = "blue",type= "o")

legend("bottomright", legend = c("T = 10","T = 20","T = 30"),pch = c(15,16,17),
       col = c(2,3,4),cex = 1.5)


#### Plot history length against lead time ####
# Use Hist1, Hist2, Hist 3 for history lengths. Sigma = 60
lead1 <- numeric(6) ; lead2 <- numeric(6) ; lead3 <- numeric(6)

lead3[6] <- mean(coverages)

empfr1 <- 1+ lead1 ; empfr2 <- 1 + lead2 ; empfr3 <- 1 + lead3

frseq <- c(0.75,0.8,0.85,0.9,0.95,0.99)
plot(x = frseq, y = empfr1, ylim = c(0.7,1), xlim = c(0.75, 1), pch = 15, cex = 1.5, lwd = 1.5, type = "o",col = "red",
     xlab = "Fill rate (theoretical)", ylab = "Fill rate (empirical)", main = "Fill rate (Emp. vs. Theor.) - by lead time")
abline(a = 0, b = 1, col = "black",lty = 2)
lines(x = frseq, y = empfr2, pch = 16, cex=  1.25, lwd = 2, col = "green",type= "o")
lines(x = frseq, y = empfr3, pch = 17, cex=  1.25, lwd = 2, col = "blue",type= "o")

legend("bottomright", legend = c("L = 1","L = 2","L = 3"),pch = c(15,16,17),
       col = c(2,3,4),cex = 1.5)
#### Manual k-calibration case ####
k.adj <- c(0,0.1,0.2,0.3,0.4,0.5)
#### Backorders case ####
# Use sigma1 vector from earlier plot
back.sigma1 <- numeric(6)

back.sigma1[6] <- mean(coverages)
empfr1 <- 1+ sigma1 ; empfr2 <- 1 + back.sigma1

frseq <- c(0.75,0.8,0.85,0.9,0.95,0.99)
plot(x = frseq, y = empfr1, ylim = c(0.7,1), xlim = c(0.75, 1), pch = 15, cex = 1.25, lwd = 2, type = "o",col = "red",
     xlab = "Fill rate (theoretical)", ylab = "Fill rate (empirical)", main = "Fill rate (Emp. vs. Theor.) - by sigma value")
abline(a = 0, b = 1, col = "black",lty = 2)
lines(x = frseq, y = empfr2, pch = 16, cex=  1.25, lwd = 2, col = "green",type= "o")

legend("bottomright", legend = c("Lost Sales","Back Orders"),pch = c(15,16),
       col = c(2,3))

#### CSL Plots ####

# CSLvals <- c(95,85,75) # sigma vals <- c(10,40,60)
CSL1 <- numeric(3) ; CSL2 <- numeric(3) ; CSL3 <- numeric(3)

CSL3[3] <- 100 - (100* mean(cslprops.overall))

plot(y = CSL1, x = c(10,40,60), ylim = c(70, 100), xlim = c(0,70),
     pch = 15, cex = 1.25, lwd = 2, type = "o",col = "red",
     main = "CSL - theoretical vs. empirical (varying sigma)",
     xlab = "Sigma",ylab = "CSL")
abline(h = 75, col = "black",lty = 2)
abline(h = 85, col = "black",lty = 2)
abline(h = 95, col = "black",lty = 2)
lines(x = c(10,40,60), y = CSL2, pch = 16, cex=  1.25, lwd = 2, col = "green",type= "o")
lines(x = c(10,40,60), y = CSL3, pch = 17, cex=  1.25, lwd = 2, col = "blue",type= "o")

legend("bottomright", legend = c("75%","85%","95%"),pch = c(15,16,17),
       col = c(2,3,4))

#### CSL vs. sigma - use sigma1, sigma2, sigma3 as vectors for CSL ####
# Use Sigma1, Sigma2, Sigma3 to store the average inventory carried

sigma1 <- numeric(6) ; sigma2 <- numeric(6) ; sigma3 <- numeric(6)
Sigma1 <- numeric(6) ; Sigma2 <- numeric(6) ; Sigma3 <- numeric(6)

sigma3[6] <- mean(cslprops.overall)
Sigma3[6] <- mean(inventories)

empcsl1 <- 1 - sigma1 ; empcsl2 <- 1 - sigma2 ; empcsl3 <- 1 - sigma3

# Plot - empirical vs theor FR by sigma value

frseq <- c(0.75,0.80,0.85,0.90,0.95,0.99)
plot(x = frseq, y = empcsl1, ylim = c(0.75,1.00), xlim = c(0.75, 1.00), pch = 15, cex = 1.25, lwd = 2, type = "o",col = "red",
     xlab = "CSL (theoretical)", ylab = "CSL (empirical)", main = "CSL rate (Emp. vs. Theor.) - by sigma value")
abline(a = 0, b = 1, col = "black",lty = 2)
lines(x = frseq, y = empcsl2, pch = 16, cex=  1.25, lwd = 2, col = "green",type= "o")
lines(x = frseq, y = empcsl3, pch = 17, cex=  1.25, lwd = 2, col = "blue",type= "o")

legend("bottomright", legend = c("sigma = 10","sigma = 40","sigma = 60"),pch = c(15,16,17),
       col = c(2,3,4))

# Plot - Emp. fill rate vs. average OHI

e1 <- 1 - empcsl1 ; e2 <- 1 - empcsl2 ; e3 <- 1 - empcsl3
plot(y = e1, x = Sigma1,type= "o",ylim = c(0,0.35),xlim = c(0,205),ylab = "1 - CSL",
     xlab = "Average on-hand inventory", main = "CSL vs. average on-hand inventory- by sigma value",
     col = "red", pch = 15, lwd = 1.5, cex = 1.5)
lines(y = e2, x = Sigma2, type = "o",col = "green", pch = 16, lwd = 1.5, cex = 1.5)
lines(y = e3, x = Sigma3, type = "o",col = "blue", pch = 17, lwd = 1.5, cex = 1.5)

legend("topright", legend = c("sigma = 10","sigma = 40","sigma = 60"),pch = c(15,16,17),
       col = c(2,3,4),cex = 1.5)

#### Plot history length against CSL ####
# Use Hist1, Hist2, Hist 3 for history lengths. Sigma = 60
hist1 <- numeric(6) ; hist2 <- numeric(6) ; hist3 <- numeric(6)

hist3[6] <- mean(cslprops.overall)

empcsl1 <- 1 - hist1 ; empcsl2 <- 1 - hist2 ; empcsl3 <- 1 - hist3

frseq <- c(0.75,0.8,0.85,0.9,0.95,0.99)
plot(x = frseq, y = empcsl1, ylim = c(0.7,1), xlim = c(0.75, 1), pch = 15, cex = 1.5, lwd = 1.5, type = "o",col = "red",
     xlab = "CSL (theoretical)", ylab = "CSL (empirical)", main = "CSL (Emp. vs. Theor.) - by history length")
abline(a = 0, b = 1, col = "black",lty = 2)
lines(x = frseq, y = empcsl2, pch = 16, cex=  1.25, lwd = 2, col = "green",type= "o")
lines(x = frseq, y = empcsl3, pch = 17, cex=  1.25, lwd = 2, col = "blue",type= "o")

legend("bottomright", legend = c("T = 10","T = 20","T = 30"),pch = c(15,16,17),
       col = c(2,3,4),cex = 1.5)

#### CSL plot Lead time ####

lead1 <- numeric(6) ; lead2 <- numeric(6) ; lead3 <- numeric(6)

lead1[1] <- mean(cslprops.overall)

empcsl1 <- 1 - lead1 ; empcsl2 <- 1 - lead2 ; empcsl3 <- 1 - lead3

frseq <- c(0.75,0.8,0.85,0.9,0.95,0.99)
plot(x = frseq, y = empcsl1, ylim = c(0.7,1), xlim = c(0.75, 1), pch = 15, cex = 1.5, lwd = 1.5, type = "o",col = "red",
     xlab = "CSL (theoretical)", ylab = "CSL (empirical)", main = "CSL (Emp. vs. Theor.) - by lead time")
abline(a = 0, b = 1, col = "black",lty = 2)
lines(x = frseq, y = empcsl2, pch = 16, cex=  1.25, lwd = 2, col = "green",type= "o")
lines(x = frseq, y = empcsl3, pch = 17, cex=  1.25, lwd = 2, col = "blue",type= "o")

legend("bottomright", legend = c("L = 1","L = 2","L = 3"),pch = c(15,16,17),
       col = c(2,3,4),cex = 1.5)
