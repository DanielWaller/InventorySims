# Inventory simulations

#### 1: Implementing the (R,S) method on a simple (baseline + error) sales model ####
# R is the review period, we have it as 1 here

#(a) Specific DGP (we call it model 1)

DGP_1 <- function(baseline, sigma, Length){
  
  errors <- rnorm(mean = 0, sd = sigma, n = Length)
  dataseries <- errors + baseline ; return(list(1,dataseries))
  
}

#(b) Model estimation function (general)

get.parameter.estimates <- function(datalist,estL){
  
  if(datalist[[1]] == 1){
    N <- length(datalist[[2]])
    dataseries <- datalist[[2]][(N - estL + 1):N]
    baseline.est <- mean(dataseries) ; sigma.est <- sd(dataseries)
    parameter.estimates <- c(baseline.est,sigma.est)
    return(list(1,datalist[[2]],parameter.estimates,estL))
    
  }
  
}

#(b.II) Calculating S

calculate.S <- function(mu,sigma,R,L,fill.rate){
  PI <- R+L ; sigma.RL <- sqrt(PI * sigma^2) ; xhat.RL <- PI * mu
  ESPRC <- (1 - fill.rate) * xhat.RL ; Guk <- ESPRC/sigma.RL
  
  z <- sqrt(log(25/(Guk^2)))
  
  k <- (a[1] + (a[2] * z) + (a[3] * z^2) + (a[4] * z^3)) / (b[1] + (b[2] * z) + (b[3] * z^2) + (b[4] * z^3) + (b[5] * z^4))
  
  SS <- k * sigma.RL ; S <- xhat.RL + SS
  
  return.vec <- c(S,xhat.RL,sigma.RL)
  return(return.vec)
}

#(c) Initialising the on-hand stock and inventory position

# values for calculation of k
a <- c(-5.3925569,5.6211054,-3.8836830,1.0897299)
b <- c(1,-0.72496485,0.507326622,0.0669136868,-0.00329129114)

initialise.inventory.sim <- function(datalist,review.period,lead.time,fill.rate){
  
  R <- review.period ; L <- lead.time 
  if(datalist[[1]] == 1){
    
    mu <- datalist[[3]][1] ; sigma <- datalist[[3]][2]
    S <- calculate.S(mu,sigma,R,L,fill.rate)
    
  }
  # Initialise IP and on-hand-stock
  
  order.size <- xhat.RL/2 ; num.orders <- L ; periods.next <- 1# Assumption: R = 1
  Initial.IP <- S[1] ; Initial.OHS <- S[1] - (order.size * num.orders)
  Time.since.R <- 0
  
  orders <- list() ; orders[[1]] <- order.size ; order.times <- c(1)
  
  return(list(1,datalist[[2]],datalist[[3]],estL,review.period,lead.time,fill.rate,orders,order.times,periods.next,Initial.IP,Initial.OHS,Time.since.R))
}


#(d) Loop the inventory replenishment process

simulate.inventory.process <- function(datachunk){
  
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
    
    for(i in 1:N){
      
      # Observe demand
      IP <- max(IP - data.holdout[i] , 0) ; lostdemand[i] <- min(OHS - data.holdout[i],0) ; OHS <- max(OHS - data.holdout[i], 0)
      OHS.As[i] <- OHS
      # Receive any orders
      order.times <- order.times - 1
      if(order.times[1] == 0){
        OHS <- OHS + orders[[1]] ; orders[[1]] <- NA
        orders[[1]] <- NULL ; order.times <- order.times[-1]
        OHS.Bs[i] <- OHS
      }
      
      # Consider a replenishment
      Time.since.R <- Time.since.R + 1
      if(R == Time.since.R){
        new.param.ests <- get.parameter.estimates(datalist = list(1,data[1:(estL + i)]),estL = estL)
        mu <- new.param.ests[[3]][1] ; sigma <- new.param.ests[[3]][2] ; mus[i] <- mu ; sigmas[i] <- sigma
        newS <- calculate.S(mu,sigma,R,L,fill.rate)
        S <- newS[1] ; xhat.RL <- newS[2] ; sigma.RL <- newS[3] ; Ss[i] <- S
        order.size <- S - IP ; ordersizes[i] <- order.size
        t <- length(order.times) ; orders[[(t+1)]] <- order.size ; order.times[(t+1)] <- L
        IP <- IP + order.size ; IPs[i] <- IP
        Time.since.R <- 0
      }
      
    }

  }
  
  output <- list(lostdemand,sigmas,mus,IPs,OHS.As,OHS.Bs,Ss,ordersizes)
  return(output)
}

#### 1b. Running the simulation ####

set.seed(4052)

coverages <- numeric(500)
Saverage <- numeric(500)

for(i in 1:500){

  baseline = 60 ; sigma = 20 ; Length = 52 ; estL = 20
  R = 1 ; L = 1 ; fill.rate = 0.75

  data <- DGP_1(baseline,sigma,Length)
  datalist <- get.parameter.estimates(datalist = data,estL)
  datachunk <- initialise.inventory.sim(datalist,R,L,fill.rate)
  process <- simulate.inventory.process(datachunk)
  
  coverages[i] <- mean(process[[1]])/60
  Saverage[i] <- mean(process[[7]])
  print(i)
  
}





