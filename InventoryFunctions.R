# Inventory simulations

#### 1: Implementing the (R,S) method on a simple (baseline + error) sales model ####
# R is the review period, we have it as 1 here

#(a) Specific DGP (we call it model 1)

DGP_1 <- function(baseline, sigma, L){
  
  errors <- rnorm(mean = 0, sd = sigma, n = L)
  dataseries <- errors + baseline ; return(list(1,dataseries))
  
}

#(b) Model estimation function (general)

get.parameter.estimates <- function(datalist,estL){
  
  if(datalist[[1]] == 1){
    
    dataseries <- datalist[[2]][1:estL]
    baseline.est <- mean(dataseries) ; sigma.est <- sd(dataseries)
    parameter.estimates <- c(baseline.est,sigma.est)
    return(list(1,datalist[[2]],parameter.estimates,estL))
    
  }
  
}

#(c) Initialising the on-hand stock and inventory position

# values for calculation of k
a <- c(-5.3925569,5.6211054,-3.8836830,1.0897299)
b <- c(1,-0.72496485,0.507326622,0.0669136868,-0.00329129114)

initialise.inventory.sim <- function(datalist,review.period,lead.time,fill.rate){
  
  R <- review.period ; L <- lead.time 
  if(datalist[[1]] == 1){
    
    mu <- datalist[[3]][1] ; sigma <- datalist[[3]][2]
    
    # Variance of demand in interval R+L
    
    PI <- R+L ; sigma.RL <- sqrt(PI * sigma^2) ; xhat.RL <- PI * mu
    ESPRC <- (1 - fill.rate) * xhat.RL ; Guk <- ESPRC/sigma.RL
    
    z <- sqrt(log(25/(Guk^2)))
    
    k <- (a[1] + (a[2] * z) + (a[3] * z^2) + (a[4] * z^3)) / (b[1] + (b[2] * z) + (b[3] * z^2) + (b[4] * z^3) + (b[5] * z^4))
    
    SS <- k * sigma.RL ; S <- xhat.RL + SS
    
  }
  # Initialise IP and on-hand-stock
  
 order.size <- xhat.RL ; num.orders <- L - 1 ; periods.next <- 1# Assumption: R = 1
  Initial.IP <- S ; Initial.OHS <- S - (order.size * num.orders)
  Time.since.R <- 0
  
  orders <- list() ; orders[[1]] <- order.size ; order.times <- c(1)
}

#(d) Loop the inventory replenishment process






# some graphs that 




