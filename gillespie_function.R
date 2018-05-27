#x1 is mRNA which has low amount, 10
#x2 is protein which is abundant, 1000
##### Initialisation #####
library(parallel)
set.seed(2018)

x1 <- 10
x2 <- 1000
iteration <- 1e6

beta1 <- 1
beta2 <- 100
lamda1 <- 100
lamda2 <- 100

createParameters <- function(){
  lambda1 <- seq(10,100,by=10)
  lambda2 <- 110 - lambda1 
  beta1 <- seq(10,100,by=10)
  beta2 <- 110 - beta1
  parameter_matrix <- matrix(nrow = 100, ncol = 4,
                     dimnames = list(1:100, c("lambda1", "lambda2", "beta1", "beta2")))
  parameter_matrix[,"lambda1"] <- lambda1
  parameter_matrix[,"lambda2"] <- lambda2
  counter <- 1
  for (row in seq(1,nrow(parameter_matrix), by =10)){
    parameter_matrix[row:(row+9),"beta1"] <- beta1[counter]
    parameter_matrix[row:(row+9),"beta2"] <- beta2[counter]
    counter <- counter + 1
  }
  return(as.data.frame(parameter_matrix))
}

parameters <- createParameters()

#auto regulation, suppose reduction of protein could lead to increase of mRNA production rate
gillespie <- function(x1, x2, iteration, beta1, beta2, lamda1, lamda2){
    
    time_keeper <- c()
    x1_storage <- x2_storage <- rep(0, iteration)
    r_plus_1 <- r_minus_1 <- r_plus_2 <- r_minus_2 <- rep(0, iteration)
    
    time_keeper[1] <- 0
    x1_storage[1] <- x1
    x2_storage[1] <- x2
    r_plus_1[1] <- lamda1
    r_minus_1[1] <- x1 * beta1
    r_plus_2[1] <- x2 * beta2
    r_minus_2[1] <- x1 * lamda2
    
    epoch_size <- round(iteration/10)
    epoch_test_1 <- epoch_test_2 <- rep(0, 10)

    #calculate rate of change in x1 and x2, both birth and death
    for (i in 2:iteration){
        x1_birth <- lamda1 
        x1_death <- x1 * beta1
        x2_death <- x2 * beta2
        x2_birth <- x1 * lamda2
    
        #total rate of change of x_1 and x_2
        Tot_rate <- sum(x1_birth, x1_death, x2_birth, x2_death)
        time_keeper[i] <- rexp(1,Tot_rate) + time_keeper[i-1]
    
        #choose which event takes place
        #based on which event takes place, update x1 or x2 accordingly
        u_event <- runif(1,0,1)
        
        stick_frac <- c(0, x1_birth, x1_death, x2_birth, x2_death)
        stick <- cumsum(stick_frac)/Tot_rate
        u_event <- runif(1)
        which_event <- tail(which(u_event > stick), 1) #Choose the falling region
        
        if (which_event == 1){ #x1_birth
          x1 <- x1 + 1
        } 
        if (which_event == 2) { #x1_death
          x1 <- x1 - 1
        }
        if (which_event == 3) { #x2_birth
          x2 <- x2 + 1
        }
        if (which_event == 4) { #x2_death
          x2 <- x2 - 1
        }
        
        #store new x1 and x2 values
        x1_storage[i] <- x1
        x2_storage[i] <- x2
        
        # store new postive and negative fluxes for x1 and x2
        r_plus_1[i] <- x1_birth
        r_minus_1[i] <- x1_death
        r_plus_2[i] <- x2_birth
        r_minus_2[i] <- x2_death
        if (i %% epoch_size == 0){
          mult <- i%/% epoch_size
          epoch_test_1[i] <- abs((lamda1/mean(x1_storage[((mult-1)*epoch_size+1):(mult*epoch_size)]))-beta1)/beta1
          epoch_test_2[i] <- abs(mean(x1_storage[((mult-1)*epoch_size+1):(mult*epoch_size)])/
            mean(x2_storage[((mult-1)*epoch_size+1):(mult*epoch_size)])*lamda2-beta2)/beta2
        }
    }
    flux <- data.frame(cbind(r_plus_1, r_minus_1, r_plus_2, r_minus_2))
    names(flux) <- c("R_plus_x1", "R_minus_x1", "R_plus_x2", "R_minus_x2")
    epochs <- data.frame(cbind(epoch_test_1, epoch_test_2))
    names(epochs) <- c("x1", "x2")
    return(list("parm" = c(beta1=beta1, beta2=beta2, lamda1=lamda1, lamda2=lamda2), 
                "x1" = x1_storage,"x2" = x2_storage, "time" = time_keeper, "flux" = flux,
                "epochs" = epcohs))
}

runSimulation <- function(parameters){
  cl <- makeCluster( 4 )
  clusterExport(cl, c("x1", "x2", "iteration", "gillespie", "parameters"))
  rows <- 1:100
  
  results <- parLapply(cl = cl,  X = rows, fun =  function(X){
    gillespie(x1, x2, iteration, parameters[X,1],parameters[X,2],parameters[X,3],parameters[X,4])
  })
  stopCluster(cl = cl)
  return(results)
}


simulation_results <- runSimulation(parameters)

check_simulation <- function(parameters){
  results <- runSimulation(parameters)
  relative_error_1 <- relative_error_2 <- rep(0, 100)
  check_stationary <- rep(0, 100)
  for (i in 1:100){
    result <- results[[i]]
    relative_error_1[i] <- result$epochs$x1[10]
    relative_error_2[i] <- result$epochs$x2[10]
    r_plus_x1 <- mean(result$flux$R_plus_x1)
    r_minus_x1 <- mean(result$flux$R_minus_x1)
    r_plus_x2 <- mean(result$flux$R_plus_x2)
    r_minus_x2 <- mean(result$flux$R_minus_x2)
    if ((all.equal(r_plus_x1, r_minus_x1)==TRUE) && all.equal(r_plus_x2, r_minus_x2)==TRUE){
      check_stationary[i] <- TRUE
    }
  }
  par(mfrow=c(1, 2))
  plot(hist(relative_error_1), main = "Histogram of relative error in flux balance relation
       for x1")
  plot(hist(relative_error_2), main = "Histogram of relative error in flux balance relation
       for x2")
  par(mfrow=c(1, 1))
  if (all(check_stationary==1)){
    cat(print("All simulations have ran long enough to describe stationarity"))
  }
  else{
    message("Simulations that have not ran long enough to decribe stationarity: ", 
            which(check_stationary==0))
  }
}


check_accuracy <- function(parameters){
  results <- runSimulation(parameters)
  eta11 <- eta12 <- eta22 <- rep(0, 100)
  eta11_analytic <- eta12_analytic <- eta22_analytic <- rep(0, 100)
  for (i in 1:100){
    result <- results[[i]]
    x1 <- result$x1
    x2 <- result$x2
    tau1 <- 1/result$parm["beta1"]
    tau2 <- 1/result$parm["beta2"]
    var_x1 <- var(x1)
    var_x2 <- var(x2)
    cov_x1_x2 <- cov(x1, x2)
    x1_mean <- mean(x1)
    x2_mean <- mean(x2)
    eta11[i] <- var_x1/(x1_mean^2)
    eta22[i] <- var_x2/(x2_mean^2)
    eta12[i] <- cov_x1_x2/(x1_mean*x2_mean)
    eta11_analytic[i] <- 1/x1_mean
    eta12_analytic[i] <- tau1/((tau1+tau2)*x1_mean)
    eta22_analytic[i] <- 1/x2_mean + tau1/((tau1+tau2)*x1_mean)
  }
  dev11 <- abs(eta11 - eta11_analytic)/eta11_analytic
  dev12 <- abs(eta12 - eta12_analytic)/eta12_analytic
  dev22 <- abs(eta22 - eta22_analytic)/eta22_analytic
  par(mfrow=c(1, 3))
  plot(hist(dev11), main="relative deviation of eta11")
  plot(hist(dev12), main="relative deviation of eta12")
  plot(hist(dev22), main="relative deviation of eta22")
}




noise_calculator <- function(dat1, dat2, parm){
    tau1 <- 1/parm["beta1"]
    tau2 <- 1/parm["beta2"]
    mRNA_noise <- 1/mean(dat1)
    ex_noise <- mRNA_noise + tau1/(tau1+tau2)
    in_noise <- 1/mean(dat2)
    return(list("Extrinsic" = ex_noise, "Intrinsic" = in_noise))
}

results <- gillespie(x1, x2, iteration, beta1, parameters$beta1, parameters$lambda1, parameters$lambda2)













par(mfrow=c(1,1),mar=c(4,3,1,1))
plot(results$time_keeper[900000:1000000], results$x2_storage[900000:1000000],type="l")
par(new=FALSE)
plot(results$time_keeper[900000:1000000], results$x1_storage[900000:1000000],type="l", col = "red")

##start with small beta1 and large beta2

