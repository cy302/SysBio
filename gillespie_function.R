#x1 is mRNA which has low amount, 10
#x2 is protein which is abundant, 1000
##### Initialisation #####
library(parallel)
set.seed(2018)

#initial x1 and x2
x1 <- 10
x2 <- 1000

#maximal number of steps
iteration <- 1e6


beta1 <- 1
beta2 <- 100
lambda1 <- 100
lambda2 <- 100

# lambda1 <- seq(1,20,by=1)
# lambda2 <- seq(0,500,by=50)[-1]
# beta1 <- seq(1,20,by=1)
# beta2 <- seq(2,30,by=2)

#This function creates combination of 500 different parameters
createParameters <- function(){
    lambda1 <- seq(2,500,by=2)
    lambda2 <- seq(2,500,by=2)
    beta1 <- seq(2,500,by=2)
    beta2 <- seq(2,500,by=2)
    check_dup <- c(1,1)
    while(any(table(check_dup)>1)){
        parameter_matrix <- sapply(1:500,function(x){
            c(sample(lambda1,1),sample(lambda2,1),sample(beta1,1),sample(beta2,1))
        })
        parameter_matrix <- t(parameter_matrix)
        check_dup <- apply(parameter_matrix,1,function(x)paste(x,collapse = "_"))
    }
    parameters <- as.data.frame(parameter_matrix)
    colnames(parameters) <- c("lambda1","lambda2","beta1","beta2")
    parameters_list <- split(parameters, seq(nrow(parameters)))
    parameters_list
}

# Core function
gillespie <- function(x1, x2, iteration, lambda1, beta1, lambda2, beta2, check_interval, ar_func=function(x){1}){
    
    ## ar_func: autorepression function of x2, the birth rate of x1
    time_keeper <- x1_storage <- x2_storage <- rep(0, iteration)
    r_plus_1 <- r_minus_1 <- r_plus_2 <- r_minus_2 <- rep(0, iteration)
    
    time_keeper[1] <- 0
    x1_storage[1] <- x1
    x2_storage[1] <- x2
    r_plus_1[1] <- lambda1
    r_minus_1[1] <- x1 * beta1
    r_plus_2[1] <- x1 * lambda2
    r_minus_2[1] <- x2 * beta2
    
    epoch_test_1 <- epoch_test_2 <- c()
    indi <- FALSE
    
    #calculate rate of change in x1 and x2, both birth and death
    for (i in 2:iteration){
        x1_birth <- lambda1 * ar_func(x2)
        x1_death <- x1 * beta1
        x2_death <- x2 * beta2
        x2_birth <- x1 * lambda2
        
        #total rate of change of x_1 and x_2
        Tot_rate <- sum(x1_birth, x1_death, x2_birth, x2_death)
        time_keeper[i] <- rexp(1,Tot_rate) + time_keeper[i-1]
        
        #Choose a event
        u_event <- runif(1,0,1)
        stick_frac <- c(0, x1_birth, x1_death, x2_birth, x2_death)
        stick <- cumsum(stick_frac)/Tot_rate
        u_event <- runif(1)
        which_event <- tail(which(u_event > stick), 1) #Choose the falling region
        
        #Execute an event
        if (which_event == 1)x1 <- x1 + 1  #x1_birth
        if (which_event == 2)x1 <- x1 - 1  #x1_death
        if (which_event == 3)x2 <- x2 + 1  #x2_birth
        if (which_event == 4)x2 <- x2 - 1  #x2_death
        
        #Check stationery state
        if (i %% check_interval == 0){
            check_range <- (i-check_interval+1):i
            last_x1 <- x1_storage[check_range]
            last_x2 <- x2_storage[check_range]
            epoch_test_1 <- c(epoch_test_1, abs((lambda1/mean(last_x1)-beta1)/beta1))
            epoch_test_2 <- c(epoch_test_2, abs(mean(last_x1)/mean(last_x2)*lambda2-beta2)/beta2)
            R_plus_1 <- mean(r_plus_1[check_range])
            R_plus_2 <- mean(r_plus_2[check_range])
            R_minus_1 <- mean(r_minus_1[check_range])
            R_minus_2 <- mean(r_minus_2[check_range])
            relative_R_err1 <- abs(R_plus_1-R_minus_1)/mean(c(R_plus_1, R_minus_1))
            relative_R_err2 <- abs(R_plus_2-R_minus_2)/mean(c(R_plus_2, R_minus_2))
            indi <- (relative_R_err1<0.01)&&(relative_R_err2<0.01)
        }
        
        #Store values
        x1_storage[i] <- x1
        x2_storage[i] <- x2
        r_plus_1[i] <- x1_birth
        r_minus_1[i] <- x1_death
        r_plus_2[i] <- x2_birth
        r_minus_2[i] <- x2_death
        
        if (indi){
            break
            message("The random walk has reached stationarity at iteration: ", i)
        }
    }
    
    
    #Output result
    result <- cbind(time = time_keeper, x1 = x1_storage, x2 = x2_storage,
                    R_plus_x1 = r_plus_1, R_minus_x1 = r_minus_1,
                    R_plus_x2 = r_plus_2, R_minus_x2 = r_minus_2)
    result <- head(result, i)
    epochs <- cbind(x1 = epoch_test_1, x2 = epoch_test_2)
    gc()
    print("done")
    return(list("parm" = c(beta1=beta1, beta2=beta2, lambda1=lambda1, lambda2=lambda2), 
                result = result, "epochs" = epochs, 
                "check_interval" = check_interval))
}

#create parameters
parameters <- createParameters()


#run simulation 500 times, with different parameters
system.time(all_parm_result <- mclapply(parameters,function(x){
    gillespie(10,1000, 5e5, x$lambda1, x$beta1, x$lambda2, x$beta2, 100000)}
    ,mc.cores=7L))

##check if all simulations ran until the stationary state was reached.
#also plots histograms of relative errors in flux balance relations
check_simulation <- function(simulation_results){

    relative_diff_x1 <- sapply(simulation_results, function(x){
        last_dat <- as.data.frame(tail(x$result,10000))
        abs(mean(last_dat$R_plus_x1) - mean(last_dat$R_minus_x1))/
            mean(c(mean(last_dat$R_plus_x1),mean(last_dat$R_minus_x1)))
    })
    
    relative_diff_x2 <- sapply(simulation_results, function(x){
        last_dat <- as.data.frame(tail(x$result,10000))
        abs(mean(last_dat$R_plus_x2) - mean(last_dat$R_minus_x2))/
            mean(c(mean(last_dat$R_plus_x2),mean(last_dat$R_minus_x2)))
    })
    
    
    reach_ss <- which(relative_diff_x1<0.01 & relative_diff_x2<0.01)
    
    pdf("plot2B_1.pdf", width = 8, height = 6)
    par(mfrow=c(1, 2))
    hist(relative_diff_x1, main = "Histogram of relative error in flux balance relation for x1")
    hist(relative_diff_x2, main = "Histogram of relative error in flux balance relation for x2")
    par(mfrow=c(1, 1))
    dev.off()
    reach_ss
    
}

##checks how do the simulated and theoretical etas differ. Also plots out
##simulated vs theoretical etas.
check_accuracy <- function(simulation_results){
    eta11 <- eta12 <- eta22 <- rep(0, 100)
    eta11_analytic <- eta12_analytic <- eta22_analytic <- rep(0, 100)
    for (i in 1:100){
        result <- simulation_results[[i]]
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
    pdf("plot2B_2.pdf", width = 8, height = 6)
    par(mfrow=c(1, 3))
    plot(hist(dev11), main="relative deviation of eta11")
    plot(hist(dev12), main="relative deviation of eta12")
    plot(hist(dev22), main="relative deviation of eta22")
    dev.off()
    eta <- data.frame(cbind(eta11, eta12, eta22))
    names(eta) <- c("eta11", "eta12", "eta22")
    eta_analytic <- data.frame(cbind(eta11_analytic, eta12_analytic, eta22_analytic))
    names(eta_analytic) <- c("eta11", "eta12", "eta22")
    return(list("observedEta" = eta, "theoreticalEta" = eta_analytic))
}


## 2.C
D <- check_accuracy(simulation_results)
observed_eta22 <- D$observedEta$eta22
theoretical_eta22 <- D$theoreticalEta$eta22

noise_plot <- function(observed, theoretical){
    pdf("2Cplot.pdf", width = 8, height = 6)
    plot(observed, theoretical, type="l", xlab="Theoretical noise", ylab="Observed noise", main="Plot of observed
         theoretical protein noise")
    dev.off()
}

noise_plot(observed_eta22, theoretical_eta22)

noise_calculator <- function(dat1, dat2, parm){
    tau1 <- 1/parm["beta1"]
    tau2 <- 1/parm["beta2"]
    mRNA_noise <- 1/mean(dat1)
    ex_noise <- mRNA_noise + tau1/(tau1+tau2)
    in_noise <- 1/mean(dat2)
    return(list("Extrinsic" = ex_noise, "Intrinsic" = in_noise))
}

results <- gillespie(x1, x2, iteration, beta1, parameters$beta1, parameters$lambda1, parameters$lambda2)


plot2D <- function(results, x_lim){
    pdf(width = 8, height = 6, "2D.pdf")
    #plot only last epoch
    indeces <- (length(results$x1)- 100000):length(results$x1) 
    
    par(mfrow = c(1,1), mar = c(5,5,5,5))
    plot(results$time[indeces], results$x2[indeces], col = "red", 
         xlim = x_lim, type = "l", ylab = "Protein molecules",
         xlab = "Arbitrary time units",main = "Protein and mRNA levels over time")
    par(new = TRUE)
    plot(results$time[indeces], results$x1[indeces], col = "blue",
         xlim = x_lim, type = "l",
         axes = FALSE, xlab = NA, ylab = NA)
    axis(side = 4)
    mtext(side = 4, line = 3, "mRNA molecules")
    dev.off()
}

plot2D(results)

par(mfrow=c(1,1),mar=c(4,3,1,1))
plot(results$time_keeper[900000:1000000], results$x2_storage[900000:1000000],type="l")
par(new=FALSE)
plot(results$time_keeper[900000:1000000], results$x1_storage[900000:1000000],type="l", col = "red")

##start with small beta1 and large beta2
