#x1 is mRNA which has low amount, 10
#x2 is protein which is abundant, 1000
##### Initialisation #####
library(parallel)
set.seed(2018)

#This function creates combination of 500 different parameters
createParameters <- function(l1_min, l1_max, l2_min, l2_max,
                             b1_min,b1_max, b2_min, b2_max, number=40){
    lambda1 <- seq(l1_min,l1_max,by=2)
    lambda2 <- seq(l2_min,l2_max,by=2)
    beta1 <- seq(b1_min,b1_max,by=2)
    beta2 <- seq(b2_min,b2_max,by=2)
    check_dup <- c(1,1)
    while(any(table(check_dup)>1)){
        parameter_matrix <- sapply(1:number,function(x){
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



###################
# Calculate noise #
###################
noise_calculator <- function(x1_average, x2_average, parm){
    tau1 <- 1/parm["beta1"]
    tau2 <- 1/parm["beta2"]
    mRNA_noise <- 1/x1_average
    ex_noise <- mRNA_noise * tau1/(tau1+tau2)
    in_noise <- 1/x2_average
    names(ex_noise) <- NULL
    return(c("Extrinsic" = ex_noise, "Intrinsic" = in_noise))
}


getMean <- function(params){
    ##<x1> is always 5, based on parameters
    mean_x1 <- unname(unlist(params[1]/params[3]))
    mean_x2 <- unname(unlist((mean_x1*params[2])/params[4]))
    return(floor(c("x1"= mean_x1, "x2" = mean_x2)))
}

#################
# Core function #
#################
gillespie <- function(x1, x2, iteration, lambda1, beta1,
                      lambda2, beta2, check_interval=100000, ar_func=function(x){1}){
    # lambda1 <- parameters1[[1]]$lambda1
    # beta1 <- parameters1[[1]]$beta1
    # lambda2 <- parameters1[[1]]$lambda2
    # beta2 <- parameters1[[1]]$beta2
    # iteration <-5e4
    # x1 <- 100
    # x2 <- 100
    # check_interval <- 10000
    # 
    ## ar_func: autorepression function of x2, the birth rate of x1
    time_keeper <- x1_storage <- x2_storage <- rep(0, iteration)
    time_keeper[1] <- 0
    x1_storage[1] <- x1
    x2_storage[1] <- x2
    indi <- FALSE
    stationary_track <- c()
    #ar_func <- function(x)1
    #calculate rate of change in x1 and x2, both birth and death
    for (i in 2:iteration){
        x1_birth <- lambda1 * ar_func(x2)
        x1_death <- x1 * beta1
        x2_birth <- x1 * lambda2
        x2_death <- x2 * beta2
        
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
        
        #Store values
        x1_storage[i] <- x1
        x2_storage[i] <- x2
        
        #Explicitly check simulation
        
        if(i %% check_interval == 0){
            check_range <- (i-check_interval+1):i
            time_duration <- diff(time_keeper[c(check_range,i)])
            x1_average <- weighted.mean(x=x1_storage[check_range], w=time_duration)
            x2_average <- weighted.mean(x=x2_storage[check_range], w=time_duration)
            R_plus_1 <- lambda1 * ar_func(x2_average)
            R_minus_1 <- beta1 * x1_average
            R_plus_2 <- lambda2 * x1_average
            R_minus_2 <- beta2 * x2_average
            
            relative_R1_diff <- abs(R_plus_1-R_minus_1)/mean(c(R_plus_1, R_minus_1))
            relative_R2_diff <- abs(R_plus_2-R_minus_2)/mean(c(R_plus_2, R_minus_2))
            indi <- (relative_R1_diff <0.01 && relative_R2_diff<0.01)
            stationary_track <- c(stationary_track, indi)
        }
        
    }
    
    
    #Output result
    result <- cbind.data.frame(time = time_keeper, x1 = x1_storage, x2 = x2_storage)
    
    print("done")
    return(list("parm" = c(beta1=beta1, beta2=beta2, lambda1=lambda1, lambda2=lambda2),
                result = result,
                "check_interval" = check_interval,
                "stationary_reached" = stationary_track, 
                "Relative_diff"=c(relative_R1_diff, relative_R2_diff),
                "R" = c(R_plus_1, R_minus_1, R_plus_2, R_minus_2),
                "means" = c(x1_average, x2_average)))
    gc()
}

#####################
# create parameters #
#####################
parameters1 <- createParameters(50,200,50,500,10,80,10,100, number=75)
parameters2 <- createParameters(3000,4000,150,200,800,1200,10,30, number=25)
parameters <- c(parameters1, parameters2)
names(parameters) <- 1:100

##############
# Simulation #
##############
system.time(all_parm_result <- mclapply(parameters,function(x){
    gillespie(100,100, 5e5, x$lambda1, x$beta1, x$lambda2, x$beta2, 500000)}
    ,mc.cores=6L))


###############################
# Calculate noise and plot 2A #
###############################
noise_calculator <- function(dat1, dat2, parm){
    tau1 <- 1/parm["beta1"]
    tau2 <- 1/parm["beta2"]
    mRNA_noise <- 1/mean(dat1)
    ex_noise <- mRNA_noise + tau1/(tau1+tau2)
    in_noise <- 1/mean(dat2)
    return(list("Extrinsic" = ex_noise, "Intrinsic" = in_noise))
}


all_noise <- lapply(all_parm_result,function(x){
    check_data <- x$means
    noise_calculator(check_data[1],check_data[2],x$parm)
})

all_noise <- do.call(rbind.data.frame,all_noise)
names(all_noise) <- c("Extrinsic", "Intrinsic")
plot_col <- rep("black",length(all_noise$Extrinsic))
plot_col[all_noise$Extrinsic/all_noise$Intrinsic > 1.5] <- "blue"
plot_col[all_noise$Intrinsic/all_noise$Extrinsic > 1.5] <- "green"

coverged_loc <- sapply(all_parm_result, function(x)x$stationary_reached)
plot_col[coverged_loc] <- "red"

pdf("/home/skllr-b/Desktop/Temp/sba1/noise_plot_2A.pdf",width=7,height=7)
plot(log10(all_noise),col=plot_col,
     ylab="log10(Intrinsic noise)",
     xlab="log10(Extrinsic noise)",pch=19,cex=0.5)
lines(c(-1000,1000),c(-1000,1000))
legend("topleft", legend=c("Dominated by intrinsic noise",
                           "Dominated by extrinsic noise",
                           "In between",
                           "Reached stationary state"),
       col=c("green","blue","black","red"),pch=19)
dev.off()

############################################
############## 2B: flux error ##############
############################################
R_diff <- do.call(rbind,lapply(all_parm_result,function(x)x$Relative_diff))
pdf("/home/skllr-b/Desktop/Temp/sba1/flux_error_2B.pdf",width=10,height=5)
par(mfrow=c(1,2),mar=c(4,5,1,1))
hist(log10(R_diff[,1]),breaks=20,
     xlab="log10(Relative error in flux balance relation for x1)",
     main="",col="grey80"); box()
abline(v=log10(0.01),col="red")
hist(log10(R_diff[,2]),breaks=20,
     xlab="log10(Relative error in flux balance relation for x2)",
     main="",col="grey80"); box()
abline(v=log10(0.01),col="red")
dev.off()

#############################################
############## 2B: co-variance ##############
#############################################

##checks how do the simulated and theoretical etas differ. Also plots out
##simulated vs theoretical etas.
check_accuracy <- function(simulation_results_raw){
    dat_length <- length(simulation_results_raw)
    eta11 <- eta12 <- eta22 <- rep(0, dat_length)
    eta11_analytic <- eta12_analytic <- eta22_analytic <- rep(0, dat_length)
    
    #all_parm_result
    for (i in 1:100){
        simulation_results <- simulation_results_raw[[i]]
        result <- simulation_results$result
        x1 <- result$x1
        x2 <- result$x2
        tau1 <- 1/simulation_results$parm["beta1"]
        tau2 <- 1/simulation_results$parm["beta2"]
        var_x1 <- var(tail(x1,100000))
        var_x2 <- var(tail(x2,100000))
        cov_x1_x2 <- cov(tail(x1, 100000), tail(x2, 100000))
        x1_mean <- simulation_results$means[1]
        x2_mean <- simulation_results$means[2]
        eta11[i] <- var_x1/(x1_mean^2)
        eta22[i] <- var_x2/(x2_mean^2)
        eta12[i] <- cov_x1_x2/(x1_mean*x2_mean)
        theo_x1 <- simulation_results$parm["lambda1"]/simulation_results$parm["beta1"]
        theo_x2 <- simulation_results$parm["lambda2"]*theo_x1/simulation_results$parm["beta2"]
        eta11_analytic[i] <- 1/theo_x1
        eta12_analytic[i] <- tau1/((tau1+tau2)*theo_x1)
        eta22_analytic[i] <- 1/theo_x2 + tau1/((tau1+tau2)*theo_x1)
    }
    dev11 <- abs(eta11 - eta11_analytic)/eta11_analytic
    dev12 <- abs(eta12 - eta12_analytic)/eta12_analytic
    dev22 <- abs(eta22 - eta22_analytic)/eta22_analytic
    dev_result <- cbind.data.frame(dev11,dev12,dev22)
    all_eta <- cbind.data.frame(eta11,eta12,eta22)
    list(all_eta, dev_result)
}


dev_result <- check_accuracy(all_parm_result)[[2]]
pdf("/home/skllr-b/Desktop/Temp/sba1/covariance_plot2B_2.pdf", width = 18, height = 6)
par(mfrow=c(1, 3))
hist(dev_result[,1],breaks=100, main="relative deviation of eta11",col="grey80",xlab="relative deviation");box()
hist(dev_result[,2],breaks=100, main="relative deviation of eta12",col="grey80",xlab="relative deviation");box()
hist(dev_result[,3],breaks=100, main="relative deviation of eta22",col="grey80",xlab="relative deviation");box()
dev.off()

#########################
########## 2.C ##########
#########################
pdf("/home/skllr-b/Desktop/Temp/sba1/noise_correlation_2C.pdf",height=6,width=6)
plot(eta22,eta22_analytic,col=plot_col,pch=19,
     ylab="theoretical noise",xlab="observed noise",asp=1)
lines(-1000:1000,-1000:1000)

legend("topleft", legend=c("Dominated by intrinsic noise",
                           "Dominated by extrinsic noise",
                           "In between",
                           "Reached stationary state"),
       col=c("green","blue","black","red"),pch=19)


############################
# Plot 2D ##################
############################
last_stable_dat <- tail(all_parm_result[[which(coverged_loc)[20]]]$result,5000)

last_stable_x1 <- last_stable_dat$x1
last_stable_x2 <- last_stable_dat$x2
last_stable_time <- last_stable_dat$time
pdf("/home/skllr-b/Desktop/Temp/sba1/example_plot1_2D.pdf",width=10,height=5)
par(mfrow=c(1,1),mar=c(4,5,4,5))
plot(last_stable_time, last_stable_x1,type="l",lwd=1,
     xlab="Time", ylab="The number of x1")
par(new=TRUE)
plot(last_stable_time, last_stable_x2,type="l",col="red",lwd=1,yaxt="n",xlab="",ylab="")
axis(4,col="red")
mtext("The number of x2", side=4, line=3, col="red")
dev.off()

pdf("/home/skllr-b/Desktop/Temp/sba1/example_plot_2D.pdf",width=8,height=5)
par(mfrow=c(2,1),mar=c(4,5,1,5))
all_dat <- all_parm_result[[which(coverged_loc)[20]]]$result
plot(all_dat[seq(1,nrow(all_dat),len=50000),]$time[-(1:1000)],
     all_dat[seq(1,nrow(all_dat),len=50000),]$x1[-(1:1000)],
     type="l",lwd=1, xlab="Time",ylab="The number of x1")

plot(all_dat[seq(1,nrow(all_dat),len=50000),]$time[-(1:1000)],
     all_dat[seq(1,nrow(all_dat),len=50000),]$x2[-(1:1000)],
     type="l",lwd=1, xlab="Time",ylab="The number of x2",col="red")

mtext("The number of x2", side=4, line=3, col="red")
legend("topleft", legend=c("x1","x2"),col=c("black","red"),lty=1)
dev.off()

#########################
######### 2E ############
#########################
parameters2E <- createParameters(900,1000,1000,2000,100,200,20,50, number=5)
x1_birth_fun <- function(x){(800)/(800+x)}
#x1_birth_fun <- function(x){1}
system.time(all_parm_result_2E <- mclapply(parameters2E,function(x){
    gillespie(100,100, 1e7, lambda1=x$lambda1, beta1=x$beta1, lambda2=x$lambda2, beta2=x$beta2, 200000, ar_func = x1_birth_fun)}
    ,mc.cores=5L))
noise_calculator(parameters2E)
mean(tail(all_parm_result_2E$`1`$result$x2,500000))
all_parm_result_2E$`1`$means
all_noise <- lapply(all_parm_result_2E,function(x){
    check_data <- x$means
    noise_calculator(check_data[1],check_data[2],x$parm)
})
do.call(rbind,all_noise)
check_accuracy(all_parm_result_2E)
all_parm_result_2E$`4`$means
lapply(stationary_reached2E,function(x))

tail(all_parm_result_2E$`3`$result)
