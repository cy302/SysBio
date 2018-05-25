#x1 is mRNA which has low amount, 10
#x2 is protein which is abundant, 1000
##### Initialisation #####
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
  return(parameter_matrix)
}

parameters <- createParameters()

#auto regulation, suppose reduction of protein could lead to increase of mRNA production rate
gillespie <- function(x1, x2, iteration, K, beta1, beta2, lamda1, lamda2){
    
    time_keeper <- c()
    x1_storage <- x2_storage <-c()
    
    time_keeper[1] <- 0
    x1_storage[1] <- x1
    x2_storage[1] <- x2

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
    }
    return(list("parm" = c(K=K, beta1=beta1, beta2=beta2, lamda1=lamda1, lamda2=lamda2), 
                "x1" = x1_storage,"x2" = x2_storage, "time" = time_keeper))
}


noise_calculator <- function(dat1, dat2, parm){
    tau1 <- 1/parm["beta1"]
    tau2 <- 1/parm["beta2"]
    mRNA_noise <- 1/mean(dat1)
    ex_noise <- mRNA_noise + tau1/(tau1+tau2)
    in_noise <- 1/mean(dat2)
}

results <- gillespie(x1, x2, iteration, K, beta1, beta2, lamda1, lamda2)
par(mfrow=c(1,1),mar=c(4,3,1,1))
plot(results$time_keeper[900000:1000000], results$x2_storage[900000:1000000],type="l")
par(new=FALSE)
plot(results$time_keeper[900000:1000000], results$x1_storage[900000:1000000],type="l", col = "red")

##start with small beta1 and large beta2

