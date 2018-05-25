#x1 is mRNA which has low amount, 10
#x2 is protein which is abundant, 1000
##### Initialisation #####
x1 <- 10
x2 <- 1000
iteration <- 1e6

beta1 <- 1
beta2 <- 1
lamda1 <- 100
lamda2 <- 100

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
        if (u_event < x1_birth/Tot_rate){ #x1_birth
          x1 <- x1 + 1
        } else if ((x1_birth/Tot_rate < u_event) &
                   (u_event < (x1_birth+x1_death)/Tot_rate)) { #x1_death
          x1 <- x1 - 1
        } else if (((x1_birth+x1_death)/Tot_rate < u_event) &
                  (u_event < (x1_birth+x1_death+x2_birth)/Tot_rate)) { #x2_birth
          x2 <- x2 + 1
        } else if (u_event > (x1_birth+x1_death+x2_birth)/Tot_rate) { #x2_death
          x2 <- x2 - 1
        }
        
        #store new x1 and x2 values
        x1_storage[i] <- x1
        x2_storage[i] <- x2
    }
    return(list("x1_storage" = x1_storage,"x2_storage" = x2_storage,
                "time_keeper" = time_keeper))
}
results <- gillespie(x1, x2, iteration, K, beta1, beta2, lamda1, lamda2)
par(mfrow=c(1,1),mar=c(4,3,1,1))
plot(results$time_keeper, results$x2_storage,type="l")
par(new=FALSE)
plot(results$time_keeper, results$x1_storage,type="l", col = "red")


