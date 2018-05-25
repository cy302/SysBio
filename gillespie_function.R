#x1 is mRNA which has low amount, 10
#x2 is protein which is abundant, 1000
##### Initialisation #####
x1 <- 10
x2 <- 1000
iteration <- 1e4
K <- 500
beta1 <- 1
beta2 <- 1
lamda1 <- 100
lamda2 <- 100

#auto regulation, suppose reduction of protein could lead to increase of mRNA production rate
gillespie <- function(){
    x2_death_rate <- x2 * beta2
    x2_birth_rate <-  x1 * lamda2
    
    x1_birth_rate <-  lamda1*(K/(K+x2)) 
    x1_death_rate <- x1 * beta1
    
    Tot_rate1 <- sum(x1_birth_rate,x1_death_rate)
    Tot_rate2 <- sum(x2_birth_rate,x2_death_rate)
    time_keeper1 <- time_keeper2 <- c()
    x1_stroage <- x2_stroage <-c()
    counter1 <- 1
    counter2 <- 1
    x1_stroage[counter1] <- x1
    x2_stroage[counter2] <- x2
    time_keeper1[counter1] <- rexp(1,Tot_rate1)
    time_keeper2[counter2] <- rexp(1,Tot_rate2)
    
    for (i in 1:iteration){
        
        #The timeline of protein catches up the timeline of mRNA
        while(time_keeper1[counter1] > time_keeper2[counter2]){ 
            x2_death_rate <- x2 * beta2
            x2_birth_rate <-  x1 * lamda2
            Tot_rate2 <- sum(x2_birth_rate,x2_death_rate)
            x2 <- ifelse(runif(1) < x2_birth_rate/Tot_rate2, x2+1, x2-1)
            x2_stroage[counter2 + 1] <- x2
            time_keeper2[counter2 + 1] <- rexp(1,Tot_rate2) + time_keeper2[counter2]
            counter2 <- counter2 + 1
        }
        
        #The timeline of mRNA catches up the timeline of protein
        while(time_keeper1[counter1] < time_keeper2[counter2]){
            x1_birth_rate <-  lamda1*(K/(K+x2)) 
            x1_death_rate <- x1 * beta1
            Tot_rate1 <- sum(x1_birth_rate,x1_death_rate)
            x1 <- ifelse(runif(1) < x1_birth_rate/Tot_rate1, x1+1, x1-1)
            x1_stroage[counter1 + 1] <- x1
            time_keeper1[counter1 + 1] <- rexp(1,Tot_rate1) + time_keeper1[counter1]
            counter1 <- counter1 + 1
        }
    }
}

par(mfrow=c(1,1),mar=c(4,3,1,1))
plot(time_keeper2[seq(0,length(time_keeper2),1000)],
     x2_stroage[seq(0,length(time_keeper2),1000)],type="l")
par(new=TRUE)
plot(time_keeper1[seq(0,length(time_keeper1),1)],
     x1_stroage[seq(0,length(time_keeper1),1)],col="red",type="l")


