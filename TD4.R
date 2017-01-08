################################ Problem 7.1 ################################
#############################################################################
library(ggplot2)


#### Question (a) : Accept-Reject

#Gamma simulation with alpha an integer
transform_exp <- function(u, beta){
  epsilon <- - log(u)*beta
  return(epsilon)
}

gamma_simul_integer <- function(n_sim = 500, alpha,beta){
  gamma_simul = rep(0, n_sim)
  for (k in 1:alpha){
    u <- runif(n_sim, 0, 1)
    epsilon <- sapply(u,transform_exp, beta = beta)
    gamma_simul <- gamma_simul + epsilon
  }
  return(gamma_simul)
}


quotient_f_g <- function(y,alpha, beta, a, b){
  result = y**(alpha - a) * exp(alpha - a)
  result = result * exp(y * ((beta - b)/ (b*beta)))
  result = result * (((b - beta) /((alpha - a)*(b*beta)) )**(alpha - a))
  return(result)
}



n_sim = 10000


y <- gamma_simul_integer(n_sim = n_sim, alpha = 4, beta = 7) 
#Check if it really looks like a Gamma density:
qplot(y , geom = 'density')

unif <- runif(n_sim, min = 0, max = 1)
accepted <- ifelse(unif <= sapply(y,quotient_f_g,alpha = 4.3, beta = 6.2, a = 4, b = 7),TRUE, FALSE)
print(length(accepted))
summary(accepted)

#Check if the empirical acceptance rate is close to the unconditional acceptance probability (equal to 1 / M)
compute_inv_M <- function(alpha,beta,a,b){
  result <- exp(alpha - a) * (((b - beta) /((alpha - a)*(b*beta)) )**(alpha - a))
  result <- result * gamma(alpha) * (beta**(alpha)) * (1 / (gamma(a)*(b**(a)) ))
  return(result)
}

compute_inv_M(alpha = 4.3, beta = 6.2, a = 4, b = 7)

qplot(y[accepted], geom = 'histogram', xlab='Accepted values', bins = 100)
print(length(y[accepted])) 

#Compute the mean and monitor the convergence
est_vector_AR = c()
est_mean = 1 #arbitrary value initialization
for (k in 1:n_sim){
  est_mean <- ifelse(accepted[k] == TRUE,mean(y[1:k][accepted[1:k]]),est_mean)
  est_vector_AR <- c(est_vector_AR, est_mean)
}


q_AR <- qplot(seq(1,n_sim,1),est_vector_AR[1:n_sim], geom = 'line', ylab = 'Mean of Gamma(4.3,6.2)', xlab = 'Overall number of simulations')
q_AR
qplot(seq(1,length(accepted[accepted == TRUE]),1),est_vector_AR[accepted == TRUE], geom = 'line', ylab = 'Mean of Gamma(4.3,6.2)', xlab = 'Number of simulations')

#final value of the mean
mean(y[accepted])

#### Question (b) : Metropolis-Hastings with candidate = Gamma(4,7)

rho_accept <- function(x,y, alpha, beta, a ,b){
  quotient <- ((y**(alpha-a))/(x**(alpha-a))) * exp(((beta -b)/(b*beta))*(y - x))
  return(min(1,quotient))
}

independent_metropolis_hastings_gamma <-function(n_sim, alpha, beta, a, b){
  y <- gamma_simul_integer(n_sim, a, b)
  #arbitrary initialization
  x <- c(20)
  acceptance <- 0
  for (k in 1:n_sim){
    u <- runif(1,0,1)
    rho <- rho_accept(x[k],y[k], alpha, beta, a ,b)
    if (u < rho){
      x <- c(x,y[k])
      acceptance <- acceptance + 1
    } else{
      x<-c(x, x[k])
    }
  }
  acceptance <- acceptance / n_sim
  return(list(x, acceptance))
}

n_sim= 10000
results <- independent_metropolis_hastings_gamma(n_sim, alpha = 4.3, beta=6.2, a=4, b=7)
x_markov <- unlist(results[1])
acceptance_rate <- unlist(results[2])
print(acceptance_rate)

#Plot the autocorrelations to monitor convergence towards the stationary distribution
bacf1 <- acf(x_markov, plot = FALSE)
bacfdf1 <- with(bacf1, data.frame(lag, acf))
q1 <- ggplot(data=bacfdf1, mapping=aes(x=lag, y=acf)) +
  geom_area(fill="grey") +
  geom_hline(yintercept=c(0.05, -0.05), linetype="dashed") +
  theme_bw() + ggtitle("ACF")

q1

#Plot the Markov chain
qplot(seq(1,n_sim + 1,1), x_markov, geom = 'line', ylab = 'Value' )

#Compute the mean and monitor the convergence
est_vector_MH1 = c()
for (k in 1:(n_sim+1)){
  est_mean <- mean(x_markov[1:k])
  est_vector_MH1 <- c(est_vector_MH1, est_mean)
}
q_MH1 <- qplot(seq(1,n_sim,1),est_vector_MH1[2:(n_sim + 1)], geom = 'line', ylab = 'Mean of Gamma(4.3,6.2)', xlab = 'Number of simulations')
q_MH1

#Final value of the mean
est_vector_MH1[(n_sim + 1)]

#### Question (c) : Metropolis-Hastings with candidate = Gamma(5,6)
n_sim= 10000
results <- independent_metropolis_hastings_gamma(n_sim, alpha = 4.3, beta=6.2, a=5, b=6)
x_markov <- unlist(results[1])
acceptance_rate <- unlist(results[2])
print(acceptance_rate)

#Plot the autocorrelations to monitor convergence towards the stationary distribution
bacf2 <- acf(x_markov, plot = FALSE)
bacfdf2 <- with(bacf2, data.frame(lag, acf))
q2 <- ggplot(data=bacfdf2, mapping=aes(x=lag, y=acf)) +
  geom_area(fill="grey") +
  geom_hline(yintercept=c(0.05, -0.05), linetype="dashed") +
  theme_bw() + ggtitle("ACF")

q2

#Plot the Markov chain
qplot(seq(1,n_sim + 1,1), x_markov, geom = 'line', ylab = 'Value' )

#Compute the mean and monitor the convergence
est_vector_MH2 = c()
for (k in 1:(n_sim+1)){
  est_mean <- mean(x_markov[1:k])
  est_vector_MH2 <- c(est_vector_MH2, est_mean)
}
q_MH2 <- qplot(seq(1,n_sim,1),est_vector_MH2[2:(n_sim + 1)], geom = 'line', ylab = 'Mean of Gamma(4.3,6.2)', xlab = 'Number of simulations')
q_MH2

#Final value of the mean
est_vector_MH2[(n_sim + 1)]

#Compare algorithms
mdf <- data.frame(AR = est_vector_AR, MH1 = est_vector_MH1[2:(n_sim + 1)], MH2 =est_vector_MH2[2:(n_sim + 1)] , Simulation = seq(1,n_sim,1))

library("reshape2")
mdf2 <- melt(mdf, id="Simulation") 
ggplot(data=mdf2,
       aes(x=Simulation, y=value, colour=variable)) +
  geom_line()


#Candidate densities and "target" (using rgamma)
x <- data.frame(target=rgamma(100000, shape = 4.3, scale = 6.2 ),candidate1=rgamma(100000, shape = 4, scale = 7 ),candidate2=rgamma(100000, shape = 5, scale = 6 ))
data<- melt(x)
ggplot(data,aes(x=value, colour=variable)) + geom_density(alpha=0.25)

################################ Problem 7.21 ################################
#############################################################################

#Create the dataframe containing the braking data
x <- c(rep(4,2), rep(7,2), 8, 9, rep(10, 3), rep(11,2), rep(12,4),rep(13,4), 
       rep(14,4), rep(15,3), rep(16,2), rep(17,3), rep(18,4), rep(19,3), 
       rep(20, 5), 22, 23, rep(24,4), 25)
y<- c(2,10, 4,22, 16, 10, 18,26, 34, 17,28, 14,20, 24,28, 26,34, 34,46, 26,36, 60,80, 20,26, 54,
      32,40, 32,40, 50, 42,56, 76,84, 36,46, 68, 32,48, 52,56,64, 66, 54, 70, 92, 93, 120, 85) 

braking_data <- data.frame(x = x, y = y, x2 = x**2)

#### Question (a) : Linear regression

lin_reg <- lm(y ~ x + x2, data = braking_data)
summary(lin_reg)

coefficients_vector <- coefficients(lin_reg)
sigma <- summary(lin_reg)$sigma

theta_est_linreg <- c(coefficients_vector,sigma**2)
theta_est_linreg

### Question (b) & (c) : Metropolis-Hastings

#compute the log likelihood on the braking data for a given parameter theta (since the likelihood is too small for R !)
log_likelihood <- function(theta){
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  sigma2 <- theta[4]
  N = dim(braking_data)[1]
  vector_for_sum <- (braking_data$y - a - b*braking_data$x -c*braking_data$x2)**2
  result = (N/2)*log(1/sigma2) -(1/(2*sigma2))*sum(vector_for_sum)
  return(result)
}

#compute the acceptance probability
rho_accept_log <- function(x,y){
  quotient <- log_likelihood(y) - log_likelihood(x)
  result <- min(0, quotient)
  return(result)
}


#Independent Metropolis Hastings with candidate = prior (here 3 Normal distributions and 1 Inverse gamma)
independent_metropolis_hastings_posterior <-function(n_sim, theta_est_linreg, var_a, var_b, var_c, beta_sigma2){
  a <- rnorm(n_sim, theta_est_linreg[1], var_a)
  b <- rnorm(n_sim, theta_est_linreg[2], var_b)
  c <- rnorm(n_sim, theta_est_linreg[3], var_c)
  shape_inv_gamma = 1 + (beta_sigma2/theta_est_linreg[4])
  sigma2 <- 1/rgamma(n_sim, shape = shape_inv_gamma, rate = beta_sigma2)
  y <-cbind(a,b,c,sigma2)
  #arbitrary initialization
  x <- t(as.matrix(theta_est_linreg))
  acceptance <- 0
  for (k in 1:n_sim){
    u <- runif(1,0,1)
    rho <- rho_accept_log(x[k,],y[k,])
    if (log(u) < rho){
      x <- rbind(x,y[k,])
      acceptance <- acceptance + 1
    } else{
      x<-rbind(x, x[k,])
    }
  }
  acceptance <- acceptance / n_sim
  return(list(x, acceptance))
}

n_sim= 50000
results <- independent_metropolis_hastings_posterior(n_sim = n_sim, theta_est_linreg = theta_est_linreg, var_a= 0.5, var_b= 0.05, var_c = 0.01, beta_sigma2 = 250)
theta_markov <- matrix(unlist(results[1]), ncol = 4, byrow = FALSE)
acceptance_rate <- unlist(results[2])
print(acceptance_rate)

library(ggplot2)
#Plot the posterior
#Posterior of a
qplot( theta_markov[,1], geom = 'histogram', bins = 100, fill = I("#FF6666"), xlab = 'Value of a')  
qplot(seq(1,n_sim + 1,1), theta_markov[,1], geom = 'line', ylab = 'Value' ,xlab = "", colour =I("#FF6666"))
#Posterior of b
qplot( theta_markov[,2], geom = 'histogram', bins = 100, fill = I("#33CC33"), xlab = 'Value of b' )
qplot(seq(1,n_sim + 1,1), theta_markov[,2], geom = 'line', ylab = 'Value', xlab = "", colour =  I("#33CC33")  )
#Posterior of c
qplot( theta_markov[,3], geom = 'histogram', bins = 100, fill = I("#0099FF"), xlab = 'Value of c' )
qplot(seq(1,n_sim + 1,1), theta_markov[,3], geom = 'line', ylab = 'Value', xlab = "", colour = I("#0099FF") )
#Posterior of sigma2
qplot( theta_markov[,4], geom = 'histogram', bins = 100 , xlab = 'Value of sigma2')
qplot(seq(1,n_sim + 1,1), theta_markov[,4], geom = 'line', ylab = 'Value', xlab = "" )

### Question (d) : Metropolis-Hastings for the second model

#compute the log likelihood on the braking data for a given parameter theta (since the likelihood is too small for R !)
log_likelihood2 <- function(theta, nu = 4){
  a <- theta[1]
  b <- theta[2]
  c <- theta[3]
  sigma2 <- theta[4]
  N = dim(braking_data)[1]
  vector_for_sum <- 1 + ((braking_data$y - a - b*braking_data$x -c*braking_data$x2)**2)/nu
  result = (N/2)*log(1/sigma2) -((nu + 1)/2)*sum(sapply(vector_for_sum, log))
  return(result)
}

log_likelihood2(theta_est_linreg)

#compute the acceptance probability
rho_accept_log2 <- function(x,y){
  quotient <- log_likelihood2(y) - log_likelihood2(x)
  result <- min(1, quotient)
  return(result)
}


#Independent Metropolis Hastings with candidate = prior
#NORMAL INVERSE GAMMA PRIOR
independent_metropolis_hastings_posterior2 <-function(n_sim, theta_est_linreg, var_a, var_b, var_c, beta_sigma2){
  a <- rnorm(n_sim, theta_est_linreg[1], var_a)
  b <- rnorm(n_sim, theta_est_linreg[2], var_b)
  c <- rnorm(n_sim, theta_est_linreg[3], var_c)
  shape_inv_gamma = 1 + (beta_sigma2/theta_est_linreg[4])
  sigma2 <- 1/rgamma(n_sim, shape = shape_inv_gamma, rate = beta_sigma2)
  y <-cbind(a,b,c,sigma2)
  #arbitrary initialization
  x <- t(as.matrix(theta_est_linreg))
  acceptance <- 0
  for (k in 1:n_sim){
    u <- runif(1,0,1)
    rho <- rho_accept_log2(x[k,],y[k,])
    if (log(u) < rho){
      x <- rbind(x,y[k,])
      acceptance <- acceptance + 1
    } else{
      x<-rbind(x, x[k,])
    }
  }
  acceptance <- acceptance / n_sim
  return(list(x, acceptance))
}


#We try to adjust the parameters...
n_sim= 10000
results <- independent_metropolis_hastings_posterior2(n_sim = n_sim, theta_est_linreg = theta_est_linreg, var_a= 2, var_b= 1, var_c = 0.25, beta_sigma2 = 250)
theta_markov <- matrix(unlist(results[1]), ncol = 4, byrow = FALSE)
acceptance_rate <- unlist(results[2])
print(acceptance_rate)

library(ggplot2)
#Plot the posterior (NORMAL INVERSE GAMMA PRIOR)
#Posterior of a
qplot( theta_markov[,1], geom = 'histogram', bins = 100, fill = I("#FF6666"), xlab = 'Value of a')  
qplot(seq(1,n_sim + 1,1), theta_markov[,1], geom = 'line', ylab = 'Value' ,xlab = "", colour =I("#FF6666"))
#Posterior of b
qplot( theta_markov[,2], geom = 'histogram', bins = 100, fill = I("#33CC33"), xlab = 'Value of b' )
qplot(seq(1,n_sim + 1,1), theta_markov[,2], geom = 'line', ylab = 'Value', xlab = "", colour =  I("#33CC33")  )
#Posterior of c
qplot( theta_markov[,3], geom = 'histogram', bins = 100, fill = I("#0099FF"), xlab = 'Value of c' )
qplot(seq(1,n_sim + 1,1), theta_markov[,3], geom = 'line', ylab = 'Value', xlab = "", colour = I("#0099FF") )
#Posterior of sigma2
qplot( theta_markov[,4], geom = 'histogram', bins = 100 , xlab = 'Value of sigma2')
qplot(seq(1,n_sim + 1,1), theta_markov[,4], geom = 'line', ylab = 'Value', xlab = "" )


                      ######################

#Independent Metropolis Hastings with candidate = prior
#t and half-t candidates
#Note that here we don't not use the linear regression anymore to choose the prior parameters since the mean of a t-distribution is 0
#We only use the OLS estimators as an initialization
independent_metropolis_hastings_posterior3 <-function(n_sim, theta_est_linreg, nu_a, nu_b, nu_c, nu_sigma2){
  a <- rt(n_sim, nu_a)
  b <- rt(n_sim, nu_b)
  c <- rt(n_sim, nu_c)
  sigma2 <- sapply(rt(n_sim, nu_sigma2), abs)
  y <-cbind(a,b,c,sigma2)
  #arbitrary initialization
  x <- t(as.matrix(theta_est_linreg))
  acceptance <- 0
  for (k in 1:n_sim){
    u <- runif(1,0,1)
    rho <- rho_accept_log2(x[k,],y[k,])
    if (log(u) < rho){
      x <- rbind(x,y[k,])
      acceptance <- acceptance + 1
    } else{
      x<-rbind(x, x[k,])
    }
  }
  acceptance <- acceptance / n_sim
  return(list(x, acceptance))
}


#We try to adjust the parameters...
n_sim= 15000
results <- independent_metropolis_hastings_posterior3(n_sim = n_sim, theta_est_linreg = theta_est_linreg, nu_a = 3.5, nu_b= 10, nu_c= 50, nu_sigma2 = 0.95)
theta_markov <- matrix(unlist(results[1]), ncol = 4, byrow = FALSE)
acceptance_rate <- unlist(results[2])
print(acceptance_rate)
mean(theta_markov[,4])
mean(theta_markov[,3])
mean(theta_markov[,2])
mean(theta_markov[,1])