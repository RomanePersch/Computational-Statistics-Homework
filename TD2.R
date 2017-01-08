################################ Problem 3.1 ################################
#############################################################################

#### Question (a)

#Plot the integrand
integrand_num <- function(theta,x){
  result <- (theta / (1 + theta**2))*exp(-0.5*((x - theta)**2))
  return(result)
}

integrand_den <- function(theta,x){
  result <- (1 / (1 + theta**2))*exp(-0.5*((x - theta)**2))
  return(result)
}

x = 0

library("reshape2")
library("ggplot2")


num_df <-
  data.frame(
    theta = seq(-10,10,0.1),
    numerator1 = sapply(seq(-10,10,0.1),integrand_num, x = 0),
    numerator2 = sapply(seq(-10,10,0.1),integrand_num, x = 1),
    numerator3 = sapply(seq(-10,10,0.1),integrand_num, x = 5)
  )

den_df <-
  data.frame(
    theta = seq(-10,10,0.1),
    denominator1 = sapply(seq(-10,10,0.1),integrand_den, x = 0),
    denominator2 = sapply(seq(-10,10,0.1),integrand_den, x = 1),
    denominator3 = sapply(seq(-10,10,0.1),integrand_den, x = 5)
  )


num_df_long <- melt(num_df, id="theta")  # convert to long format
den_df_long <- melt(den_df, id="theta") 

ggplot(data=num_df_long,
       aes(x=theta, y=value, colour=variable, xlab='theta')) +
  geom_line() + scale_colour_discrete(name  ="Value of x",
                                     breaks=c("numerator1", "numerator2", "numerator3"),
                                     labels=c("x = 0", "x = 1", "x = 5"))

ggplot(data=den_df_long,
       aes(x=theta, y=value, colour=variable, xlab='theta')) +
  geom_line() + scale_colour_discrete(name  ="Value of x",
                                      breaks=c("denominator1", "denominator2", "denominator3"),
                                      labels=c("x = 0", "x = 1", "x = 5"))

# Monte Carlo Simulation of Delta(x)
transform_num <- function(u){
  return(u / (1 + u**2))
}

transform_den <- function(u){
  return(1 / (1 + u**2))
}

theta_sim <- function(m,x){
  return(rnorm(m, mean = x, sd = 1))
}

delta_estimator <-function(theta){
  #Returns Delta(x) estimator using the vector of theta simulations
  #Note that we use the same theta simulations to compute the numerator and the denominator
  delta_x_num <- sum(sapply(theta, transform_num))
  delta_x_den <- sum(sapply(theta, transform_den))
  return(delta_x_num / delta_x_den)
}

m_max = 1000

x = 5
replications_sim = list()
replications_est = list()
for (i in 1:100){
  replications_sim[[i]] <- theta_sim(m_max,x)
  est_vector = c()
  for (k in 1:m_max){
    delta <- delta_estimator(replications_sim[[i]][1:k])
    est_vector <- c(est_vector, delta)
  }
  replications_est[[i]] <- est_vector
}

replications_df <- h <- do.call(cbind, replications_est)
print(dim(replications_df))

replications_max <- apply(replications_df, 1, max)
replications_min <- apply(replications_df, 1, min)

df_plot <- data.frame(
  m = seq(1,m_max,1),
  replications_max = replications_max,
  replications_min = replications_min,
  replications_mean <- apply(replications_df, 1, mean)
)

ggplot(data = df_plot,aes(m, replications_mean))+
  geom_ribbon(aes(x=m, ymax=replications_max, ymin=replications_min), fill="pink", alpha=.5) +
  geom_line(aes(y = replications_max), colour = 'black') +
  geom_line(aes(y = replications_min), colour = 'black')+
  geom_line() + xlab('m') + ylab('Delta(x) estimator')

#### Question (b)


standard_error_estimation <- function(first_part,second_part){
  delta <- sum(first_part) / sum(second_part)
  sum_elements <- (first_part - delta * second_part)**2
  m = length(second_part)
  alpha <- mean(second_part)
  delta_variance <- (1/ m**2)*(1/alpha**2)*sum(sum_elements)
  delta_sd <- sqrt(delta_variance)
  return(c(delta, delta_sd))
}


m_max = 20000
x = 0

theta <- theta_sim(m_max,x)
first_part <- sapply(theta, transform_num)
second_part <- sapply(theta, transform_den)
evolving_var = c()
evolving_delta = c()
for (k in 2:m_max){
  results <- standard_error_estimation(first_part[1:k], second_part[1:k])
  evolving_var <- c(evolving_var, results[2])
  evolving_delta <- c(evolving_delta, results[1])
}

qplot(seq(2,m_max,1), evolving_var, xlab='m', ylab='Standard Error') 
qplot(seq(2,m_max,1), evolving_delta, xlab='m', ylab='Delta(0)') 



### /!\ We only make a two-digits approximation for computational reasons
#Check when the standard error becomes lower than 0.01/1.96
test <- (evolving_var < (0.01/qnorm(0.975)))
min_m <- min(which(test == TRUE))
print(min_m + 1) #result : 13963
print(evolving_var[min_m])
print( evolving_delta[min_m])

#0.95 CI is then :
print(evolving_delta[min_m] + evolving_var[min_m]*qnorm(0.975))
print(evolving_delta[min_m] - evolving_var[min_m]*qnorm(0.975))

################################ Problem 3.2 ################################
#############################################################################

#### Question (a)

quotient_f_g <- function(theta,x){
  return(exp(-(x - theta)**2/2))
}


n_sim = 100000
x = 5

theta <- rcauchy(n_sim, location = 0, scale = 1)
unif <- runif(n_sim, min = 0, max = 1)
accepted <- ifelse(unif <= sapply(theta,quotient_f_g,x = x),TRUE, FALSE)
print(length(accepted))
summary(accepted)

qplot(theta[accepted], geom = 'histogram', xlab='Accepted theta', main = 'Theta posterior distribution when x = 1', bins = 100)
print(length(theta[accepted])) 

#Estimation of Delta(x):
est_vector = c()
delta = 10 #arbitrary value initialization
for (k in 1:n_sim){
  delta <- ifelse(accepted[k] == TRUE,mean(theta[accepted[1:k]]),delta)
  est_vector <- c(est_vector, delta)
}

qplot(seq(1,1000,1),est_vector[1:1000], geom = 'line')

qplot(seq(1,100000,1),est_vector[1:100000], geom = 'line', ylab = 'Delta(5)', xlab = 'Number of prior simulations')
qplot(seq(1,length(accepted[accepted == TRUE]),1),est_vector[accepted == TRUE], geom = 'line')

#### Question (b)

m_max = 1000
x = 0

theta <- theta_sim(m_max,x)
first_part <- sapply(theta, transform_num)
second_part <- sapply(theta, transform_den)
evolving_var_samesim = c()
for (k in 2:m_max){
  delta_sd <- standard_error_estimation(first_part[1:k], second_part[1:k])[2]
  evolving_var_samesim <- c(evolving_var_samesim, delta_sd)
}

theta2 <- theta_sim(m_max,x)
second_part2 <- sapply(theta2, transform_den)
evolving_var_diffsim = c()
for (k in 2:m_max){
  delta_sd <- standard_error_estimation(first_part[1:k], second_part2[1:k])[2]
  evolving_var_diffsim <- c(evolving_var_diffsim, delta_sd)
}

compare_var_df <- data.frame(
  m = seq(2,m_max,1),
  same_sim = evolving_var_samesim,
  diff_sim = evolving_var_diffsim)

compare_var_df_long <- melt(compare_var_df, id="m") 

ggplot(data=compare_var_df_long,
       aes(x=m, y=value, colour=variable, xlab='m')) +
  geom_line() + ylab('Standard Error Estimation') + scale_colour_discrete(name  ="",
                                      breaks=c("same_sim", "diff_sim"),
                                      labels=c("Same simulations", "Different simulations"))

################################ Problem 3.10 ################################
##############################################################################
library(ggplot2)

### Classic Monte Carlo estimation
n_sim = 1000

h <- function(x){
  res <- exp(-(1/2)*(x**2))/sqrt(2*pi)
  return(res)
}


classic_mc <- function(n_sim){
  x <- runif(n_sim, min = 1, max = 2)
  h_x <- sapply(x, h)
  I_MC_evol <- c()
  for (i in 1:n_sim){
    I_MC <- mean(h_x[1:i])
    I_MC_evol <- c(I_MC_evol, I_MC)
  }
  return(I_MC_evol)
}

results <- classic_mc(n_sim)
qplot(seq(1,n_sim,1), results, geom = 'line', xlab = 'Number of simulations', ylab = 'I estimation')

### Importance Sampling using optimal instrumental distribution
