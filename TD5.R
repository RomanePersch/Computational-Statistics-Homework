######################## Problem 10.10 #####################
############################################################

### Question (b) : Gibbs-Sampler with alpha = 0.1 and a = b = 1

data_herds = c(rep(0,7), rep(1, 12), rep(2,8), rep(3,9), rep(4,8), rep(5,8), rep(6,9), rep(7,6),
               rep(8,5), rep(9,3), rep(10,4), rep(11,7), rep(12,4), rep(13,5), rep(14,2), rep(15,1), 
               rep(16,4), rep(17,3), rep(18,3), rep(19,4), rep(20,2), rep(21,2), rep(22,4), rep(23,1), 
               rep(25,6))

length(data_herds)

gibbs <- function(n_sim, alpha, a, b, x){
  #initialize
  result <- t(as.matrix(c(1, 1)))
  for (t in 2:n_sim){
    given_beta = result[(t-1),2] + 1
    lambda_t <- rgamma(1, shape = alpha + x, rate = given_beta)
    beta_t <- rgamma(1, shape = alpha + a, lambda_t + b)
    vector_t <- c(lambda_t, beta_t)
    result <- rbind(result, vector_t)
  }
  return(result)
}

#test with x = data_herds[1]
#test <- gibbs(100, 0.1, 1, 1,  data_herds[1])

gibbs_full_data <- function(n_sim, alpha, a, b, data_herds){
  parameters = list()
  for (i in 1:127){
    x = data_herds[i]
    parameters[[i]] <- gibbs(n_sim, alpha, a, b, x)
  }
  return(parameters)
}

#Results
n_sim = 10000
alpha = 0.1 
a = 1
b = 1

gibbs_questionb = gibbs_full_data(n_sim, alpha, a, b, data_herds)

result_1 <- matrix(unlist(gibbs_questionb[1]), ncol = 2, byrow = FALSE)

#Posterior of lambda 1
library(ggplot2)
qplot( result_1[,1], geom = 'histogram', bins = 100, fill = I("#FF6666"), xlab = 'Value of lambda 1')

#lambda 5
result_5 <- matrix(unlist(gibbs_questionb[5]), ncol = 2, byrow = FALSE)
qplot( result_5[,1], geom = 'histogram', bins = 100, fill = I("#FF6666"), xlab = 'Value of lambda 5')
#It is similar to lambda 1 since herd 5 has the same number of occurrences of clinical mastisis as herd 1...
#Compute the mean and monitor the convergence
est_vector_lambda5 = c()
for (k in 1:(n_sim)){
  est_mean <- mean(result_5[1:k,1])
  est_vector_lambda5 <- c(est_vector_lambda5, est_mean)
}
qplot(seq(1,n_sim-99,1),est_vector_lambda5[100:(n_sim)], geom = 'line', ylab = 'Mean of lambda 5', xlab = 'Number of simulations', colour = I("#FF6666"))


########### Herd 15 : 1 occurrence of mastisis in the herd ###################
#lambda 15
result_15 <- matrix(unlist(gibbs_questionb[15]), ncol = 2, byrow = FALSE)
qplot( result_15[,1], geom = 'histogram', bins = 100, fill = I("#33CC33"), xlab = 'Value of lambda 15')
#It is similar to lambda 1 since herd 5 has the same number of occurrences of clinical mastisis as herd 1...
#Compute the mean and monitor the convergence
est_vector_lambda15 = c()
for (k in 1:(n_sim)){
  est_mean <- mean(result_15[1:k,1])
  est_vector_lambda15 <- c(est_vector_lambda15, est_mean)
}
qplot(seq(1,n_sim-99,1),est_vector_lambda15[100:n_sim], geom = 'line', ylab = 'Mean of lambda 15', xlab = 'Number of simulations', colour = I("#33CC33"))

#beta 15 
qplot( result_15[,2], geom = 'histogram', bins = 100, fill = I("#33CC33"), xlab = 'Value of beta 15')
#Compute the mean and monitor the convergence
est_vector_beta15 = c()
for (k in 1:(n_sim)){
  est_mean <- mean(result_15[1:k,2])
  est_vector_beta15 <- c(est_vector_beta15, est_mean)
}
qplot(seq(1,n_sim-4,1),est_vector_beta15[5:n_sim], geom = 'line', ylab = 'Mean of beta 15', xlab = 'Number of simulations', colour = I("#33CC33"))


########### Herd 45 : 5 occurrences of mastisis in the herd ###################
#lambda 45
result_45 <- matrix(unlist(gibbs_questionb[45]), ncol = 2, byrow = FALSE)
qplot( result_45[,1], geom = 'histogram', bins = 100, fill = I("#0099FF"), xlab = 'Value of lambda 45')
#Compute the mean and monitor the convergence
est_vector_lambda45 = c()
for (k in 1:n_sim){
  est_mean <- mean(result_45[1:k,1])
  est_vector_lambda45 <- c(est_vector_lambda45, est_mean)
}
qplot(seq(1,n_sim-99,1),est_vector_lambda45[100:10000], geom = 'line', ylab = 'Mean of lambda 45', xlab = 'Number of simulations', colour = I("#0099FF"))

#beta 45 
qplot( result_45[,2], geom = 'histogram', bins = 100, fill = I("#0099FF"), xlab = 'Value of beta 45')
#Compute the mean and monitor the convergence
est_vector_beta45 = c()
for (k in 1:(n_sim)){
  est_mean <- mean(result_45[1:k,2])
  est_vector_beta45 <- c(est_vector_beta45, est_mean)
}
qplot(seq(1,n_sim-4,1),est_vector_beta45[5:n_sim], geom = 'line', ylab = 'Mean of beta 45', xlab = 'Number of simulations', colour = I("#0099FF"))


########### Herd 127 : 25 occurrences of mastisis in the herd ###################
#lambda 127
result_127 <- matrix(unlist(gibbs_questionb[127]), ncol = 2, byrow = FALSE)
qplot( result_127[,1], geom = 'histogram', bins = 100, xlab = 'Value of lambda 127')
#Compute the mean and monitor the convergence
est_vector_lambda127 = c()
for (k in 1:n_sim){
  est_mean <- mean(result_127[1:k,1])
  est_vector_lambda127 <- c(est_vector_lambda127, est_mean)
}
qplot(seq(1,n_sim-99,1),est_vector_lambda127[100:n_sim], geom = 'line', ylab = 'Mean of lambda 127', xlab = 'Number of simulations')

#beta 127 
qplot( result_127[,2], geom = 'histogram', bins = 100, xlab = 'Value of beta 127')
#Compute the mean and monitor the convergence
est_vector_beta127 = c()
for (k in 1:(n_sim)){
  est_mean <- mean(result_127[1:k,2])
  est_vector_beta127 <- c(est_vector_beta127, est_mean)
}
qplot(seq(1,n_sim,1),est_vector_beta127, geom = 'line', ylab = 'Mean of beta 127', xlab = 'Number of simulations')


### Question (c) : Impact of the hyperparameters a, b and alpha 
#We will check first Herd 15 and Herd 45 (instead of herd 5 and herd 15, which are too similar)

#################### Impact of alpha ############
n_sim = 2500
a = 1
b = 1

#On herd 15 (1 occurrenc of mastisis)
alpha_to_test = c(0.1,1,10)
results_alpha = list()
for (k in 1:length(alpha_to_test)){
  results_alpha[[k]] <- gibbs(n_sim, alpha_to_test[k], a, b, data_herds[15])
}

#lambda 15
mdf_hist <- data.frame(nb_sim = seq(1,n_sim,1))
mdf_mean <- data.frame(nb_sim = seq(1,n_sim,1))
for (k in 1:length(alpha_to_test)){
  result_15 <- matrix(unlist(results_alpha[k]), ncol = 2, byrow = FALSE)
  mdf_hist <- cbind(mdf_hist, result_15[,1])
  est_vector_lambda15 = c()
  for (j in 1:(n_sim)){
    est_mean <- mean(result_15[1:j,1])
    est_vector_lambda15 <- c(est_vector_lambda15, est_mean)
  }
  mdf_mean <- cbind(mdf_mean, est_vector_lambda15)
}

#Compare the mean
library("reshape2")
colnames(mdf_mean) <- c("nb_sim", "alpha = 0.1", "alpha = 1", "alpha = 10")
mdf2_mean <- melt(mdf_mean, id="nb_sim") 
ggplot(data=mdf2_mean,
       aes(x=nb_sim, y=value, colour=variable)) +
  geom_line() +ylab("Mean of lambda 15") +xlab("Simulations")

#Compare the histograms of the posterior
colnames(mdf_hist) <- c("nb_sim", "alpha1", "alpha2", "alpha3")
mdf2_hist <- melt(mdf_hist, id="nb_sim") 
ggplot(data=mdf2_hist,
       aes(x=value, fill=variable)) +
  geom_histogram(alpha = 0.25) 

#beta 15 
mdf_hist <- data.frame(nb_sim = seq(1,n_sim,1))
mdf_mean <- data.frame(nb_sim = seq(1,n_sim,1))
for (k in 1:length(alpha_to_test)){
  result_15 <- matrix(unlist(results_alpha[k]), ncol = 2, byrow = FALSE)
  mdf_hist <- cbind(mdf_hist, result_15[,2])
  est_vector_beta15 = c()
  for (j in 1:(n_sim)){
    est_mean <- mean(result_15[1:j,2])
    est_vector_beta15 <- c(est_vector_beta15, est_mean)
  }
  mdf_mean <- cbind(mdf_mean, est_vector_beta15)
}

#Compare the mean
library("reshape2")
colnames(mdf_mean) <- c("nb_sim", "alpha = 0.1", "alpha = 1", "alpha = 10")
mdf2_mean <- melt(mdf_mean, id="nb_sim") 
ggplot(data=mdf2_mean,
       aes(x=nb_sim, y=value, colour=variable)) +
  geom_line() +ylab("Mean of beta 15") +xlab("Simulations")

#Compare the histograms of the posterior
colnames(mdf_hist) <- c("nb_sim", "alpha1", "alpha2", "alpha3")
mdf2_hist <- melt(mdf_hist, id="nb_sim") 
ggplot(data=mdf2_hist,
       aes(x=value, fill=variable)) +
  geom_histogram(alpha = 0.25) 


#################### Impact of a ############
n_sim = 2500
alpha = 0.1
b = 1

#On herd 15 (1 occurrenc of mastisis)
a_to_test = c(0.1,1,10)
results_a = list()
for (k in 1:length(a_to_test)){
  results_a[[k]] <- gibbs(n_sim, alpha, a_to_test[k], b, data_herds[15])
}

#lambda 15
mdf_hist <- data.frame(nb_sim = seq(1,n_sim,1))
mdf_mean <- data.frame(nb_sim = seq(1,n_sim,1))
for (k in 1:length(a_to_test)){
  result_15 <- matrix(unlist(results_a[k]), ncol = 2, byrow = FALSE)
  mdf_hist <- cbind(mdf_hist, result_15[,1])
  est_vector_lambda15 = c()
  for (j in 1:(n_sim)){
    est_mean <- mean(result_15[1:j,1])
    est_vector_lambda15 <- c(est_vector_lambda15, est_mean)
  }
  mdf_mean <- cbind(mdf_mean, est_vector_lambda15)
}

#Compare the mean
library("reshape2")
colnames(mdf_mean) <- c("nb_sim", "a = 0.1", "a = 1", "a = 10")
mdf2_mean <- melt(mdf_mean, id="nb_sim") 
ggplot(data=mdf2_mean,
       aes(x=nb_sim, y=value, colour=variable)) +
  geom_line() +ylab("Mean of lambda 15") +xlab("Simulations")

#Compare the histograms of the posterior
colnames(mdf_hist) <- c("nb_sim", "a= 0.1", "a = 1", "a = 10")
mdf2_hist <- melt(mdf_hist, id="nb_sim") 
ggplot(data=mdf2_hist,
       aes(x=value, fill=variable)) +
  geom_histogram(alpha = 0.25) 

#beta 15 
mdf_hist <- data.frame(nb_sim = seq(1,n_sim,1))
mdf_mean <- data.frame(nb_sim = seq(1,n_sim,1))
for (k in 1:length(a_to_test)){
  result_15 <- matrix(unlist(results_a[k]), ncol = 2, byrow = FALSE)
  mdf_hist <- cbind(mdf_hist, result_15[,2])
  est_vector_beta15 = c()
  for (j in 1:(n_sim)){
    est_mean <- mean(result_15[1:j,2])
    est_vector_beta15 <- c(est_vector_beta15, est_mean)
  }
  mdf_mean <- cbind(mdf_mean, est_vector_beta15)
}

#Compare the mean
library("reshape2")
colnames(mdf_mean) <- c("nb_sim", "a = 0.1", "a = 1", "a = 10")
mdf2_mean <- melt(mdf_mean, id="nb_sim") 
ggplot(data=mdf2_mean,
       aes(x=nb_sim, y=value, colour=variable)) +
  geom_line() +ylab("Mean of beta 15") +xlab("Simulations")

#Compare the histograms of the posterior
colnames(mdf_hist) <- c("nb_sim", "a = 0.1", "a = 1", "a = 10")
mdf2_hist <- melt(mdf_hist, id="nb_sim") 
ggplot(data=mdf2_hist,
       aes(x=value, fill=variable)) +
  geom_histogram(alpha = 0.25) 


#################### Impact of b ############
n_sim = 2500
alpha = 0.1
a = 1

#On herd 15 (1 occurrenc of mastisis)
b_to_test = c(0.1,1,10)
results_b = list()
for (k in 1:length(b_to_test)){
  results_b[[k]] <- gibbs(n_sim, alpha, a, b_to_test[k], data_herds[15])
}

#lambda 15
mdf_hist <- data.frame(nb_sim = seq(1,n_sim,1))
mdf_mean <- data.frame(nb_sim = seq(1,n_sim,1))
for (k in 1:length(b_to_test)){
  result_15 <- matrix(unlist(results_b[k]), ncol = 2, byrow = FALSE)
  mdf_hist <- cbind(mdf_hist, result_15[,1])
  est_vector_lambda15 = c()
  for (j in 1:(n_sim)){
    est_mean <- mean(result_15[1:j,1])
    est_vector_lambda15 <- c(est_vector_lambda15, est_mean)
  }
  mdf_mean <- cbind(mdf_mean, est_vector_lambda15)
}

#Compare the mean
library("reshape2")
colnames(mdf_mean) <- c("nb_sim", "b = 0.1", "b = 1", "b = 10")
mdf2_mean <- melt(mdf_mean, id="nb_sim") 
ggplot(data=mdf2_mean,
       aes(x=nb_sim, y=value, colour=variable)) +
  geom_line() +ylab("Mean of lambda 15") +xlab("Simulations")

#Compare the histograms of the posterior
colnames(mdf_hist) <- c("nb_sim", "b= 0.1", "b = 1", "b = 10")
mdf2_hist <- melt(mdf_hist, id="nb_sim") 
ggplot(data=mdf2_hist,
       aes(x=value, fill=variable)) +
  geom_histogram(alpha = 0.25) 

#beta 15 
mdf_hist <- data.frame(nb_sim = seq(1,n_sim,1))
mdf_mean <- data.frame(nb_sim = seq(1,n_sim,1))
for (k in 1:length(b_to_test)){
  result_15 <- matrix(unlist(results_b[k]), ncol = 2, byrow = FALSE)
  mdf_hist <- cbind(mdf_hist, result_15[,2])
  est_vector_beta15 = c()
  for (j in 1:(n_sim)){
    est_mean <- mean(result_15[1:j,2])
    est_vector_beta15 <- c(est_vector_beta15, est_mean)
  }
  mdf_mean <- cbind(mdf_mean, est_vector_beta15)
}

#Compare the mean
library("reshape2")
colnames(mdf_mean) <- c("nb_sim", "b = 0.1", "b = 1", "b = 10")
mdf2_mean <- melt(mdf_mean, id="nb_sim") 
ggplot(data=mdf2_mean,
       aes(x=nb_sim, y=value, colour=variable)) +
  geom_line() +ylab("Mean of beta 15") +xlab("Simulations")

#Compare the histograms of the posterior
colnames(mdf_hist) <- c("nb_sim", "b = 0.1", "b = 1", "b = 10")
mdf2_hist <- melt(mdf_hist, id="nb_sim") 
ggplot(data=mdf2_hist,
       aes(x=value, fill=variable)) +
  geom_histogram(alpha = 0.25,bins = 40) 



