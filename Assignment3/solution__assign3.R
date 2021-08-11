#############
#assignment3#
#############
library(ggplot2)
library(gplots)

#############
#question 1 to 4 are on the pdf and below it is the code for question5
#############

Y_i <- c(162, 267, 271, 185, 111, 61, 27, 8, 3, 1)
i <- c(0:9)

# Initialize setting
p <- 0
lambda_1 <- 0
lambda_2 <- 0

# Set initial values
p <- c(p, 0.01)
lambda_1 <- c(lambda_1, 0.001)
lambda_2 <- c(lambda_2, 0.001)

p_prim <- rep(0, 10)

# Run algorithm
k <- 3

# E step (computing p_prim)
  for (j in 1:10){
    p_prim[j] <- p[k-1] * dpois(j-1, lambda_1[k-1]) / ((p[k-1] * dpois(j-1, lambda_1[k-1]) + (1 - p[k-1]) * dpois(j-1, lambda_2[k-1])))
  }
  
# M step
  p[k] <- sum(Y_i * p_prim) / sum(Y_i)
  lambda_1[k] <- sum(Y_i * i * p_prim) / (1 + sum(Y_i * p_prim))
  lambda_2[k] <- sum(Y_i * i * (1 - p_prim)) / (1 + sum(Y_i * (1 - p_prim)))               
  
theta <- c(p,lambda_1,lambda_2)

##### Plots #####

# MLE of lambda
lambda_est <- sum(Y_i*i) / sum(Y_i)
n <- sum(Y_i)
est_vals <- n * dpois(i, lambda_est)

# Counts from model fitted in assignment 3

extended_model <- p[3] * dpois(i, lambda_1[3]) + (1 - p[3]) * dpois(i, lambda_2[3])

ggplot() +
  geom_line(
    aes(x = i, y =  n * extended_model), colour = "red") +#Extended model
  geom_line(
    aes(x = i, y = Y_i), colour = "black") +#Actual counts
  labs(title = " Extended model (red line) and actual counts (black line)",
       x = "i",
       y = "Density")
