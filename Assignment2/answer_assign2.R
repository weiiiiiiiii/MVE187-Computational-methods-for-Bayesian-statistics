#############
#assignment2#
#############
library(ggplot2)
library(gplots)

#1a)
#Pr(p(v))>0=Pr(4<v<25)
prob <- pgamma(25, shape = 0.8, scale = 6)- pgamma(4, shape = 0.8, scale = 6)
prob #[1] 0.4019568

#1b)
#define the p function
P <- function(v){
  n <- length(v)
  P_temp <- rep(0, n)
  for (i in 1:n){
    if (v[i] >= 4 && v[i] < 15 ){
      P_temp[i] <- cos(v[i]*pi/11+7*pi/11) + 1}
    if (v[i] >= 15 && v[i] < 25){
      P_temp[i] <- 7/4 + v[i]/30 - v[i]^2/900}}
  return(P_temp)
}
# Compute approximate expected value
N <- 10000
sample <- rgamma(N, shape = 0.8, scale = 6)
P_sample <- P(sample)
approx <- 1/N * sum(P_sample)
# Compute approximate 95% CI
s <- sqrt(1 / (N-1)  * sum((P_sample - approx)^2))
ci_lower <- approx - 1.96 * s / sqrt(N)
ci_upper <- approx + 1.96 * s / sqrt(N)
print(mean(approx))
#0.3084502
print(c(ci_lower,ci_upper))
#[1] 0.2883652 0.3114523

#1c)
# Define integrand
integrand <- function(v){
  p <- P(v)
  gamma <- dgamma(v, shape = 0.8, scale = 6)
  return(p*gamma)
}
# Compute integral numerically
E <- integrate(integrand, 4, 25)
#0.3057371 with absolute error < 3.3e-05
x <- seq(4, 25, length.out=10001)
y <- P(x)*dgamma(x, 0.8, 1/6)*(x[2]-x[1])
print(sum(y))
#[1] 0.3057338

#1d)
v <- seq(0, 30, length.out = 10000)
d1 <- P(v)*dgamma(v, 0.8, 1/6)
d2 <- dnorm(v, mean=11, sd=6)
df <- data.frame(x,d1,d2)
ggplot(df, aes(v)) +                    # basic graphical object
  geom_line(aes(y=d1), colour="red") +  # first layer
  geom_line(aes(y=d2), colour="black")+  # second layer
  labs(
    y="density",
    x="v",
    caption = "Product of densities of P(v) and v and choosen proposal density.")






#1e)
# Define importance sampling function
P_importance <- function(v){
  p_v <- P(v)
  pi_v <- dgamma(v, shape = 0.8, scale = 6)
  chosen_v <- dnorm(v, mean = 11, sd = 6, log = FALSE)
  return(p_v * pi_v / chosen_v)
}

N <- 10000
importance_sample <- rnorm(N, mean = 11, sd = 6)
P_importance_sample <- P_importance(importance_sample)
importance_approx <- 1/N * sum(P_importance_sample)

# 95% confidence interval
s_importance <- sqrt(1 / (N-1)  * sum((P_importance_sample - importance_approx)^2))

ci_importance_lower <- importance_approx - 1.96 * s_importance / sqrt(N)
ci_importance_upper <- importance_approx + 1.96 * s_importance / sqrt(N)
mean(importance_approx)
#0.3054301
print(c(ci_importance_lower,ci_importance_upper))
#[1]0.3015001 0.3076036

#1f)

#choose a new integral
new_v <-  rnorm(N, 15, 7)
val <- P(new_v)^2*dgamma(new_v, 0.8, 1/6)/dnorm(new_v, 15, 7)
print(mean(val))
# 0.4543129
# 0.3084502 (get from b)
var  <-  0.4543129- 0.3084502^2
#0.3591714

var(P(rgamma(N, 0.8, 1/6))) # 0.3527757


######################################################################################################
#2a)
Y_i <- c(162, 267, 271, 185, 111, 61, 27, 8, 3, 1)
I <- c(0:9)

# MLE of lambda
M <- sum(Y_i*I) / sum(Y_i)

# Compute values from distribution
n <- sum(Y_i)
M_val <- n * dpois(I, M, log = FALSE)

# Plot distribution with predicted counts and actual counts.
ggplot() +
  geom_line(
    aes(x = I, y = m_val), colour = "red") +
  geom_point(
    aes(x = I, y = Y_i), colour = "black") +
  
  labs(title = " Predicted values and Actual values",
       x = "i",
       y = "Density")




#2b,2c,2d,2e on the pdf
#2f)

# Initialize parameters
N <- 10000
p <- rep(0.5, N)
lambda_1 <- rep(2, N)
lambda_2 <- rep(2, N)
n <- length(Y_i)

p_bin <- rep(0, n)
Z_i <- rep(0, n)

for (i in 2:N) {
  
  for (j in 1:10){
    p_bin[j] <- p[i-1] * dpois(j-1, lambda_1[i-1]) / ((p[i-1] * dpois(j-1, lambda_1[i-1]) + (1-p[i-1]) * dpois(j-1, lambda_2[i-1])))
    Z_i[j] <- rbinom(1, Y_i[j], p_bin[j])
  }
  
  p[i] <- rbeta(1, (1 + n * mean(Z_i)), 1 + 1 + n * mean(Y_i - Z_i))
  
  for (j in 1:10){
    p_bin[j] <- p[i] * dpois(j-1, lambda_1[i-1]) / ((p[i] * dpois(j-1, lambda_1[i-1]) + (1-p[i]) * dpois(j-1, lambda_2[i-1])))
    Z_i[j] <- rbinom(1, Y_i[j], p_bin[j])
  }
  
  lambda_1[i] <- rgamma(1, (1 + sum(I * Z_i)), (1 + n * mean(Z_i)))
  
  for (j in 1:10){
    p_bin[j] <- p[i] * dpois(j-1, lambda_1[i]) / ((p[i] * dpois(j-1, lambda_1[i]) + (1-p[i]) * dpois(j-1, lambda_2[i-1])))
    Z_i[j] <- rbinom(1, Y_i[j], p_bin[j])
  }
  
  lambda_2[i] <- rgamma(1, (1 + sum(I * (Y_i - Z_i))), (1 + n * mean(Y_i - Z_i)))
}


#plot
ggplot() +
  geom_histogram(aes(x = p)
  )

ggplot() +
  geom_histogram(aes(x = lambda_1 )
  )

ggplot() +
  geom_histogram(aes(x = lambda_2 )
  )
