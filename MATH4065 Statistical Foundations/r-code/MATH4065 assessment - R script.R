

##### Q1 a) Read the data file into R and investigate whether or not there is any graphical evidence for relationships between 𝑌 and the predictor variables. 

data = read.csv("~/Documents/uon/MATH4065/assessment/ExperimentData.csv")

par(mfrow=c(1,2))
plot(Y~X1, data = data, main="Y vs X1")
plot(Y~X2, data = data, main="Y vs X2")
par(mfrow=c(1,1))




##### Q1 b) Using R, fit a linear model (model 1) in which 𝑌 is the response variable, 𝑋1 is the predictor variable and there is also an intercept term. Find the least-squares estimates of the model parameters. Is there evidence that 𝑋1 can explain the variation in 𝑌, or not?

# use R build-in function lm(...)
model1 = lm(Y ~ X1, data = data)

# see the summary of model1
summary(model1)






##### Q1 c) Using R, fit a new model (model 2) in which 𝑌 is the response variable, both 𝑋1 and 𝑋2 are predictor variables, and there is an intercept term. Is there evidence that model 2 fits the data significantly better than model 1, or not?

# use R build-in function lm(...)
model2 = lm(Y ~ X1 + X2, data = data)

# see the summary of model2
summary(model2)




##### Q1 d)
residu = model2$residuals



## i) find qq plot for residuals to see the residuals like normal
qqnorm(residu)

## ii) residuals against fitted values
plot(fitted(model2), residu, xlab="fitted values", ylab="Residuals", 
     main="Residuals vs fitted values")
abline(h=0, col="red", lty=2)


## iii) residuals against predictors
par(mfrow=c(1, 2))
plot(data$X1, residu, xlab="X1", ylab="Residuals", 
     main="Residuals vs X1")
abline(h=0, col="red", lty=2)

plot(data$X2, residu, xlab="X2", ylab="Residuals", 
     main="Residuals vs X2")
abline(h=0, col="red", lty=2)


par(mfrow=c(1, 1))



## 2) get residuals sum of square, R^2 and adjusted R^2

RSS = sum(residu^2)


summary(model2)$r.squared
# 0.94795

summary(model2)$adj.r.squared
# 0.9442322






#####################################################################################
# remove all declared variables
rm(list=ls())

#####################################################################################

##### Q2 b)
## given data of X
data = c(
  0.337, 0.507, 0.250, 3.131, 6.908, 
  0.195, -3.097, 1.002, 2.011, 1.710
)

## next define the loglikelihood function
loglikelihood_function = function(sigma){
  
  n = length(data)
  sum_of_x = sum(data)
  
  # first term
  term1 = -1 * n * log(sigma)
  
  # second term
  term2 = (1/sigma) * sum_of_x
  
  # third term
  term3 = -2 * sum(log(1 + exp(data / sigma)))
  
  return (term1 + term2 + term3)
}

## to ensure global max will be reached, first plot the graph to see if there are any local max
## 0 is excluded because we know that sigma > 0
sigmas = seq(1e-9, 500, by=.1)
logliks = sapply(sigmas, FUN=function(sigma){
  loglikelihood_function(sigma)
})

par(mfrow=c(2, 1))
plot(sigmas, logliks, type="l", main="Plot for looking at global max")
plot(sigmas, logliks, type="l", main="Plot for looking at global max for sigma between 0 and 3", xlim=c(0, 3))
par(mfrow=c(1, 1))
## from the plot, we can see that there are only one maximum, 
## and it is located between 0 and 3
## so we are safe to use optimize in the interval (1e-9, 3)

## because we are optimizing one single parameter
## so we need to use optimize

optimize_result = optimize(loglikelihood_function, lower = 1e-9, upper = 3, maximum = TRUE)

mle_sigma = optimize_result$maximum
mle_sigma    # 1.437023






##### Q2 c)
max_loglik = optimize_result$objective


## plot the loglikelihood function, with sigma between 0 to 5
sigmas = seq(1e-9, 3, length=1000)
logliks = sapply(sigmas, FUN=function(sigma){
  loglikelihood_function(sigma)
})


# because we know that loglik value goes to -Inf when sigma tends to 0
# and there are only one Global Max
# so we can just look at the range of y to (-50, 0), which will also contain the location of MLE
plot(sigmas, logliks, type="l", ylim=c(-50, 0),
     xlab=expression(sigma),
     ylab=expression(paste("\u2113", "(", sigma, ")")),
     main="log-likelihood function plot")

points(mle_sigma, max_loglik, col="red")
text(mle_sigma, max_loglik, col="red", labels=parse(text = paste0('hat(sigma) == ', round(mle_sigma, 2))), pos = 3)

lines(c(mle_sigma, mle_sigma), c(-100, max_loglik), lty = 2, col="red")

legend(x = "topright",
       inset = .05,
       legend = c(expression(paste("\u2113", "(", sigma, ")")), 
                  expression(paste("location of MLE: ", hat(sigma)))),
       col=c("black", "red"),
       lty=c(1, 2))



#####################################################################################
# remove all declared variables
rm(list=ls())
#####################################################################################


##### Q3 a) Write a function in R whose input is a positive integer 𝑇 and whose output is a single realisation of
#{𝑋1, 𝑋2, ... , 𝑋𝑇}.

# because T is a reserved keyword in R
# so I use another name
single_realization = function(int.T){

  # by the question requirement, input can only be positive integer
  # so check if the input is an integer AND positive, or else exit with message
  if(int.T <= 0 | (int.T != round(int.T))){
    stop(paste("input is not an positive integer: ", int.T))
  } else {
    
    # confirmed positive integer, can go on the logic
    
    # if the input is exactly 1, then just return a vector that containing X1, i.e., 0
    if(int.T == 1){
      return (c(0))
    } else {
      # execute when the input is positive integer greater than 1
      
      # first store the initial value
      process = as.numeric(0)
      
      
      for(i in 2:int.T){
        # note that the length of process is always >= 1
        # so getting the last element is always valid
        x_current = tail(process, 1)    
        
        # epsilon can either 1 or -1, and is equally likely
        epsilon = sample(x = c(-1, 1), size = 1, prob = c(0.5, 0.5))
        
        # the given relation
        x_next = 2 - x_current + epsilon
        
        # append the new item x_next at the end of the process vector
        process = c(process, x_next)
      }
     
      return (process)
       
    }
  }
}


##### Q3 b) Produce a plot showing a single realisation of {𝑋1, 𝑋2, ... , 𝑋10}, with the axes suitably labelled.
n = 10
process = single_realization(n)

# plot the line graph as requested
plot(1:n, process, type="l", lty=2, col="red",
     xlim=c(1, n),
     ylim=c(min(process) - 1, max(process) + 1),    # make the graph a bit wider along y-axis, so an extra +-1 in ylim
     xlab="n",
     ylab=expression(X[n]),
     main=expression(paste("Plot of single realization of {", X[1], ", ", X[2], ", ..., ", X[10], "}")))
points(1:n, process)
axis(1, at = 1:n)
axis(2, at = seq(min(process)-1, max(process)+1, by=1))

legend(x = "topleft", inset = .05,
       legend=c("transition", expression(X[n])),
       lty=c(2, NA),
       pch=c(NA, 1), 
       col=c("red", "black"))



##### Q3 c) Use your function to find numerical estimates of the variance of 𝑋100 and 𝑃[𝑋21 = 0]. You will not be given marks for using any other method.


### first part: find numerical estimates of the variance of X100
### i.e. simulate a couple of times, with input = 100, and then find the last element
### finally calculate the variance

# number of simulation
N = 10000

# we want to find input n = 100
n = 100


# we will then obtain 10000 rows and 100 columns matrix

# N for number of simulation
# n is the last time step
sim = function(N, n){
  m = as.numeric()
  
  for(i in 1:N){
    realization = single_realization(n)
    m = rbind(m, realization)
  }
  
  colnames(m) = 1:n
  rownames(m) = 1:N
  
  m
}

# after that, perform the simulation as per request
sim_result = sim(N, n)

# find the variance of X100, i.e. the last columns
vector_X100 = sim_result[, 100]
var(vector_X100)
# 99.68786


# to find P(X21 = 0), we first obtain a vector of X21, then sum if X21 == 0, and finally divided by N

vector_X21 = sim_result[, 21]

# find the ratio
sum(vector_X21 == 0) / N
# 0.1787

