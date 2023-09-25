

### Question 1

## Q1 a) Describe the rejection algorithm using a uniform proposal distribution
##> see the handwritten work pdf

## Q1 b) Write an R function to implement your rejection method. The function should take
## as inputs the parameters r and R, and also store the number of attempts needed to 
## successfully produce one accepted sample of X


# define the pdf of X
f = function(x, r, R){
  (1 / (2 * pi)) * (1 + (r / R) * cos(x) )
}

# implement the rejection method, using Y ~ U(0,2*pi) uniform proposal distribution
# parameters: r, R are parameter of X
rejection_algorithm = function(r, R){
  
  # step0: define local variables
  # from the hand written work, we know that M = (1 + r/R)
  M = 1 + r / R
  
  # and we also need to know the umber of attempts
  n = 0
  
  repeat {
    
    # increase the counter
    n = n + 1
    
    # step1: simulate Y ~ U(0,2*pi) and U ~ U(0,1)
    Y = runif(1, 0, 2*pi)
    U = runif(1)
    
    # avoid manually type the pdf of g, just use dunif
    ratio = f(Y, r, R) /  ( M * dunif(Y, 0, 2*pi) )
    
    # step2: determine of U < f(Y) / Mg(Y), if not, repeat the same code block
    if(U < ratio){
      X = Y
      
      break
    }
  }
  
  return (data.frame(num_attempts = n, X = X))
}

# just try one to see if that works
rejection_algorithm(r = 0.9, R = 1)


## Q1 c) For the case r = 0.9, R = 1, use your function to sample 10000 values of X
## Plot a histogram of your samples and overlay the true pdf, to check the sampler
# is working

# simulate 10000
sim_xs = function(r, R, N = 10000){
  
  # initialize a "empty" data frame
  df = data.frame()
  
  # write a for loop, for each iteration, find the newly calculated dataframe using Q1b
  # and then "append" to the data frame df
  for(i in 1:N){
    
    local_df = rejection_algorithm(r, R)
    df = rbind(df, local_df)
  }
  
  df
}

r0 = 0.9
R0 = 1

history = sim_xs(r = r0, R = R0, N = 10000)
# just take a look
dim(history)
head(history)

# plot the histogram
hist(history$X, freq=FALSE, xlab="x", ylab="density", 
     main = "histogram of samples with fitted true pdf",
     breaks=50,
     ylim=c(0, 0.35))

xs = seq(0, 2 * pi, length=1000)
fxs = sapply(xs, FUN=function(x){f(x, r = r0, R = R0)})
lines(xs, fxs, col="red")

legend(x="topright", 
       legend=c("true pdf"),
       col=c("red"),
       lty=c(1))


## Q1 d) How does the theoretical efficiency (in terms of expected
## number of attempts needed to successfully produce one accepted
## sample of X) depend on the parameters r and R?

# sample mean attempts
sample_mean_attempts = mean(history$num_attempts)
sample_mean_attempts

# because the expected number of attempts needed 
# to successfully produce one accepted sample of X  = M = 1 + r / R
# so for smaller r and larger R will produce smaller M
# i.e. less attempts to accept one sample of X


## Q1 e) For various values of r and R, use your R function
## to illustrate empirically your answer in (d)


## for numerical illustration,
# calculate the mean number of attempts for different combinations of r and R
rR_combinations = matrix(c(0.1, 1,
                           0.1, 10, 
                           0.1, 100, 
                           0.1, 500,
                           0.1, 1000,
                           0.01, 1,
                           0.05, 1,
                           0.1, 1,
                           0.5, 1, 
                           0.9, 1), nc=2, byrow=T)

# run the simulation and find the sample mean number of trials
m = apply(rR_combinations, 1, FUN=function(row){
  r = row[1]
  R = row[2]
  mean(sim_xs(r = r, R = R, N = 10000)$num_attempts)
})

# calculate the theoretical mean number of trials
theoretical_m = apply(rR_combinations, 1, FUN=function(row){
  r = row[1]
  R = row[2]
  1 + r / R
})

# combine those, rename the columns and see the result
result = cbind(rR_combinations, m, theoretical_m)
colnames(result) = c('r', 'R', 'sample_mean_trials', 'theoretical_mean_trials')
result



# Graphical Illustration, see f(x) and M*g(x)

# first note that the upper bound of f is 
# M * g = (1 + r / R) * dunif(x, 0, 2*pi)
# I will show some r and R
# note that this "does not depend upon x" because dunif(x, 0, 2*pi) = 1/(2*pi)
upper_f = function(x, r, R){
  (1 + r / R) * (dunif(x, 0, 2*pi))
}

# declares the x-axis, which is the [0, 2* pi]
xs = seq(0, 2 * pi, length=1000)

# define a helper function to plot the graph
graph_plotter = function(r, R){
  
  # find f
  fxs = sapply(xs, FUN=function(x){f(x, r = r, R = R)})
  
  # plot f against x
  plot(xs, fxs, type="l", main=title(paste("r=", r, ", R = ", R)),
       ylim=c(0, max(fxs)))
  
  # find the upper bound of f
  upper_f1 = sapply(xs, FUN=function(x){upper_f(x, r, R)})
  
  # plot the lines upper_f against x
  lines(xs, upper_f1, col="red", lty=2)
  
  # add the legend
  legend(x="bottomright", 
         legend=c("f(x)", "Mg(x)"),
         lty=c(1,2),
         col=c("black", "red"))
  
  
}



par(mfrow=c(2, 3))

# fixed r, increasing R
# so r / R is decreasing
# expect f(x) will be closer to the bound Mg(x)
graph_plotter(r = 0.1, R = 1)
graph_plotter(r = 0.1, R = 10)
graph_plotter(r = 0.1, R = 100)

# fixed R, increasing r (r<R)
# so r / R is increasing
# expect f(x) will be more away from the bound Mg(x)
graph_plotter(r = 0.1, R = 1)
graph_plotter(r = 0.5, R = 1)
graph_plotter(r = 0.9, R = 1)

# choose r < R such that r close to R, R = r + 0.01
# so r / R close to 1
# expect the shape are highly similar

par(mfrow=c(1,3))
graph_plotter(r = 0.1, R = 0.11)
graph_plotter(r = 1, R = 1.01)
graph_plotter(r = 10, R = 10.01)

par(mfrow=c(1,1))    







##############################################################################################
## remove all defined variables
rm(list=ls())


##############################################################################################











### Question 2



# define the data in matrix format
data = matrix(c(1:10, 1.3, 1.9, 2.4, 2.5, 2.4, 2.5, 2.6, 2.7, 2.7, 2.7), nc=2, byrow=F)
data

# define the constants for us to use
sigma_square = 0.01
lambda = 0.01

## Q2 b) Give details of, and write R code to implement,
## the 2-d steepest ascent algorithm to find the mode of the posterior distribution
## using log

# define the function, which takes log of the posterior distribution (ignored normalizing constant)
# f = log(posterior) ignored proportional constant
f = function(par){
  b1 = par[1]
  b2 = par[2]
  
  term1 = (b1 / sigma_square) * sum(apply(data, 1, FUN=function(row){
    x = row[1]
    y = row[2]
    
    y * exp(1 - exp(-1 * x / b2))
  }))
  
  term2 = (b1^2 / 2 / sigma_square) * sum(apply(data, 1, FUN=function(row){
    x = row[1]
    y = row[2]
    
    exp(2 - 2* exp(-1 * x / b2))
  }))
  
  term3 = lambda * (b1 + b2)
  term1 - term2 - term3
}


# function of making contour plot
f2 = function(b1, b2){
  f(c(b1, b2))
}

# plotted a contour plot
b1s = seq(0, 500, by=1)
b2s = seq(0, 500, by=1)
z = outer(b1s, b2s, FUN=Vectorize(f2))

# raw contour plot
contour(b1s, b2s, z)

# see only for large z
contour(b1s, b2s, z, zlim=c(max(z) * 0.9, max(z)))

# take a closer look to the contour plot
contour(b1s, b2s, z, zlim=c(max(z) * 0.99, max(z)), xlim=c(0.95,1.05), ylim=c(1,3.5))
abline(h=2, v=1)
# can expect that optimal point around b1 = 1, b2 = 2

# define the gradient
dfdb1 = function(par){
  b1 = par[1]
  b2 = par[2]
  
  term1 = (1/sigma_square) * sum(apply(data, 1, FUN=function(row){
    x = row[1]
    y = row[2]
    
    y * exp(1 - exp(-1 * x / b2))
  }))
  
  term2 = (b1 / sigma_square) * sum(apply(data, 1, FUN=function(row){
    x = row[1]
    y = row[2]
    
    exp(2 - 2 * exp(-1 * x / b2))
  }))
  
  term1 - term2 - lambda
}

dfdb2 = function(par){
  b1 = par[1]
  b2 = par[2]
  
  term1 = (b1/sigma_square) * sum(apply(data, 1, FUN=function(row){
    x = row[1]
    y = row[2]
    
    y * exp(1 - exp(-1 * x / b2)) * (-1 * exp(-1 * x / b2)) * (x / b2^2)
  }))
  
  term2 = (b1^2/ (2 * sigma_square)) * sum(apply(data, 1, FUN=function(row){
    x = row[1]
    y = row[2]
    
    exp(2 - 2 * exp(-1 * x / b2)) * (-2 * exp(-1 * x / b2)) * (x / b2^2)
  }))
  
  term1 - term2 - lambda
}

# using adaptive step size (not fixed delta)
steepest_ascent = function(f, dfdb1, dfdb2, x0, tol = 1e-6, max_steps = 1000){
  
  # declare initial diff, any large value is fine
  diff = 9999
  history = data.frame(n = 0, b1=x0[1], b2=x0[2], f = f(x0), dfdb1 = dfdb1(x0), dfdb2 = dfdb2(x0), diff = diff, converged = (diff < tol))
  
  # define the total runs
  n = 0
  
  # initial guess, which will be changed in subsequent runs
  x = x0
  
  # runs only when diff >= tolerence and n <= max allowed runs
  while(diff >= tol & n <= max_steps ){
    n = n + 1
    delta = 1e-3  # adaptive step size
    
    f_val = f(x)
    dfdb1_val = dfdb1(x)
    dfdb2_val = dfdb2(x)
    
    x_new = NA
    f_new_val = NA
    trial = 0
    while(trial < 100){
      trial = trial + 1
      
      # this is the x(delta) from lecture notes and my handwritten work pdf
      x_new = x + delta * c(dfdb1_val, dfdb2_val)
      
      # calculate f using this new value
      f_new_val = f(x_new)
      
      # check if the new value is less than old value (not maximizing)
      # so divide the delta and re-calculate the above
      # it is confident that we can find a "good" f_new within 100 trials
      if(f_new_val < f_val){
        delta = delta / 2
      } else {
        break;
      }
    }
    
    dfdb1_new_val = dfdb1(x_new)
    dfdb2_new_val = dfdb2(x_new)
    
    diff = abs(f_new_val - f_val)
    
    # put everything into a single dataframe
    df = data.frame(n = n, b1 = x_new[1], b2=x_new[2], f = f_new_val, dfdb1 = dfdb1(x_new), dfdb2 = dfdb2(x_new), diff = diff, converged = (diff < tol))
    #    print(df)
    history = rbind(history, df)
    
    x = x_new
  }
  
  
  history  
  
}

# initial guess, using my value from the contour plot
x0 = c(b1 = 1, b2 = 2)

history = steepest_ascent(f, dfdb1, dfdb2, x0 = x0, tol = 1e-9, max_steps = 100000)
tail(history, 100)
#n       b1       b2        f       dfdb1        dfdb2         diff converged
#641 0.987296 1.819924 2891.782 -0.04893951 -0.009577987 9.958967e-11      TRUE



## Q2 c) Give full details of Laplace’s method to compute at a particular point b2
## see handwritten work




## Q2 d) Write a function in R to compute π(b2|x,y) at a particular point b2 
## using Laplace’s method derived in (c).

# because our f above is working only up to proportionality
# as the question (other than 2b) did not mention that I can use up to proportionality, I will just go for FULL pdf
# however, we do not really need to explicit write down the posterior because we can make use of R built-in function
# using dnorm and dexp is enough
full_posterior = function(b1, b2){
  # joint distribution of the data, which is independent normal
  t1 = prod(apply(data, 1, FUN=function(row){
    x = row[1]
    y = row[2]
    
    d = dnorm(y, mean=b1 * exp(1 - exp(-1 * x / b2)), sd = sqrt(sigma_square))
  }))
  
  # prior distribution
  t2 = dexp(b1, lambda) * dexp(b2, lambda)
  
  # multiply them
  t1 * t2
}

# first we need the expression of b1_hat
b1_hat = function(b2){
  
  term1 = sum(apply(data, 1, FUN=function(row){
    x = row[1]
    y = row[2]
    
    y * exp(1 - exp(-1 * x / b2))
  }))
  
  term2 = lambda * sigma_square
  
  term3 = sum(apply(data, 1, FUN=function(row){
    x = row[1]
    y = row[2]
    
    exp(2 - 2 * exp(-1 * x / b2))
  }))
  
  (term1 - term2) / term3
}

# next we need the hessian of g
H = function(b2){
  (1/sigma_square) * sum(apply(data, 1, FUN=function(row){
    x = row[1]
    y = row[2]
    
    exp(2 - 2 * exp(-1 * x / b2))
  }))
}

# the we can put all things together and find the marginal posterior distribution
marginal_posterior_b2 = function(b2){
  
  b1_hat_val = b1_hat(b2)
  H_val = H(b2)                     # our H_val is a real number
  post = full_posterior(b1_hat_val, b2)  # use our FULL posterior
  
  post * (sqrt(2 * pi) / sqrt(H_val)) * pnorm(b1_hat_val * sqrt(H_val))
}

## Q2 e) Plot π(b2|x,y) using your function from (d).

# we know that b2 > 0 as b2 ~ exp(lambda)
# define the axis for plotting
# Because I have tried the range of b2 [0, 300], but clearly the density is close to 0 for b2 > 3
# so below I will just use [0, 5]
b2s = seq(0, 5, length=1000)

marginal_posterior_b2_val = sapply(b2s, FUN=function(b2){
  marginal_posterior_b2(b2)
})
plot(b2s, marginal_posterior_b2_val, type="l", main="Laplace Approx",
     xlab="b2", ylab="density")





## Q2 f) Write a function in R to perform the Golden-ratio method to find the mode of
## using your R function from part (d) as the function to optimize.

# my function to-be-optimized is marginal_posterior_b2
# from the graph, we can limit our range to search for max be [1,3]

# f - function to be maximize
# [xmin, xmax] - search region
golden_ratio_method_max = function(f, xmin, xmax, tol=1e-4){
  x1 = xmin
  x3 = xmax
  x2 = as.numeric(0)
  x4 = as.numeric(0)
  n = 0
  r= (sqrt(5)-1)/2
  while(abs(x3-x1)>= tol){
    n = n+1
    x2 = r * x1 + (1-r) * x3
    x4 = (1-r) * x1 + r * x3
    
    if(f(x2) < f(x4)){
      # change the interval to (x2,x3)
      x1 = x2
    } else {
      # case when f(x2) > f(x4)
      # change the interval to (x1,x4)
      x3 = x4
    }
  }  
  
  return (list(xmin=x1, xmax=x3, iteration=n, diff=(x3-x1), y=f(x1)))
}


## Q2 g) Hence, find the mode of p(b2|y) to an accuracy of 1 decimal place.

# I guess the question has typo, should be finding mode of π(b2|x,y)
gr_method_result = golden_ratio_method_max(marginal_posterior_b2, 1, 3, tol = 1e-12)

# because I used a relatively small tolerance, so for 1 d.p. numerical result must be enough
round(gr_method_result$xmin, 1)
# 1.8

