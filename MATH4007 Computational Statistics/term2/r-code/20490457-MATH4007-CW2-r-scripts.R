

###########################################
### Q1 #####################

### using the given instruction of f(x; mu0, mu1, p)
f = function(x){
  
  mu0 = 0
  mu1 = 3
  p = 0.3
  
  (1-p) * dnorm(x, mean=mu0, sd=1) + 
    p * dnorm(x, mean=mu1, sd=1)
  
}

### although this is not required, just 
### see how the mixture pdf looks like
xs = seq(-6, 9, length=1000)
ys = sapply(xs, FUN=function(x){f(x)})
plot(xs, ys, type="l", ylim=c(0,0.5), 
     main="Mixture pdf")

lines(xs, dnorm(xs, 0, 1), col="red")
lines(xs, dnorm(xs, 3, 1), col="blue")

legend(x="topright", legend=c("mixture", "N(0,1)", "N(3,1)"),
       col=c("black", "red", "blue"),
       lty=1)



## true value of theta
theta_true = f(0)   
# [1] 0.2805892

lines(c(-10,0), c(theta_true, theta_true), lty=2)
lines(c(0,0), c(0, theta_true), lty=2)

#################################################

### prepare a method to simulate f directly
### without the use of any simulation techniques in Ch2
### just directly use the definition of mixtures of 2 distributions
sim_f = function(N){
  
  U = runif(N, 0, 1)
  sapply(U, function(u){
    if(u<=0.3){
      rnorm(1, mean=3, sd=1)
    } else {
      rnorm(1, mean=0, sd=1)
    }
  })
}


### Check if this method applicable
fs = sim_f(100000)
hist(fs, freq=FALSE, 
     main="Histogram of the simulation", ylim=c(0,0.3))

points(fs, rep(0, length(fs)), col="green")
xs = seq(-6, 9, length=1000)
ys = sapply(xs, FUN=function(x){f(x)})
lines(xs, ys, col="red")
legend(x="topright", legend=c("Theoretical f"),
       col="red", lty=1)

### from this plot, we can see that 
### our sim_f works


#################################################


### implement Kernel Density Estimation
### because I need to get the estimated density at x=0
### but I am not quite sure how to get this from built-in density(...) command
### so implement from scratch is more convenience
kde = function(data_points, at, h){
  
  n = length(data_points)
  
  Ks = sapply(data_points, FUN=function(xi){
    # using standard normal kernel
    dnorm((at - xi) / h)
  })
  
  sum(Ks) / (n*h)
}

### try to see if this work
### first get a set of points
### then compare it to build-in density() method
f100 = sim_f(100)

buildin_kde = density(f100)
buildin_kde$bw
plot(buildin_kde)

xs = seq(-10,10,by=.01)
ys = sapply(xs, FUN=function(x){kde(f100, at=x, h=buildin_kde$bw)})
lines(xs, ys, col="red")

### we can see a perfect match
### so my implemented kde method works, give same effect of the build-in density()



#################################################

### implement nonparametric bootstrap method

### J - total replications for B times bootstrap replicates
### N - simulation size from f
### B - bootstrap replicates
### resample_size - resample size for each bootstrap
### h - bandwidth
bootstrap = function(J, N, B, resample_size, h){
  
  # stores a matrix of 95% percentile intervals; and each determine if it contains theta_true
  out = matrix(rep(NA,3*J), nc=3)
  
  for(j in 1:J){
    if(j%%10==0) print(paste("j=",j))
    
    #  Simulate n = 100 data points from f
    data_points = sim_f(N)
    
    # store all estimates of theta for this jth iteration
    # later will use this vector to get the 95% interval
    tmp = rep(NA, B)
    
    
    for(b in 1:B){
      
      
      # Bootstrapping, sample points from the data_points
      sample_data = sample(data_points, size = resample_size, replace = TRUE)
      
      # use my function to get do kde, and get the value at 0
      theta_est = kde(sample_data, at=0, h)
      
      # store this value, later to find required percentiles
      tmp[b] = theta_est
      
    }
    
    # calculate the interval's upper bound and lower bound
    LB = quantile(tmp, 0.025)
    UB = quantile(tmp, 0.975)
    
    # Store this interval.
    out[j, 1] = LB
    out[j, 2] = UB
    
    # check if the interval contains true theta
    out[j, 3] = (LB<=theta_true) & (theta_true <=UB)
  }
  
  out
}


###############
#### for any fixed h
# try h = 0.5
out_at_fixed_h = bootstrap(J=200, N=100, B=1000, resample_size=100, h=0.5)

estimated_coverage = sum(out_at_fixed_h[,3]) / length(out_at_fixed_h[,3])
estimated_coverage
# around 0.84

#######################
#### for various h

L = 100
hs2 = seq(0.05, 1.0, length=L)
bs2 = rep(NA, L)
for(i in 1:L){
  h = hs2[i]
  print(paste("i=", i, ", h=", h))
  out = bootstrap(J=200, N=100, B=1000, resample_size=100, h=h)  
  
  # count the ratio that covers theta_true
  # or simply use mean
  bs2[i] = sum(out[,3]) / length(out[,3])
  
  ## can print the estimate of the coverage variability
  ## if I can assume normality, then the LB=mean-1.96se; UB=mean+1.96*se
  ## and so (UB-LB)/4 can be approximately the se
  ## as expected, higher value for small h; smaller value for large h

  # print(paste("estimate variability at ", h, "=",sqrt(mean((out[,2]-(out[,1])/4)^2))))
  
  print(bs2[i])
}

# running very slow, I save this result first in case I need this later
df = data.frame(h = hs2, b = bs2)
write.csv(df, "MATH4007 Computational Statistics/result.csv", row.names=FALSE)


### in case required, load this back
#df = read.csv("MATH4007 Computational Statistics/L100_J200_B1000_R100.csv", header=TRUE)
#hs2 = df$h
#bs2 = df$b

# plot estimated Coverage vs bandwidth
plot(hs2, bs2,  
     type="l", main="estimated Coverage vs bandwidth", 
     xlim=c(0.05,1), ylim=c(0,1), 
     xlab="bandwidth (h)",
     ylab="estimated coverage")

target_h = hs2[which(bs2>=0.95)]
length(target_h)

# mark the reference line 0.95
abline(h=0.95, col="red")


# overlay a smoothed curve using loess
lo = loess(bs2~hs2)
hs3 = seq(0, 1, length=1000)
lines(hs3, predict(lo, hs3), col="blue", lwd=2)








####### 
# using the same simulated datapoints, draw some same KDE plots 
h1 = 0.075
h2 = 0.2419192
h3 = 0.5
h4 = 0.7
h5 = 1

sim_f_100 = sim_f(100)


f1 = sapply(xs, FUN=function(x){kde(sim_f_100, at=x, h=h1)})
f2 = sapply(xs, FUN=function(x){kde(sim_f_100, at=x, h=h2)})
f3 = sapply(xs, FUN=function(x){kde(sim_f_100, at=x, h=h3)})
f4 = sapply(xs, FUN=function(x){kde(sim_f_100, at=x, h=h4)})
f5 = sapply(xs, FUN=function(x){kde(sim_f_100, at=x, h=h5)})

plot(xs, f1, xlim=c(-5,5),
     main="KDE using different bandwidth",
     xlab="x", ylab="y",
     ylim=c(0, max(f1, f2, f3, f4, f5)), type="l")
lines(xs, f2, col=2)
lines(xs, f3, col=3)
lines(xs, f4, col=4)
lines(xs, f5, col=5)


points(sim_f_100, rep(0, length(sim_f_100)), col=2)

abline(h=theta_true, col=6, lty=2, lwd=2)

legend(x="topleft", col=c(1:5, 2), 
       legend=c("h=0.075", "h=0.2419192", "h=0.5", "h=0.7", "h=1", "data points"),
       lty=c(1,1,1,1,1,NA),
       pch=c(NA,NA,NA,NA,NA,1))


### we can see that if h is small, then 
### the y-value is highly possible to above the true theta at x = 0
### but it is also subjected to much higher then the true theta


### if h is large, then it is almost impossible to reach the true theta

### in between the best h from my simulation is h = 0.2419 can balance the two




####################################################################
rm(list=ls())
####################################################################


#######################################################3
##### Q2


##########################################################
##### Q2c Write R code to implement your Gibbs sampler

data1 = c(5,4,7,2,3)
data2 = c(2,1,5,5,5)
data3 = c(9,7,6,7,10)
data4 = c(8,4,13,11,5)

## N - number of samples
# initial_mu - any initial starting point, I choose 0.1
gibbs_sampler = function(N, initial_mu = 0.1){
  
  ## initialize mu
  ## _cur represents current value
  ## _new represent new value
  mu_cur = initial_mu
  
  ## calculate the sum for each group, and +1 to match the parameter
  S1 = sum(data1) + 1
  S2 = sum(data2) + 1
  S3 = sum(data3) + 1
  S4 = sum(data4) + 1
  
  ## store lambda1, lambda2, lambda3, lambda4, mu
  chain = matrix(rep(NA, 5*N), nc=5)
  
  for(i in 1:N){
    
    ## sample lambda1 from Gamma(sum(y1j)+1), mu_cur + 5)
    ## perform the same sampling for each lambda
    lambda1_new = rgamma(1, shape = S1, rate = mu_cur + 5)
    lambda2_new = rgamma(1, shape = S2, rate = mu_cur + 5)
    lambda3_new = rgamma(1, shape = S3, rate = mu_cur + 5)
    lambda4_new = rgamma(1, shape = S4, rate = mu_cur + 5)
    
    lambdas_new = c(lambda1_new, lambda2_new, 
                    lambda3_new, lambda4_new)
    
    ## sample new mu from lambdai_cur, for i = 1...4
    mu_new = rgamma(1, shape = 4, rate = sum(lambdas_new))
    
    
    chain[i,] = c(lambdas_new, mu_new)
    
    ## update mu
    mu_cur = mu_new
  }
  
  chain
  
}



##########################################################
##### Q2d Run your sampler and produce suitable evidence to show
##### you can be confident it is performing reasonably

# may also try different initial mu
chain = gibbs_sampler(N=100000, initial_mu = 0.1)

diagnostic_plot = function(chain, title){
  par(mfrow=c(3,1))
  for(j in 1:5){
    hist(chain, main=paste("histogram (",title,")"))
    ts.plot(chain, main="TS plot")
    acf(chain, main="ACF plot")
  }
  par(mfrow=c(1,1))
  
}

# lambda1
diagnostic_plot(chain[,1], "lambda1")

# lambda2
diagnostic_plot(chain[,2], "lambda2")

# lambda3
diagnostic_plot(chain[,3], "lambda3")

# lambda4
diagnostic_plot(chain[,4], "lambda4")

# mu
diagnostic_plot(chain[,5], "mu")

### 1) all TS plot looks reach stationary distribution
### 2) ACF close to 0 rapidly

### means, the Gibbs sampler performing well


##### Q2e use the output of your sampler to produce suitable plots
##### and numerical summaries in order to compare lambda1,2,3,4
##### and hence comment briefly on any differences between 
##### night admission rates at the various A&E departments

##### doesn't look contain burn-in period, 
##### no need to discard the first few hundred iterations

mean(chain[,1])
var(chain[,1])

mean(chain[,2])
var(chain[,2])

mean(chain[,3])
var(chain[,3])

mean(chain[,4])
var(chain[,4])

### we can see that A&E dept 4 has the largest mean and var
### whilst A&E dept 2 has the lowest mean and var

par(mfrow=c(2,2))
hist(chain[,1], main="A&E dept 1")
hist(chain[,2], main="A&E dept 2")
hist(chain[,3], main="A&E dept 3")
hist(chain[,4], main="A&E dept 4")

par(mfrow=c(1,1))
### the histogram looks similar, but with different mean


### use density() method to find the KDE
### and plot those together
kde1 = density(chain[,1])
kde2 = density(chain[,2])
kde3 = density(chain[,3])
kde4 = density(chain[,4])


plot(kde1, main="KDE plot of different lambda", ylim=c(0,0.5), xlim=c(0,16))
lines(kde2, col="red")
lines(kde3, col="blue")
lines(kde4, col="green")

legend(x="topright", inset=.05,
       legend=c("dept1", "dept2", "dept3", "dept4"),
       col=c("black", "red", "blue", "green"),
       lty=1)

### from the density estimation plot, we can see that 
### dept 1 and 2 are similar
### dept 3 and 4 are similar





##### Q2f Consider making predictions for a new night in A&E dept 4, denoted by y*
##### Explain how you can use the output of your sampler to simulate samples from the 
##### predictive distribution


### since we have samples obtained in Q2d, the marginal posterior distribution of lambda4 | y
### we can use those lambda4|y to get sample of y|lambda4 using rpoi



##### Q2g implement this to obtain samples of y* from pi(y*|y)
##### where y* is the quantity defined in (f)
##### Use your samples of y* to estimate the probability
##### that the number of admissions to A&E department 4 on this
##### new night is at least 10

### N here is the size of the samples of the predictive distribution 
predictive_dist = function(N, posterior_marginal_distribution){
  
  ### get N random samples from the posterior marginal distribution
  ### I think allow replace or not does not matter really matter in this case
  lambda_samples = sample(posterior_marginal_distribution, size=N, replace=TRUE)
  
  ### from lambda -> simulate Po(lambda)
  sapply(lambda_samples, FUN=function(lam){
    rpois(1, lambda = lam)
  })
  
}


### although there are no obvious burn-in period
### the first few value may still depends on the initial position of mu
### so for safety for prediction, just ignore the first 1000 samples
### get these sample, and visualize it
y_predict = predictive_dist(N = 10000, chain[1001:nrow(chain),4])
hist(y_predict, freq=FALSE)

# draw a line to denote y*=10
abline(v=10, col="red", lty=2)


### from the y_predict, find the proportion at least 10
sum(y_predict >= 10) / length(y_predict)

### answer is around 0.308












