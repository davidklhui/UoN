


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
plot(xs, ys, type="l", ylim=c(0,0.5))

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
            
      # store this value, later to find 2.5; 97.5 percentile
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




L = 100
hs2 = seq(0.05, 1.0, length=L)
bs2 = rep(NA, L)
for(i in 1:L){
  h = hs2[i]
  print(paste("i=", i, ", h=", h))
  out = bootstrap(J=200, N=100, B=1000, resample_size=100, h=h)  

  # count the ratio that covers theta_true
  bs2[i] = sum(out[,3]) / length(out[,3])
  
  print(bs2[i])
}

# running very slow, I save this result first in case I need this later
df = data.frame(h = hs2, b = bs2)
write.csv(df, "MATH4007 Computational Statistics/L100_J200_B1000_R100.csv", row.names=FALSE)

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

