
##########################################################
##### Q2c Write R code to implement your Gibbs sampler

data1 = c(5,4,7,2,3)
data2 = c(2,1,5,5,5)
data3 = c(9,7,6,7,10)
data4 = c(8,4,13,11,5)

## N - number of samples
gibbs_sampler = function(N){
  
  ## initialize mu
  ## _cur represents current value
  ## _new represent new value
  mu_cur = 0.01
  
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

chain = gibbs_sampler(N=100000)

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
  lambda_samples = sample(posterior_marginal_distribution, size=N, replace=TRUE)
  
  ### from lambda -> simulate Po(lambda)
  sapply(lambda_samples, FUN=function(lam){
    rpois(1, lambda = lam)
  })
  
}


# get these sample, and visualize it
y_predict = predictive_dist(N = 10000, chain[,4])
hist(y_predict, freq=FALSE)

# draw a line to denote y*=10
abline(v=10, col="red", lty=2)


### from the y_predict, find the proportion at least 10
sum(y_predict >= 10) / length(y_predict)

### answer is around 0.308














