#Function to produce P-values for the Ljung-Box test for different lags
#where an ARMA(p,q) model has been fitted.
#Note that k must be > p+q (See Lecture 9 slides)
#Number of degrees of freedom for the test = k-p-q

#Arguments for the function "LB_test"
#resid = residuals from a fitted ARMA(p,q) model.

#max.k = the maximum value of k at which we perform the test
#Note that the minimum k is set at p+q+1 (corresponding to a test with one degree
#of freedom)

#p = Order of the AR part of the model
#q = Order of the MA part of the model 

#The function returns a table with one column showing the number of degrees 
#of freedom for the test and the other the associated P-value.

LB_test<-function(resid,max.k,p,q){
 lb_result<-list()
 df<-list()
 p_value<-list()
  for(i in (p+q+1):max.k){
   lb_result[[i]]<-Box.test(resid,lag=i,type=c("Ljung-Box"),fitdf=(p+q))
   df[[i]]<-lb_result[[i]]$parameter
   p_value[[i]]<-lb_result[[i]]$p.value
  }
 df<-as.vector(unlist(df))
 p_value<-as.vector(unlist(p_value))
 test_output<-data.frame(df,p_value)
 names(test_output)<-c("deg_freedom","LB_p_value")
 return(test_output)
 }