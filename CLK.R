####Appendix S4:R code for the CLK method and the simulations


###The simulated species are very similar to Phillips and Elith (2013). 
###The model fitting for LI and LK methods adopted the R code provided by Phillips and Elith (2013)

library(nloptr)
library(msm)

setwd("H:/R/SDM/PO")

# The following function NEW provide the new CLK method proposed in the paper
# It takes one extra arguments:
# 'fun"-- the type of link function used to fit the data
# 3 functions were used in the paper, logit, log and complementary loglog functions
# 
# Optional arguments include:
#   "quadratic" -- true if quadratic terms are to be used
#   "nlines"  -- the number of replicates in the experiment
#   "xl" -- true if the x axis labels are to be included in the plot
#  


NEW<-function(prevalence,presence,background,fun) 
{
  presence <- cbind(1, as.matrix(presence))
  background <- cbind(1, as.matrix(background))
  x0<-c(rep(1,npar))
  fn<-function(x) -mean(log(fun(presence,x)))
  heq <- function(x) mean(fun(background,x))-prevalence # heq == 0
  auglag(c(rep(0.5,npar)), fn, heq = heq,localsolver="lbfgs")$par # with COBYLA  choose c(1,1) as the starting value
}  


# Simulated response functions used
constant <- function(x) 0.3*(x/x)
linear<- function(x) 0.05 + 0.2*x
exponential<-function(x) exp(-4+4*x)
Quadratic <- function(x) .5 - 1.333*((x-0.5)*(x-0.5))
Gaussian <- function(x) 0.75*exp(-(4*x-2)^2)

semilogistic <- function(x) 8*plogis(-4+2*x)
logistic_2 <- function(x) plogis(-4+8*x)
logistic_1 <- function(x) plogis(-4+2*x)

fun0<-function (x,beta) {plogis(x%*%beta)}  ##logit function
fun<-function (x,beta) {exp(x%*%beta)}  ##exponential function
fun1<-function(x,beta) {1-exp(-exp(x%*%beta))}###completentary loglog function


###############paramters that needs to be changed #########
choice<-c("constant","linear","exponential","semilogistic","logistic_2","logistic_1","Quadratic","Gaussian")
quadratic<-F
#########################################################

npar<-if (quadratic) npar<-3 else npar<-2

xl<-T
np <-1000
nb <-20000
nt<-20000
nlines<-100


xfun<-function (num) runif(num)

par(mfrow=c(2,2))

for (ii in 1:length(choice))
{ 
  pp<-get(noquote(choice[ii]))
  
  if (ii>6) quadratic<-T
  
  print (pp)
  
  #############can change testx to different distribution #####
  test<-sort(xfun(nt))
  
  
  err<-matrix(0,nlines,9)
  
  ############## draw the P picture $$$$$$$$$$$
  truth <- pp(test)
  plot(test, truth, col=1, type="l", xlab="x", ylab="P(y=1|x)", ylim=c(0,1), lwd=3, xaxt='n')
  axis(1, labels=xl)
  
  test <- if (quadratic) {cbind(test, test*test) } else test
  
  NEWbeta<-matrix(0,nlines,npar)
  NEWbeta0<-matrix(0,nlines,npar)
  NEWbeta1<-matrix(0,nlines,npar)
  
  for (i in 1:nlines) 
  {
    print (ii)
    print (i)
    ####### generate presence data, x 
    x<-xfun(nt)
    
    p <- pp(x)
    # make some presence-absence data
    y <- rbinom(nt, 1, p)
    if (sum(y) < np) stop("Prevalence is low, need more than 50000 samples") 
    # keep np presence points
    px <- x[y==1][1:np]
    
    
    presence <- if (quadratic) data.frame(x=px, xx=px*px) else data.frame(x=px)
    
    # generate random background data
    bx<-xfun(nb)
    background <- if (quadratic) data.frame(x=bx, xx=bx*bx) else data.frame(x=bx)
    
    prevalence <- mean(p)
    prevalencemin <- min(0.99, prevalence+0.1)
    prevalencemax <- max(0.01, prevalence-0.1)
    
    
    NEWbeta0[i,]<-NEW(prevalence, presence, background,fun0)
    NEWbeta[i,]<-NEW(prevalence, presence, background,fun)
    NEWbeta1[i,]<-NEW(prevalence, presence, background,fun1)
    
    makeline1(NEWbeta[i,],test,fun, col="yellow", lwd=3)
    makeline1(NEWbeta0[i,],test,fun0,col="red", lwd=3)
    makeline1(NEWbeta1[i,],test,fun1, col="blue",lwd=3)
    
    ## The model fitting for LI and LK method adopted the function LI and LI provided by Phillips and Elith (2013)
    ## with the extra argument fun
    ##  LKbeta[i,]<-LK(presence, background,fun0)
    ##  LIbeta[i,]<-LI(presence, background,fun0)  
    
  } ###end of i in nline loop
  
  new0max<-NEW(prevalencemax, presence, background,fun0)
  new0min<-NEW(prevalencemin, presence, background,fun0)
  
  makeline1(new0max,test,fun0, col="red",lty=5, lwd=3)
  makeline1(new0min,test,fun0, col="red",lty=5, lwd=3)
  
  newmax<-NEW(prevalencemax, presence, background,fun)
  newmin<-NEW(prevalencemin, presence, background,fun)
  
  makeline1(newmax,test,fun, col="yellow",lty=5, lwd=3)
  makeline1(newmin,test,fun, col="yellow",lty=5, lwd=3)
  
  new1max<-NEW(prevalencemax, presence, background,fun1)
  new1min<-NEW(prevalencemin, presence, background,fun1)
  
  makeline1(new1max,test,fun1, col="blue",lty=5, lwd=3)
  makeline1(new1min,test,fun1, col="blue",lty=5, lwd=3)
  
  lines(test,truth, col="black", type="l", lwd=3,main=choice[ii])
} ### end of ii in nfunction loops  


# plot a logistic model with parameters given by "beta" applied to predictors in "dat"
makeline1 <- function(beta, dat,fun, ...) 
{
  preds <- apply(as.matrix(dat), 1, function(x) fun(c(1,x),beta))
  lines(as.matrix(dat)[,1], preds, ...)
}



