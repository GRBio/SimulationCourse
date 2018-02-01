rm(list=ls())

##-- Function to analytical estimate the parameters of a lognormal distribution given their moments
param.log.normal <- function(m,v){
  par.loc <- log(m^4/(v+m^2))/2
  par.scale <- sqrt(log(v/m^2 + 1))
  return(c(par.loc,par.scale))
}

##-- Parameters for log-Normal function used in the script
par1 <- param.log.normal(m=1,v=4^2)
par2 <- param.log.normal(m=2,v=4^2)

# Are they skewed?
curve(dlnorm(x,par1[1],par1[2]),0,8,col=2,lwd=2)
curve(dlnorm(x,par2[2],par2[2]),0,8,col=2,lwd=2)

#-------------------------------------
#
# Simple example with separate stages
#
#-------------------------------------

######################################
# A: Aims
######################################

# Comparing the performance of the two-sample t test 
# under skew data between unbalanced groups

######################################
# D: Data generating mechanism
######################################
##-- Parameters
N <- 30                                                                              # Total sample sizes
Nsim <- 1000

##-- Set a seed for reproducibility
set.seed(2008)

##-- Factor with 2 categories
x <- sample(1:2,N*Nsim,rep=TRUE,prob=c(2/3,1/3))                                     # Factor with a P(X=1)=2/3 & P(X=2)=1/3

##-- Response
y <- numeric(N*Nsim)
# y[x==1] <- round(rlnorm(sum(x==1),meanlog = -log(2)/2, sdlog = sqrt(log(2))),3)      # E[Y|X=1] = 1 & V[Y|X=1] = 1
# y[x==2] <- round(rlnorm(sum(x==2),meanlog = log(16/5)/2, sdlog = sqrt(log(5/4))),3)  # E[Y|X=2] = 2 & V[Y|X=2] = 1
y[x==1] <- round(rlnorm(sum(x==1),meanlog = par1[1], sdlog = par1[2]),3)               # E[Y|X=1] = 1 & V[Y|X=1] = 4
y[x==2] <- round(rlnorm(sum(x==2),meanlog = par2[1], sdlog = par2[2]),3)               # E[Y|X=2] = 2 & V[Y|X=2] = 4

##-- Store in a matrix
X <- matrix(x,ncol=Nsim)
Y <- matrix(y,ncol=Nsim)

######################################
# M: T-test
######################################
TTEST <- list()
for (i in 1:Nsim) TTEST[[i]] <- t.test(Y[,i]~X[,i],var.equal=TRUE) 


######################################
# E: Estimands --> IC95
######################################
IC95 <- matrix(nrow=Nsim,ncol=2)
colnames(IC95) <- c('LL','UL')
for (i in 1:Nsim) IC95[i,1:2] <- TTEST[[i]]$conf.int


######################################
# P: Performance --> Coverage
######################################
MD <- -1                                 # True Mean Difference [MD = E(Y|X=1) - E(Y|X=2) = -1]
for (i in 1:Nsim) 
included <- (MD>IC95[,'LL'] & MD<IC95[,'UL'])

##-- Point estimate of the coverage
point.coverage <- prop.table(table(included))[2]
point.coverage

##-- CI95% of the coverage
prop.test(x=sum(included),n=length(included))







#-------------------------------------
#
# Example using Montecarlo package
#
#-------------------------------------
library(MonteCarlo)

##-- T.test function
# N: Total sample size
# SD: Standard deviation for ecah sample
# P: Proportion assigned to arm 2
ttest <- function(N,SD,P){
  
  ##-- Generate data
  par1 <- param.log.normal(m=1,v=SD^2)
  par2 <- param.log.normal(m=2,v=SD^2)
  
  # Factor
  x <- sample(1:2,N,rep=TRUE,prob=c(1-P,P))                               
  
  # Response
  y <- numeric(N)
  y[x==1] <- round(rlnorm(sum(x==1),meanlog = par1[1], sdlog = par1[2]),3)
  y[x==2] <- round(rlnorm(sum(x==2),meanlog = par2[1], sdlog = par2[2]),3)  
  

  ##-- Calculate test statistic
  MD <- (-1)
  tt <- t.test(y~x,var.equal = TRUE)
  included <- MD>tt$conf.int[1] & MD<tt$conf.int[2]
  
  # return result:
  return(list("Coverage"=included))
}


######################################
# AIO: All In ONE
######################################
##-- Parameters grid
N_grid <- c(30,50,100)
SD_grid <- 1:5
P_grid <- c(0.3,0.5)
param_list <- list("N"=N_grid,"SD"=SD_grid,"P"=P_grid)


##-- Run simulation
set.seed(2008)
MC_res <- MonteCarlo(func=ttest, nrep=1000, param_list=param_list)
MakeTable(output=MC_res, rows=c('SD','N'), cols='P', digits=2, include_meta=FALSE)


