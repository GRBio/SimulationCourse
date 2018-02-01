rm(list=ls())

# Replace the holes with ellipses

#-------------------------------------
#
# Simple example with separate stages
#
#-------------------------------------

######################################
# A: Aims
######################################

# Comparing the performance of the two-sample t test 
# under ....

######################################
# D: Data generating mechanism
######################################
##-- Parameters
N <- ...               # Total sample sizes
Nsim <- ...            # Number of simulations Nsim=1000 --> 0.7%
P <- ...               # Probability of X=1

##-- Set a seed for reproducibility (any number)
set.seed(....)

##-- Factor with 2 categories
x <- sample(1:2,N*Nsim,rep=TRUE,prob=c(P,....))                                     # Factor with a P(X=1)=2/3 & P(X=2)=1/3

##-- Response
y <- numeric(N*Nsim)
y[x==1] <- ....               # E[Y|X=1] = 1 & V[Y|X=1] = ?
y[x==2] <- ....               # E[Y|X=2] = 2 & V[Y|X=2] = ?

# Some options

# Option 1
# y[x==1] <- round(rlnorm(sum(x==1),meanlog = -log(2)/2, sdlog = sqrt(log(2))),3)      # E[Y|X=1] = 1 & V[Y|X=1] = 1
# y[x==2] <- round(rlnorm(sum(x==2),meanlog = log(16/5)/2, sdlog = sqrt(log(5/4))),3)  # E[Y|X=2] = 2 & V[Y|X=1] = 1
# boxplot(y[x==1],y[x==2])

# Option 2
# y[x==1] <- round(rweibull(sum(x==1),shape = 1, scale = 1),3)                          # E[Y|X=1] = 1 & V[Y|X=1] = 1
# y[x==2] <- round(rweibull(sum(x==2),shape = .5, scale = 1),3)                         # E[Y|X=2] = 2 & V[Y|X=1] = 20
# boxplot(y[x==1],y[x==2])

# Option 3
# y[x==1] <- round(rgamma(sum(x==1),shape = 1, rate = 1),3)                              # E[Y|X=1] = 1 & V[Y|X=1] = 1
# y[x==2] <- round(rgamma(sum(x==2),shape = 4, rate = 2),3)                              # E[Y|X=2] = 2 & V[Y|X=1] = 1
# boxplot(y[x==1],y[x==2])

##-- Store in a matrix
X <- matrix(x,ncol=Nsim)
Y <- matrix(y,ncol=Nsim)

######################################
# M: T-test
######################################
TTEST <- list()                                                     # List for store the TTEST
for (i in 1:Nsim) TTEST[[i]] <- t.test(Y[,i]~X[,i],var.equal=....)  # For 


######################################
# E: Estimands --> IC95
######################################
IC95 <- matrix(nrow=Nsim,ncol=2)                          # Matrix for store the confidence intervals
colnames(IC95) <- c('LL','UL')                            # Names of the columns
for (i in 1:Nsim) IC95[i,1:2] <- TTEST[[i]]$...           # Extract the confidence interval


######################################
# P: Performance --> Coverage
######################################
MD <- -1                                                        # True Mean Difference [MD = E(Y|X=1) - E(Y|X=2) = -1]
for (i in 1:Nsim) included <- (MD>.... & MD< ....)

##-- Point estimate of the coverage
point.coverage <- prop.table(table(included))[2]
point.coverage

##-- CI95% of the coverage
prop.test(x=....,n=....)







#-------------------------------------
#
# Example using Montecarlo package
#
#-------------------------------------
library(MonteCarlo)

##-- T.test function
# N: Total sample size
# SD: Standard deviation for each sample
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


