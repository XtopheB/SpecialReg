#  Program implementing the Special regressor method (Lewbel et al.) 
#  12/11/2014: speciareg function
#  13/12/2014: triming and winsoring function (wintrim)


rm(list=ls())
setwd("D:/progs/Celine/Special")   

## libraries

library(np)
library(foreign)
library(AER)
#library(sem)

options(np.messages=FALSE)

specialreg <- function(D,V, endo, z, x)
  # Main function for special regressor 
{
# Defining the type of vectors
D <- as.factor(D)
v <- as.numeric(v)
endo <- as.matrix(endo)
z <- as.matrix(z)
x <- as.matrix(x)

# Step 0 : Demean V
v <- v-mean(v) 

# Step 1: OLS of V on X and Z; computing residuals.

est <- lm(v~endo+x+z)
summary(est)
uhat <- v-est$fitted.value

## Step 2: Nonparametric estimation  of the density of f

##--->  Branch 1: Choice of estimation method NP or Ordered choice

  # --> if NP estimator
    ## ---> Branch 2: Either  bw is chosed by hte user, or CV or Silverman (default) 
    
    bw.cv <-npudensbw(~uhat, ckertype="epanechnikov")
    summary(bw.cv)
    
    # normal-reference rule-of-thumb (Siverman) 
    bw.sil <-npudensbw(~uhat,bwmethod="normal-reference")
    summary(bw.sil)

  # fhat with selected bandwidth
  fhat <- fitted(npudens(bw.cv))

## Step 3: Definition of T
vpos <- as.numeric(v >= 0)
T <- ( as.numeric(D) - vpos)/ fhat

# Apply trimming / Winsorization


## Step 4: 2SLS

# with AER 
spe.ivreg <- ivreg(T~endo+x , ~z+x, data= data.work)

summary(spe.ivreg)

# TODO implement manualy 2SLS to compare the results. 

## END
}

wintrim <- function(x, trim = 0.05, winsor = "")
  # function applying trimming to vector x
{ 
  trim.U <- 1 - trim/2
  trim.L <- trim/2
  x.bounds <- quantile(x, probs = c(trim.L, trim.U))
  
  if(winsor ==""){
    x.trim <- subset(x, x >= x.bounds[1] & x<= x.bounds[2] )
  }
  if(winsor =="WINSOR"){
    x.trim <- replace(x, x <= x.bounds[1], x.bounds[1])
    x.trim <- replace(x.trim, x.trim >= x.bounds[2], x.bounds[2])
  }
  return(x.trim)
}


#Testing our functions 

## Load the data.
data.all<-read.dta("../Water/data/table for ssreg.dta")
nrow(data.all)
data.work <- subset(data.all, country == "FRANCE" | country == "CANADA" | country == "AUSTRALIA")
nrow(data.work)
attach(data.work)

# Here we prepare data with standard notations 
D <- i_tap
v <- prix_2008
endo <- isatis_health
z <- cbind(itap_2008 ,iconcernwatpol_2008)
x <- cbind(i_under18 , log_income ,i_town, i_car, b08_locenv_water, a2_age, i_can, i_fra)

#Test of  main fucntion 

specialreg(D,V, endo, z, x)

# Testing wintrim 
trim = 5
x <- rnorm(1000)
summary(x)
# a 5%: trim = 0.05
Montrim <- wintrim(x)   # <- defaut
quantile(x, probs=c(0.025, 0.975))
summary(Montrim)

# a 2.5%: trim = 0.025
Montrim2 <- wintrim(x, trim = 0.025, winsor="WINSOR")  
quantile(x, probs=c(0.0125, 0.9875))
summary(Montrim2)



