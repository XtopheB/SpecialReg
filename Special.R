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

wintrim <- function(x,  trimtype = "TRIM" , trimlevel = 0.05)
  # function applying trimming to vector x
  # Note that the level of triming /winsoring is already in percent 
{ 
  trim.U <- 1 - trimlevel/2
  trim.L <- trimlevel/2
  x.bounds <- quantile(x, probs = c(trim.L, trim.U), na.rm = TRUE)
  ##### ---- option na.rm  To be CONFIRMED -----------------------
  
  if(trimtype =="TRIM"){
    #  x.trim <- subset(x, x >= x.bounds[1] & x<= x.bounds[2] ) <-- Old version
    # Replacing low values by NAs
    x.trim <- replace(x, x <= x.bounds[1]| x >= x.bounds[2] , NA)    
    
  }
  if(trimtype =="WINSOR"){
    # Replacing low values by low bound
    x.trim <- replace(x, x <= x.bounds[1], x.bounds[1])    
    #Replacing high values by high bounds
    x.trim <- replace(x.trim, x.trim >= x.bounds[2], x.bounds[2])
  }
  # return(list(x.trim, x.bounds))
  print("number of points")
  return(x.trim)
}

# Clean version
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

est1 <- lm(v~endo+x+z)
summary(est1)
uhat <- v - est1$fitted.value   #Computing residuals

## Step 2: Nonparametric estimation  of the density of f

##--->  Branch 1: Choice of estimation method NP or Ordered choice

  # --> if NP estimator
    ## ---> Branch 2: Either  bw is chosed by the user, or CV or Silverman (default) 
    
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




specialreg2 <- function(D,V,endo, z, x,
                        trimtype = "TRIM",
                        trimlevel = 0.05, 
                        hetero = "", 
                        hetv = x, 
                        data = data.work
                        )
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
  
  # Step 1: OLS of V on X and Z
  
  est1 <- lm(v~endo+x+z)
  summary(est1)
    
  # Step 1 bis: Computing residuals depending on hetero option
  uhat <- v - est1$fitted.value   #Computing residuals
  
  if(hetero=="HETERO"){
    uhat2 <- uhat^2
    est.hetero <- lm(uhat2~endo+x+z+hetv) 
    xbetahat <- est.hetero$fitted.value
    uhat <- uhat /sqrt(abs(xbetahat))
  }
  
  ## Step 2: Nonparametric estimation  of the density of f
  
  ##--->  Branch 1: Choice of estimation method NP or Ordered choice
  
  # --> if NP estimator
  ## ---> Branch 2: Either  bw is chosed by the user, or CV or Silverman (default) 
  # TO UNCOMMENT  (Computing Time ) 
  #bw.cv <-npudensbw(~uhat, ckertype="epanechnikov")
  #summary(bw.cv)
  
  # normal-reference rule-of-thumb (Siverman) 
  bw.sil <-npudensbw(~uhat,bwmethod="normal-reference")
  summary(bw.sil)
  
  # fhat with selected bandwidth
  fhat <- fitted(npudens(bw.sil))
  
  ## Step 3: Definition of T
  print("Step3") 
  
  vpos <- as.numeric(v >= 0)
  T <- ( as.numeric(D) - vpos)/ fhat
  
  # Correction if Hetero option
  if(hetero=="HETERO"){
    T <- T * sqrt(abs(xbetahat))   # From step 1bis
  }
  # Apply trimming / Winsorization
  T.trim <- wintrim(T, trimtype= trimtype, trimlevel = trimlevel)
  print(summary(T.trim))
  print(dim(T.trim))
  
  ## Step 4: 2SLS
  print("Step4")
  
  # IV regression with package AER 
  # spe.ivreg <- ivreg(T~endo+x , ~z+x, data= data.work)
  spe.ivreg <- ivreg(T.trim~endo+x , ~z+x, data= subset(data.work, T.trim !=NA ))
  
  
  summary(spe.ivreg)
  
  # TODO implement manualy 2SLS to compare the results. 
  
  ## END
}




#Testing our functions 

## Load the data.
data.all<-read.dta("../Water/data/table for ssreg.dta")
nrow(data.all)
data.work <- subset(data.all, country == "FRANCE" | country == "CANADA" | country == "AUSTRALIA")
nrow(data.work)
attach(data.work)

data.work$i_can <- (country == "CANADA")
data.work$i_fra <- (country == "FRANCE")
data.work$i_aus <- (country == "AUSTRALIA")


# Here we prepare data with standard notations 
D <- i_tap
v <- prix_2008
endo <- isatis_health
z <- cbind(itap_2008 ,iconcernwatpol_2008)
x <- cbind(i_under18 , log_income ,i_town, i_car, b08_locenv_water, a2_age, i_can, i_fra)

#Test of  main special regression function 

#specialreg(D,V, endo, z, x)

specialreg2(D,V, endo, z, x, trimtype="TRIM" )




## ---- SAND BOX  ------
# Testing wintrim 

x <- rnorm(1000)
summary(x)
# a 5%: trim = 0.05
Montrim <- wintrim(x)   # <- defaut
quantile(x, probs=c(0.025, 0.975))
summary(Montrim)

# a 2.5%: trim = 0.025
Montrim2 <- wintrim(x, trimlevel = 0.025, trimtype="WINSOR")  
quantile(x, probs=c(0.0125, 0.9875))
summary(Montrim2)



