#  Program implementing the Special regressor method (Lewbel et al.) 
#  12/11/2014: speciareg function
#  13/11/2014: triming and winsoring function (wintrim)
#  16/11/2014 : Creation of specialreg2 that seems to work quite well
#  17/11/2014 : Creation of specialreg3 with Formula (not working)
#  18/11/2014 : Eureka: Creating ivreg.fit with new def of x and Z (exo added to z)
#             :  use new syntax of ivreg and compares to ivreg.fit

rm(list=ls())
setwd("D:/progs/Celine/Special")   

## libraries

library(np)
library(foreign)
library(AER)
library(Formula)
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




specialreg2 <- function(D,v,endo, z, exo,
                        trimtype = "TRIM",
                        trimlevel = 0.05, 
                        hetero = "", 
                        hetv = exo, 
                        data = data.work
                        )
  # Main function for special regressor 
{
  # Defining the type of vectors
  D <- as.factor(D)
  v <- as.numeric(v)
  endo <- as.matrix(endo)
  z <- as.matrix(z)
  exo <- as.matrix(exo)
  
  # Step 0 : Demean V
  v <- v-mean(v) 
  
  # Step 1: OLS of V on X and Z
  
  est1 <- lm(v~endo+exo+z)
  print(summary(est1))
    
  # Step 1 bis: Computing residuals depending on hetero option
  uhat <- v - est1$fitted.value   #Computing residuals
  
  if(hetero=="HETERO"){
    uhat2 <- uhat^2
    est.hetero <- lm(uhat2~endo+exo+z+hetv) 
    xbetahat <- est.hetero$fitted.value
    uhat <- uhat /sqrt(abs(xbetahat))
  }
  print(summary(uhat))
  
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
  print(summary(T))
  
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
  # spe.ivreg <- ivreg(T~endo+exo , ~z+exo, data= data.work)
  spe.ivreg <- ivreg(T.trim~endo+exo , ~z+exo, data= subset(data.work, T.trim !=NA ))
  
  summary(spe.ivreg)
  return(spe.ivreg)
  # TODO implement manualy 2SLS to compare the results. 
  
  ## END
}


specialreg.fit <- function(y,v, x, z,
                        trimtype = "TRIM",
                        trimlevel = 0.05, 
                        hetero = "", 
                        hetv = x, 
                        data = data.work
)
  # Main function for special regressor 
{
  # Defining the type of vectors
  y <- as.factor(y)
  v <- as.numeric(v)
  x <- as.matrix(x)
  z <- as.matrix(z)
  
  # Step 0 : Demean V
  v <- v-mean(v) 
  
  # Step 1: OLS of V on X and Z
  
  est1 <- lm(v~x+z)
  print(summary(est1))
  
  # Step 1 bis: Computing residuals depending on hetero option
  uhat <- v - est1$fitted.value   #Computing residuals
  
  if(hetero=="HETERO"){
    uhat2 <- uhat^2
    est.hetero <- lm(uhat2~x+z+hetv) 
    xbetahat <- est.hetero$fitted.value
    uhat <- uhat /sqrt(abs(xbetahat))
  }
  print(summary(uhat))
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
  T <- ( as.numeric(y) - vpos)/ fhat
  print(summary(T))
  
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
  # spe.ivreg <- ivreg(T~x , ~z+x, data= data.work)
  #spe.ivreg <- ivreg(T.trim~x , ~z, data= subset(data.work, T.trim !=NA ))
  spe.ivreg <- ivreg(T.trim~x | z, data= subset(data.work, T.trim !=NA ))
  summary(spe.ivreg)
  return(spe.ivreg)
  # TODO implement manualy 2SLS to compare the results. 
  
  ## END
}


specialreg3 <- function(formula, data, subset, na.action=, 
                        trimtype = "TRIM",
                        trimlevel = 0.05, 
                        hetero = "", 
                        hetv = "",
                        ...
)
{
## Application of Formula as in Achim Zeileis, Yves Croissant (2010). 
## Extended Model Formulas in R: Multiple Parts and Multiple Responses.
## Journal of Statistical Software 34(1), 1-13. URL http://www.jstatsoft.org/v34/i01/.  
    
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  
  f <- Formula(formula)
  mf[[1]] <- as.name("model.frame")
  mf$formula <- f
  mf <- eval(mf, parent.frame())
  
  y <- model.response(mf)
  x <- model.matrix(f, data = mf, rhs = 1)
  z <- model.matrix(f, data = mf, rhs = 2)
  
  # Construction of the variables used in the latter function 
  #Decison variable : y
  y <- model.response(mf)
  # Explanatory variables
  x <- model.matrix(f, data = mf, rhs = 1)
#head(x)
  # Special regressor
  v <- x[,2]
#head(v)
  # Extracting special regressor from x
  x <- x[, -2]
#head(x)
  # Instruments
  z <- model.matrix(f, data = mf, rhs = 2)
  
## -------- START of the Specialreg function itself -------
  # Step 1: OLS of V on X and Z
  
  est1 <- lm(v~x+z)
  summary(est1)
  
  # Step 1 bis: Computing residuals depending on hetero option
  uhat <- v - est1$fitted.value   #Computing residuals
  
  if(hetero=="HETERO"){
    uhat2 <- uhat^2
    est.hetero <- lm(uhat2~x+z+hetv) 
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
  T <- ( as.numeric(y) - vpos)/ fhat
  
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
  #spe.ivreg <- ivreg(T.trim~x , ~z, data= subset(data.work, T.trim !=NA ))
  spe.ivreg <- ivreg(T.trim~x | z, data= subset(data.work, T.trim !=NA ))
  summary(spe.ivreg)
  return(spe.ivreg)
  # TODO implement manualy 2SLS to compare the results. 
  
  ## END
  }


#Testing our functions 

## Load the data.
data.all<-read.dta("../Water/data/table for ssreg.dta")
nrow(data.all)
data.work <- subset(data.all, country == "FRANCE" | country == "CANADA" | country == "AUSTRALIA")
nrow(data.work)

data.work$i_can <- (data.work$country == "CANADA")
data.work$i_fra <- (data.work$country == "FRANCE")
data.work$i_aus <- (data.work$country == "AUSTRALIA")

attach(data.work)
# Here we prepare data with standard notations 
D <- i_tap
special <- prix_2008
Myendo <- isatis_health
iv <- cbind(itap_2008 ,iconcernwatpol_2008)
exo <- cbind(i_under18 , log_income ,i_town, i_car, b08_locenv_water, a2_age, i_can, i_fra)

#Test of  main special regression function 

#specialreg(D,V, endo, z, x)

specialreg2(D,special, Myendo, iv, exo, trimtype="TRIM" )

# Attention, nombre et ordre des variables changé : 
#  : On mets tous les X
# On ajoute les X exo aux instruments

NewX <- cbind(isatis_health, i_under18 , log_income ,i_town, i_car, b08_locenv_water, a2_age, i_can, i_fra)
Newiv <- cbind(iv, exo)
specialreg.fit(D,special, NewX, Newiv,  trimtype="TRIM" )





# specialreg3(D~prix_2008 + x |z, data = data.work )
specialreg3(i_tap~prix_2008 + i_under18  |itap_2008, data = data.work )



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



