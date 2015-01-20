# 20/01/2015  : Creation d'un fichier unique avec toutes les fonctions  à tester
#             : Récupere les fonction des programmes SpecialReg et TestFunction.old




rm(list=ls())
#setwd("D:/progs/Celine/Special")   
setwd("c:/Chris/progs/Celine/Special")   

## libraries

library(np)
library(foreign)
library(AER)
library(Formula)
#library(sem)


## Fonctions de trimming  

wintrim <- function(x,  trimtype = "TRIM" , trimlevel = 0.05)
  # function applying trimming to vector x
  # Note that the level of triming /winsoring is already in percent 
{ 
  trim.U <- 1 - trimlevel/2
  trim.L <- trimlevel/2
  x.bounds <- quantile(x, probs = c(trim.L, trim.U))
  ##### ---- option na.rm  To be CONFIRMED -----------------------
  
  if(trimtype =="TRIM"){
    #  x.trim <- subset(x, x >= x.bounds[1] & x<= x.bounds[2] ) <-- Old version
    # Replacing low values by NAs
    x.trim <- replace(x, x <= x.bounds[1] |  x >= x.bounds[2] , NA)    
    
  }
  if(trimtype =="WINSOR"){
    # Replacing low values by low bound
    x.trim <- replace(x, x <= x.bounds[1], x.bounds[1])    
    #Replacing high values by high bounds
    x.trim <- replace(x.trim, x.trim >= x.bounds[2], x.bounds[2])
  }
  # return(list(x.trim, x.bounds))
  #print("number of points trimed out")
  return(x.trim)
}


wintrim2 <- function(x,  trimtype = "TRIM" , trimlevel = 0.05)
  # function applying trimming to vector x
  # Note that the level of triming /winsoring is already in percent
  # Version that REMOVE points (instead of imputing NA)
{ 
  trim.U <- 1 - trimlevel/2
  trim.L <- trimlevel/2
  x.bounds <- quantile(x, probs = c(trim.L, trim.U), na.rm = TRUE)
  ##### ---- option na.rm  To be CONFIRMED -----------------------
  
  if(trimtype =="TRIM"){
    x.trim <- subset(x, x >= x.bounds[1] & x<= x.bounds[2] ) 
    
  }
  if(trimtype =="WINSOR"){
    # Replacing low values by low bound
    x.trim <- replace(x, x <= x.bounds[1], x.bounds[1])    
    #Replacing high values by high bounds
    x.trim <- replace(x.trim, x.trim >= x.bounds[2], x.bounds[2])
  }
  # return(list(x.trim, x.bounds))
  #print("number of points trimed out")
  return(x.trim)
}


wintrim.stata <- function(x,  trimtype = "TRIM" , trimlevel = 0.05)
  # function applying trimming to vector x as in STATA program
  # Note that the level of triming /winsoring is already in percent 
{ 
  
  x.bound <- quantile(abs(x), 1-trimlevel)
  ##### ---- option na.rm  To be CONFIRMED -----------------------
  
  if(trimtype =="TRIM"){
    #  x.trim <- subset(x, x >= x.bounds[1] & x<= x.bounds[2] ) <-- Old version
    # Replacing low values with NAs
    x.trim <- replace(x, abs(x)> x.bound , NA)    
    
  }
  if(trimtype =="WINSOR"){
    #Replacing high values by high bounds
    x.trim <- replace(x, abs(x) > x.bound, x.bound)
    # Replacing low valueswith low bound as in STATA 
    x.trim <- replace(x.trim, x <= -1*x.bound , -1*x.bound)    
  }
  # return(list(x.trim, x.bounds))
  print(paste("Number of trimed out/winsorized: ",length(which(abs(x) >= x.bound)), "obs."))
  return(x.trim)
}


wintrim.test <- function(x,  trimtype = "TRIM" , trimlevel = 0.05)
  # function applying trimming to vector x
  # Note that the level of triming /winsoring is already in percent 
{ 
  # Points to be trimmed-out 
  to.trim <-  which(abs(x) >= quantile(abs(x), probs = 1-trimlevel))
  if(trimtype =="TRIM"){
    # Replacing high values by NAs
    x.trim <- replace(x, to.trim , NA)    
  
  }
  if(trimtype =="WINSOR"){
    # Replacing high values with bound
    x.trim <- replace(x,  x > 0 & to.trim, quantile(x, probs = 1-trimlevel))
    x.trim <- replace(x.trim,  x < 0 & to.trim, -1*quantile(x, probs = 1-trimlevel))
    
  }
  # return(list(x.trim, x.bounds))
  print(paste("Number of trimed out.winsorized: ",length(to.trim), "obs."))
  return(x.trim)
}



## Fonctions  SpecialReg 



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
  
  ## Step 3: Definition of T1
  vpos <- as.numeric(v >= 0)
  T1 <- ( as.numeric(D) - vpos)/ fhat
  
  # Apply trimming / Winsorization
  
  
  ## Step 4: 2SLS
  
  # with AER 
  spe.ivreg <- ivreg(T1~endo+x , ~z+x, data= data.work)
  
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
  
  ## Step 3: Definition of T1
  print("Step3") 
  
  vpos <- as.numeric(v >= 0)
  T1 <- ( as.numeric(D) - vpos)/ fhat
  print(summary(T1))
  
  # Correction if Hetero option
  if(hetero=="HETERO"){
    T1 <- T1 * sqrt(abs(xbetahat))   # From step 1bis
  }
  # Apply trimming / Winsorization
  T1.trim <- wintrim(T1, trimtype= trimtype, trimlevel = trimlevel)
  print(summary(T1.trim))
  print(dim(T1.trim))
  
  ## Step 4: 2SLS
  print("Step4")
  # IV regression with package AER 
  # spe.ivreg <- ivreg(T1~endo+exo , ~z+exo, data= data.work)
  spe.ivreg <- ivreg(T1.trim~endo+exo , ~z+exo, data= subset(data.work, T1.trim !=NA ))
  
  summary(spe.ivreg)
  return(spe.ivreg)
  # TODO implement manualy 2SLS to compare the results. 
  
  ## END
}


#-------------------------------

specialreg.fit <- function(y,v, x, z,
                           trimtype = "TRIM",
                           trimlevel = 0.05, 
                           hetero = "", 
                           hetv = "", 
                           udensmethod = "CV",
                           ubw= "" , 
                           data = data.work,
                           ...
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
  print(" Summary of uhat")
  print(" ---------------")
  print(summary(uhat))
  
  
  ## Step 2: Nonparametric estimation  of the density of f
  
  # normal-reference rule-of-thumb (Siverman) (always computed to define the object)
  bw.sil <-npudensbw(~uhat,  bwmethod="normal-reference")
  dens.np <- npudens(bws=bw.sil,  bandwidth.compute =FALSE )
  # summary(dens.np)
  
  if(udensmethod == "CV"){
    bw.cv <-npudensbw(~uhat, ckertype="epanechnikov", bandwidth.compute =TRUE )
    dens.np<- npudens(bws=bw.cv, ckertype="epanechnikov",bandwidth.compute = FALSE)
    summary(dens.np)
  }
  if(udensmethod == "FIXED") {
    # fhat with user-defined bandwidth
    bw.fixed <- bw.sil  # to start with an existing object
    bw.fixed$bw <- ubw
    bw.fixed$bandwidth <- ubw
    dens.np <- npudens(bws=bw.fixed,ckertype="epanechnikov", bandwidth.compute = FALSE)
    summary(dens.np)
  }
  # Computing fhat 
  fhat <-predict(dens.np)
  
  print(" Summary of fhat")
  print(" ---------------")
  summary(fhat)
  print(paste("---length of fhat:",length(fhat), "obs"))
  print(paste("--- Nb of fhat = 0 : ",length(which(abs(fhat) < 0.0000001)),"obs"))
  print(" ")
  
  ## Step 3: Definition of T1
  print("Step3 -- Computation of T1 ") 
  print("  ")
  
  vpos <- as.numeric(v >= 0)
  T1 <- ( as.numeric(y) - vpos)/ fhat
  
  # Correction if Hetero option
  if(hetero=="HETERO"){
    T1 <- T1 * sqrt(abs(xbetahat))   # From step 1bis
  }
  
  # Apply trimming / Winsorization
  T1.trim <- wintrim.stata(T1, trimtype = trimtype, trimlevel = trimlevel)
  
  print(" Summary of T1")
  print(" ---------------")
  print("summary(T1)")
  
  ## Step 4: 2SLS
  print(" ")
  print("Step4:  2SLS")
  print("=============")
  
  # IV regression with package AER 
  # spe.ivreg <- ivreg(T1~x , ~z+x, data= data.work)
  #spe.ivreg <- ivreg(T1.trim~x , ~z, data= subset(data.work, T1.trim !=NA ))
  
  spe.ivreg <- ivreg(T1.trim~x | z, subset = T1.trim != "NA")
  print(summary(spe.ivreg))
  
  # print(paste("Formula in ivreg : ",spe.ivreg$call, "."))
  # 
  print(paste("Number of points : ",spe.ivreg$nobs,"obs."))
  # print(paste("Sample used in ivreg : ",length(spe.ivreg$residuals), "points."))
  
  
  return(spe.ivreg)
  # TODO implement manualy 2SLS to compare the results. 
  
  ## END
}



#### 


specialreg3 <- function(formula, data, subset, na.action=, 
                        trimtype = "TRIM",
                        trimlevel = 0.05, 
                        hetero = "", 
                        hetv = "",
                        ...
)
{
#   ## Application of Formula as in Achim Zeileis, Yves Croissant (2010). 
#   ## Extended Model Formulas in R: Multiple Parts and Multiple Responses.
#   ## Journal of Statistical Software 34(1), 1-13. URL http://www.jstatsoft.org/v34/i01/.  
#   
#   mf <- match.call(expand.dots = FALSE)
#   m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
#   mf <- mf[c(1, m)]
#   
#   f <- Formula(formula)
#   mf[[1]] <- as.name("model.frame")
#   mf$formula <- f
#   mf <- eval(mf, parent.frame())
#   
#   y <- model.response(mf)
#   x <- model.matrix(f, data = mf, rhs = 1)
#   z <- model.matrix(f, data = mf, rhs = 2)
#   
#   # Construction of the variables used in the latter function 
#   #Decison variable : y
#   y <- model.response(mf)
#   # Explanatory variables
#   x <- model.matrix(f, data = mf, rhs = 1)
#   #head(x)
#   # Special regressor
#   v <- x[,2]
#   #head(v)
#   # Extracting special regressor from x
#   x <- x[, -2]
#   #head(x)
#   # Instruments
#   z <- model.matrix(f, data = mf, rhs = 2)
#   
#   ## -------- START of the Specialreg function itself -------
#   # Step 1: OLS of V on X and Z
#   
#   est1 <- lm(v~x+z)
#   summary(est1)
#   
#   # Step 1 bis: Computing residuals depending on hetero option
#   uhat <- v - est1$fitted.value   #Computing residuals
#   
#   if(hetero=="HETERO"){
#     uhat2 <- uhat^2
#     est.hetero <- lm(uhat2~x+z+hetv) 
#     xbetahat <- est.hetero$fitted.value
#     uhat <- uhat /sqrt(abs(xbetahat))
#   }
#   
#   ## Step 2: Nonparametric estimation  of the density of f
#   
#   ##--->  Branch 1: Choice of estimation method NP or Ordered choice
#   
#   # --> if NP estimator
#   ## ---> Branch 2: Either  bw is chosed by the user, or CV or Silverman (default) 
#   # TO UNCOMMENT  (Computing Time ) 
#   #bw.cv <-npudensbw(~uhat, ckertype="epanechnikov")
#   #summary(bw.cv)
#   
#   # normal-reference rule-of-thumb (Siverman) 
#   bw.sil <-npudensbw(~uhat,bwmethod="normal-reference")
#   summary(bw.sil)
#   
#   # fhat with selected bandwidth
#   fhat <- fitted(npudens(bw.sil))
#   
#   ## Step 3: Definition of T1
#   print("Step3") 
#   
#   vpos <- as.numeric(v >= 0)
#   T1 <- (as.numeric(y) - vpos)/ fhat
#   
#   # Correction if Hetero option
#   if(hetero=="HETERO"){
#     T1 <- T1 * sqrt(abs(xbetahat))   # From step 1bis
#   }
#   # Apply trimming / Winsorization
#   T1.trim <- wintrim(T1, trimtype= trimtype, trimlevel = trimlevel)
#   print(summary(T1.trim))
#   print(dim(T1.trim))
#   
#   ## Step 4: 2SLS
#   print("Step4")
#   # IV regression with package AER 
#   # spe.ivreg <- ivreg(T1~endo+x , ~z+x, data= data.work)
#   #spe.ivreg <- ivreg(T1.trim~x , ~z, data= subset(data.work, T1.trim !=NA ))
#   spe.ivreg <- ivreg(T1.trim~x | z, data= subset(data.work, T1.trim !=NA ))
#   summary(spe.ivreg)
#   return(spe.ivreg)
#   # TODO implement manualy 2SLS to compare the results. 
#   
#   ## END
}



specialreg.fitN <- function(y,v, endo, exo,  iv,
                           trimtype = "TRIM",
                           trimlevel = 0.05, 
                           hetero = "", 
                           hetv = "", 
                           udensmethod = "CV",
                           ubw= "" , 
                           data = data.work,
                           ...
)
  # Main function for special regressor 
{
  # Defining the type of vectors
  y <- as.factor(y)
  v <- as.numeric(v)
  endo <- as.matrix(endo)
  exo <- as.matrix(exo)
  iv <- as.matrix(iv)
  # Step 0 : Demean V
  v <- v-mean(v) 
  
  # Step 1: OLS of V on exo and iv
  
  est1 <- lm(v~endo +exo+ + iv)
  print(summary(est1))
  
  # Step 1 bis: Computing residuals depending on hetero option
  uhat <- v - est1$fitted.value   #Computing residuals
  
  if(hetero=="HETERO"){
    uhat2 <- uhat^2
    est.hetero <- lm(uhat2~exo+endo + iv+hetv) 
    xbetahat <- est.hetero$fitted.value
    uhat <- uhat /sqrt(abs(xbetahat))
  }
  print(" Summary of uhat")
  print(" ---------------")
  print(summary(uhat))
  
  
  ## Step 2: Nonparametric estimation  of the density of f
  
  # normal-reference rule-of-thumb (Siverman) (always computed to define the object)
  bw.sil <-npudensbw(~uhat, ckertype="epanechnikov", bwmethod="normal-reference")
  print(" Summary of bw")
  print(" ---------------")
  summary(bw.sil)
  
  ### TO REMOVE ONE DAY !!!!!
  # One may change the bandwidth here to the one found in stata for the homo case !!: .07154983
  bw.sil$bw <- .07154983
  bw.sil$x <- .07154983
  
    
  dens.np <- npudens(bws=bw.sil, ckertype="epanechnikov", bandwidth.compute =FALSE )
  summary(dens.np)
  
  if(udensmethod == "CV"){
    bw.cv <-npudensbw(~uhat, ckertype="epanechnikov", bandwidth.compute =TRUE )
    dens.np<- npudens(bws=bw.cv, ckertype="epanechnikov",bandwidth.compute = FALSE)
    summary(dens.np)
  }
  if(udensmethod == "FIXED") {
    # fhat with user-defined bandwidth
    bw.fixed <- bw.sil  # to start with an existing object
    bw.fixed$bw <- ubw
    bw.fixed$bandwidth <- ubw
    dens.np <- npudens(bws=bw.fixed,ckertype="epanechnikov", bandwidth.compute = FALSE)
    summary(dens.np)
  }
  # Computing fhat 
  print(" ---Computing fhat ---")
  print(" ---------------")
  fhat <-predict(dens.np)
  
  print(" --- Summary of fhat ---")
  print(" ---------------")
  print(summary(fhat))
  print(" ")
  
  print(paste("---length of fhat:",length(fhat), "obs"))
  print(paste("--- Nb of fhat = 0 : ",length(which(abs(fhat) < 0.0000001)),"obs"))
  print(" ")
  
  ## Step 3: Definition of T1
  print("Step3 -- Computation of T1 ") 
  print("  ")
  
  vpos <- as.numeric(v >= 0)
  T1 <- ( as.numeric(y) - vpos)/ fhat
  
  # Correction if Hetero option
  if(hetero=="HETERO"){
    T1 <- T1 * sqrt(abs(xbetahat))   # From step 1bis
  }
  print(" --- Summary of T1 ---")
  print(" ---------------")
  print(summary(T1))
  print(" -------------------")
  
  # Apply trimming / Winsorization
  T1.trim <- wintrim.test(T1, trimtype = trimtype, trimlevel = trimlevel)
  
  print(" Summary of T1")
  print(" ---------------")
  print(summary(T1.trim))
  
  ## Step 4: 2SLS
  print(" ")
  print("Step4:  2SLS")
  print("=============")
  
  # IV regression with package AER 
  # spe.ivreg <- ivreg(T1~exo , ~iv+exo, data= data.work)
  #spe.ivreg <- ivreg(T1.trim~exo , ~iv, data= subset(data.work, T1.trim !=NA ))
  
  # spe.ivreg <- ivreg(T1.trim~ endo + exo | iv+exo, subset = T1.trim != "NA")
  
  
  spe.ivreg <- ivreg(T1.trim ~ endo + exo | iv +exo)
  
  
  print(summary(spe.ivreg))
  
  # print(paste("Formula in ivreg : ",spe.ivreg$call, "."))
  # 
  print(paste("Number of points : ",spe.ivreg$nobs,"obs."))
  # print(paste("Sample used in ivreg : ",length(spe.ivreg$residuals), "points."))
  
  
  return(spe.ivreg)
  # TODO implement manualy 2SLS to compare the results. 
  
  ## END
}



  
  
  
  
  
  