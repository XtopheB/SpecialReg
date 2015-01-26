# 26/01/2015  : Fonctions validées 


library(np)
library(foreign)
library(AER)
library(Formula)
#library(sem)
library(ks)  # for kernel estimation ...




wintrim.stata <- function(x,  trimtype = "TRIM" , trimlevel = 0.05)
  # function applying trimming to vector x as in STATA program
  # Note that the level of triming /winsoring is already in percent 
  # Quantile is computed for abs(x) and not on x  (Only upper quantile) !!
  
  
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
    # Replacing low values with low bound as in STATA 
    x.trim <- replace(x.trim, x <= -1*x.bound , -1*x.bound)    
  }
  # return(list(x.trim, x.bounds))
  print(paste("Number of trimed out/winsorized: ",length(which(abs(x) >= x.bound)), "obs."))
  return(x.trim)
}


# Specialred.fitN : Uses wintrim.stata 
# requires package np

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
    hetv <- as.matrix(hetv)
    uhat2 <- uhat^2
    summary(uhat2)
    est.hetero <- lm(uhat2~endo+ exo + iv+hetv) 
    summary(est.hetero)   # OK with STATA   !!!
    xbetahat <- est.hetero$fitted.value
    summary(xbetahat)
    uhat <- uhat /sqrt(abs(xbetahat))
  }
  print(" Summary of uhat")
  print(" ---------------")
  print(summary(uhat))   # OK with STATA   !!!
  
  
  ## Step 2: Nonparametric estimation  of the density of f
  
  
  #  normal-reference rule-of-thumb (Siverman) (always computed to define the object)
  # Stata bw for homo = 0.07154983, Hetero = .23179664 
  # bw.sil <- npudensbw(uhat, bws=0.23179664 , ckertype="epanechnikov", bandwidth.compute=FALSE) 
  
  if(udensmethod == "CV"){
    bw.dens <-npudensbw(~uhat, ckertype="epanechnikov", bandwidth.compute =TRUE )
  }
  if(udensmethod == "FIXED") {
    # fhat with user-defined bandwidth
    bw.dens <- npudensbw(uhat, bws=ubw , ckertype="epanechnikov", bandwidth.compute=FALSE) 
    #     bw.fixed <- bw.sil  # to start with an existing object
    #     bw.fixed$bw <- ubw
    #     bw.fixed$bandwidth <- ubw
    #     dens.np <- npudens(bws=bw.fixed,ckertype="epanechnikov", bandwidth.compute = FALSE)
  }
  # Density estimation 
  dens.np <- npudens(bws=bw.dens,ckertype="epanechnikov", bandwidth.compute = FALSE)
  
  #   Computing fhat 
  #   print(" ---Computing fhat ---")
  #   print(" ---------------")
  fhat <-dens.np$dens
  
  print(" --- Summary of fhat ---")
  print(" ---------------")
  print(summary(fhat))
  print(sd(fhat))     # Still OK with stata
  
  print(paste("---length of fhat:",length(fhat), "obs"))
  print(paste("--- Nb of fhat = 0 : ",length(which(abs(fhat) < 0.0000001)),"obs"))
  print(" ")
  
  ## Step 3: Definition of T1
  print("Step3 -- Computation of T1 ") 
  print("  ")
  
  vpos <- as.numeric(v >= 0)
  T1 <- ( as.numeric(y==1) - vpos)/ fhat   
  
  # Correction if Hetero option
  if(hetero=="HETERO"){
    T1 <- T1 * sqrt(abs(xbetahat))   # From step 1bis
  }
  print(" --- Summary of T1 ---")
  print(" ---------------")
  print(summary(T1))
  print(" -------------------")
  
  # Apply trimming / Winsorization
  T1.trim <- wintrim.stata(T1, trimtype = trimtype, trimlevel = trimlevel)
  
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
  
  
  ## END
}

