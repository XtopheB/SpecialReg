# Test direct de nos fonction avec les données...


rm(list=ls())
#setwd("D:/progs/Celine/Special")   
setwd("c:/Chris/progs/Celine/Special")   


## libraries

library(np)
library(foreign)
library(AER)
library(Formula)
#library(sem)

options(np.messages=FALSE)

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

wintrim.stata <- function(x,  trimtype = "TRIM" , trimlevel = 0.05)
  # function applying trimming to vector x as in STATA program
  # Note that the level of triming /winsoring is already in percent 
{ 
  
  x.bound <- quantile(abs(x), 1-trimlevel)
  ##### ---- option na.rm  To be CONFIRMED -----------------------
  
  if(trimtype =="TRIM"){
    #  x.trim <- subset(x, x >= x.bounds[1] & x<= x.bounds[2] ) <-- Old version
    # Replacing low values with NAs
    x.trim <- replace(x, abs(x)>= x.bound , NA)    
    
  }
  if(trimtype =="WINSOR"){
    #Replacing high values by high bounds
    x.trim <- replace(x, abs(x) >= x.bound, x.bound)
    # Replacing low valueswith low bound as in STATA 
    x.trim <- replace(x.trim, x <= -1*x.bound , -1*x.bound)    
  }
  # return(list(x.trim, x.bounds))
  print(paste("Number of trimed out/winsorized: ",length(which(abs(x) >= x.bound)), "obs."))
  return(x.trim)
}



# Visual Test
x <- v
x.trim <- wintrim.stata(x, trimtype="WINSOR", trimlevel=0.05)

plot(x, x, pch=1)
points(x.trim, x, col="red", pch=3)
summary(x)
summary(x.trim)

### TEST de la fonction Directement; 

## Load the dat directly after an estimation in STATA to have same sample 2771 obs. 
data.all<-read.dta("../Water/data/FinalFile.dta")
nrow(data.all)
data.work <- data.all

attach(data.work)
# Here we prepare data with standard notations 
D <- i_tap
special <- prix_2008
Myendo <- isatis_health
iv <- cbind(itap_2008 ,iconcernwatpol_2008)
exo <- cbind(a2_age, i_under18 , log_income ,i_town, i_car, b08_locenv_water,  i_can, i_fra)

# On ajoute les X exo aux instruments

NewX <- cbind(isatis_health, i_under18 , log_income ,i_town, i_car, b08_locenv_water, a2_age, i_can, i_fra)
Newiv <- cbind(iv, exo)

retour <- specialreg.fit(D,special, NewX, Newiv,  trimtype="TRIM", trimlevel=0.025 )





########################## DEBUT de la fonction ###############################

# Paramètres de la fonctions entrés directment 
y <- as.factor(D)
v <- as.numeric(Special)
x <- as.matrix(NewX)
z <- as.matrix(Newiv)
trimtype="TRIM"
trimlevel=0.025
hetero=""

udensmethod = "FIXED" 
ubw= 0.23 

###Debut 
  
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
  print(paste("---length of uhat",length(uhat), "obs"))
  
  ## Step 2: Nonparametric estimation  of the density of f
  
  ##--->  Branch 1: Choice of estimation method NP or Ordered choice
  
  # --> if NP estimator
  ## ---> Branch 2: Either  bw is chosed by the user, or CV or Silverman (default) 
  
  # normal-reference rule-of-thumb (Siverman) (always computed to define the object)
  bw.sil <-npudensbw(~uhat,bwmethod="normal-reference")
  dens.np <- npudens(bws=bw.sil, bandwidth.compute =FALSE )
  summary(dens.np)

  if(udensmethod == "CV"){
    bw.cv <-npudensbw(~uhat, ckertype="epanechnikov", bandwidth.compute =TRUE )
    dens.np<- npudens(bws=bw.cv, bandwidth.compute = FALSE)
    summary(dens.np)
  }
  if(udensmethod == "FIXED") {
  # fhat with user-defined bandwidth
    bw.fixed <- bw.sil  # to start with an existing object
    bw.fixed$bw <- ubw
    bw.fixed$bandwidth <- ubw
    dens.np <- npudens(bws=bw.fixed, bandwidth.compute = FALSE)
    summary(dens.np)
}
  # Computing fhta 
  fhat <-predict(dens.np)
  summary(fhat)
  print(paste("---length of fhat",length(fhat), "obs"))
  print(paste("--- Nb of fhat =0",length(which(abs(fhat) < 0.0000001)),"obs"))
  
  
  ## Step 3: Definition of T1
  print("Step3 -- Computation of T1 ") 
  print("-- ")
  vpos <- as.numeric(v >= 0)
  T1 <- ( as.numeric(y) - vpos)/ fhat
  
  print("-- Summary of T1--")
  print(paste("---length of T1",length(T1), "obs"))
  print(summary(T1))
  print(paste("Statut of T1 ",class(T1),"."))
  
  # Correction if Hetero option
  if(hetero=="HETERO"){
    T1 <- T1 * sqrt(abs(xbetahat))   # From step 1bis
  }
  
  # Apply trimming / Winsorization
  T1.trim <- wintrim(T1, trimtype = trimtype, trimlevel = trimlevel)
  # PB sur wintrim !!!!!
  
  trim.U <- 1 - trimlevel/2
  trim.L <- trimlevel/2
  x.bounds <- quantile(T1, probs = c(trim.L, trim.U))
  nb.trim <- length(which(T1 <= x.bounds[1] |  T1 >= x.bounds[2]) )   
  print(paste(" Bounds: ",x.bounds,"."))
  
  print(paste(" Nb de points trimed_out: ", nb.trim,"."))
  
  print("-- Summary of T1.trim --")
  print(paste("---length of T1.trim",length(T1.trim), "obs"))
  print(summary(T1.trim))
  print(" ")
  print(" ")
  print(paste(" Quantiles",quantile(T1, probs = c(trimlevel/2, 1-trimlevel/2), na.rm = TRUE)))
  print(" -")
  print("-- ")
  
  trimmed.out <- length(T1)- length(which(T1.trim != "NA"))
  print(paste("Nb of points trimed out: ",trimmed.out, "."))
  
  
  ## Step 4: 2SLS
  print("Step4 - 2SLS")
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




  # Defining the type of vectors
  T1 <- as.numeric()
 
  # Step 0 : Demean V
  T1 <- T1-mean(T1) 
  
  ## Step 3: Definition of T1
  print("Step3 -- Computation of T1 ") 
  print("-- ")
  
  print("-- Summary of T1--")
  print(paste("---length of T1",length(T1), "obs"))
  print(summary(T1))
  print(paste("Statut of T1::::::::  ",class(T1),"."))
  
  
  # Apply trimming / Winsorization
  T1.trim <- wintrim(T1, trimtype = trimtype, trimlevel = trimlevel)
  # PB sur wintrim !!!!!
  
  trim.U <- 1 - trimlevel/2
  trim.L <- trimlevel/2
  x.bounds <- quantile(T1, probs = c(trim.L, trim.U))
  nb.trim <- length(which(T1 <= x.bounds[1] |  T1 >= x.bounds[2]) )   
  print(paste(" Bounds: ",x.bounds,"."))
  
  print(paste(" Nb de points trimed_out: ", nb.trim,"."))
  
  print("-- Summary of T1.trim --")
  print(paste("---length of T1.trim",length(T1.trim), "obs"))
  print(summary(T1.trim))
  print(" ")
  print(" ")
  print(paste(" Quantiles",quantile(T1, probs = c(trimlevel/2, 1-trimlevel/2), na.rm = TRUE)))
  print(" -")
  print("-- ")
  
  trimmed.out <- length(T1)- length(which(T1.trim != "NA"))
  print(paste("Nb of points trimed out: ",trimmed.out, "."))
  
  
}
  
  



# retour <- specialreg.fit(D,special, NewX, Newiv,  trimtype="WINSOR", trimlevel=0.025,
#                          udensmethod = "FIXED", ubw= 0.23 )
# 
# retour <- specialreg.fit(D,special, NewX, Newiv)






## ---- SAND BOX  ------
# Testing wintrim 

x.u <- rnorm(1000)
summary(x.u)
# a 5%: trim = 0.05
Montrim3 <- wintrim(x.u)   # <- defaut
quantile(x.u, probs=c(0.025, 0.975))
summary(Montrim3)

# a 2.5%: trim = 0.025
Montrim2 <- wintrim(x.u, trimlevel = 0.025, trimtype="TRIM")  
quantile(x.u, probs=c(0.0125, 0.9875))
summary(Montrim2)



  
  
  
  
  
  
  