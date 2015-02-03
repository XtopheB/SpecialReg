# 26/01/2015  : Functions tested and up-to-date
# 28/01/2015  : Silverman's rule for bandwidth
# 28/01/2015  : Specialreg with marginal effects


library(np)
library(foreign)
library(AER)
library(Formula)
#library(sem)
library(ks)  # for kernel estimation ...




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
  }
  if(udensmethod == "SILVERMAN") {
    # As in Stata : see http://fmwww.bc.edu/repec/bocode/k/kdens.pdf
    delta.k <- (3/(5*sqrt(5)))^(0.2)
    n.5 <- (length(uhat))^(-0.2)
    sigma.x <-min(sd(uhat, na.rm= FALSE), IQR(uhat, na.rm = FALSE, type = 7)/1.349)
    bw.sil <- 1.159* delta.k * sigma.x * n.5
    bw.dens <- npudensbw(uhat, bws=bw.sil , ckertype="epanechnikov", bandwidth.compute=FALSE) 
  }
  # bandwidth information   
  print(" --- bandwith used ---")
  print(" ---------------")
  print(summary(bw.dens))
  
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



### Version 02/02/2015


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

silverman.bw <- function(x){
  delta.k <- (3/(5*sqrt(5)))^(0.2)
  n.5 <- (length(x))^(-0.2)
  sigma.x <-min(sd(x, na.rm= TRUE), IQR(x, na.rm = TRUE, type = 7)/1.349)
  band <- 1.159* delta.k * sigma.x * n.5
  return(band)
}


# 30/01/2015 : Fonctionne correctement  y compris mfxs...

specialreg.fitMarg <- function(y,v, endo, exo,  iv,
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
  }
  if(udensmethod == "SILVERMAN") {
    # As in Stata : see http://fmwww.bc.edu/repec/bocode/k/kdens.pdf
    
    bw.sil <- silverman.bw(uhat)
    bw.dens <- npudensbw(uhat, bws=bw.sil , ckertype="epanechnikov", bandwidth.compute=FALSE) 
  }
  # bandwidth information   
  print(" --- bandwith used ---")
  print(" ---------------")
  print(summary(bw.dens))
  
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
  y.num <- as.numeric(y==1)
  T1 <- ( y.num - vpos)/ fhat   
  
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
  
  #   names(spe.ivreg$coefficients)<- c("Constant",colnames(endo), colnames(exo))   # marche pas  !!!!
  # print(summary(spe.ivreg))
  print(screenreg(spe.ivreg, digits = 3))
  
  # print(paste("Formula in ivreg : ",spe.ivreg$call, "."))
  # 
  print(paste("Number of points : ",spe.ivreg$nobs,"obs."))
  # print(paste("Sample used in ivreg : ",length(spe.ivreg$residuals), "points."))
  
  
  #return(spe.ivreg)
  
  ## Step 5: Computing marginal effects 
  # we keep close to stata notations 
  beta.hat <- as.matrix(spe.ivreg$coefficients)
  xbeta.hat <- spe.ivreg$fitted.values
  length(xbeta.hat) <- length(y)    # One has to have the same length  (Trimming + estimation)
  # Computing the bandwidth 
  
  aif.bw <- silverman.bw(xbeta.hat)
  
  #Computing Mi the average index function (equ. A2.3 p 25)
  
  aif.hat <- npreg(tydat=y, txdat=xbeta.hat, bws = aif.bw,
                   ckertype="epanechnikov") 
  
  
  print(" Summary of AIF hat")
  print(" ---------------")
  print(summary(aif.hat))
  
  # Estimated AIF in equation  A2.3
  Mi <- fitted(aif.hat)
  length(Mi) <- length(v)    # One has to have the same length  (Trimming + estimation)
  
  #computing equation A2.2
  DMi <- y.num - Mi
  num <- npksum(txdat=xbeta.hat,
                tydat= DMi,
                operator="derivative",
                bws = aif.bw,
                bandwidth.divide=TRUE,
                ckertype="epanechnikov",
                ckerorder=2)
  
  deno <-  npksum(txdat=xbeta.hat,
                  bws = aif.bw,
                  bandwidth.divide=FALSE,
                  ckertype="epanechnikov",
                  ckerorder=2)
  
  mfxi <- num$ksum /deno$ksum
  
  # we have to add 1 to beta for V 
  mfx <- mfxi  %*% t(beta.hat)  #  Then Special is at the end (PB with names !!)
  
  
  print(" Summary of mfx (with only Xbeta (and beta)) ")
  print(" ---------------")
  print(t(t(colMeans(mfx))))
  #return(mfx)
  
  ### NEW version :   Adding one to beta at the end  !!! 
  beta.hat.1 <- c(1, beta.hat)
  names(beta.hat.1) <- c("Special", "cons", colnames(endo), colnames(exo))
  
  mfx.1 <- mfxi  %*% t(beta.hat.1)  #  
  print(" Summary of mfx (with  Xbeta and Beta+1)")
  print(" ---------------")
  print(t(t(colMeans(mfx.1))))
  #return(mfx)
  
  
  ##### Correct ? version with Xbeta + V !!!!!! 
  
  #Computing Mi the average index function (equ. A2.3 p 25)
  
  # In here we use Xbeta +v for computation  
  xbetav.hat <- xbeta.hat +v 
  betav.hat <- c(1, beta.hat)
  length(xbetav.hat) <- length(y)    # One has to have the same length  (Trimming + estimation)
  
  #Compution the bandwidth
  aifv.bw <- silverman.bw(xbetav.hat)
  
  #Computing Mi the average index function (equ. A2.3 p 25)
  
  aifv.hat <- npreg(tydat=y, txdat=xbetav.hat, bws = aifv.bw,
                    ckertype="epanechnikov") 
  
  print(" Summary of AIF hat")
  print(" ---------------")
  print(summary(aifv.hat))
  
  # Estimated AIF in equation  A2.3
  Miv <- fitted(aifv.hat)
  length(Miv) <- length(v)    # One has to have the same length  (Trimming + estimation)
  
  #computing equation A2.2
  DMiv <- y.num - Miv
  numv <- npksum(txdat=xbetav.hat,
                 tydat= DMiv,
                 operator="derivative",
                 bws = aifv.bw,
                 bandwidth.divide=TRUE,
                 ckertype="epanechnikov",
                 ckerorder=2)
  
  denov <-  npksum(txdat=xbetav.hat,
                   bws = aifv.bw,
                   bandwidth.divide=FALSE,
                   ckertype="epanechnikov",
                   ckerorder=2)
  
  mfxiv <- numv$ksum /denov$ksum
  
  # recalling the names..
  names(betav.hat) <- c("Special", "cons", colnames(endo), colnames(exo))
  
  mfxv <- mfxiv  %*% t(betav.hat ) 
  # Estimating the average marginal effects   
  
  print(" Summary of mfx (computed with Xbeta+V (and beta+1)")
  print(" ---------------")
  print(t(t(colMeans(mfxv))))
  print(" ")
  
  # Comparaison des deux approches  
  print(" Comparing  mfx computed with Xbeta+V or Xbeta (and adding 1 to beta")
  print(" ---------------")
  print(cbind(colMeans(mfxv), colMeans(mfx.1)))
  
  return(mfxv)
  
  ## END
}

# On sépare les fonctions : 


specialreg.fit<- function(y,v, endo, exo,  iv,
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
  # 30/01/2015 : V1.0 
{
  # Defining the type of vectors
  y <- as.factor(y)
  y.num <- as.numeric(y==1)
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
  }
  if(udensmethod == "SILVERMAN") {
    # As in Stata : see http://fmwww.bc.edu/repec/bocode/k/kdens.pdf
    
    bw.sil <- silverman.bw(uhat)
    bw.dens <- npudensbw(uhat, bws=bw.sil , ckertype="epanechnikov", bandwidth.compute=FALSE) 
  }
  # bandwidth information   
  print(" --- bandwith used ---")
  print(" ---------------")
  print(summary(bw.dens))
  
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
  T1 <- ( y.num - vpos)/ fhat   
  
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
  
  #   names(spe.ivreg$coefficients)<- c("Constant",colnames(endo), colnames(exo))   # marche pas  !!!!
  # print(summary(spe.ivreg))
  print(screenreg(spe.ivreg, digits = 3))
  
  # print(paste("Formula in ivreg : ",spe.ivreg$call, "."))
  # 
  print(paste("Number of points : ",spe.ivreg$nobs,"obs."))
  # print(paste("Sample used in ivreg : ",length(spe.ivreg$residuals), "points."))
  
  
  return(spe.ivreg)
}




specialreg.mfx<- function(y,v, endo, exo,  iv,
                          trimtype = "TRIM",
                          trimlevel = 0.05, 
                          hetero = "", 
                          hetv = "", 
                          udensmethod = "CV",
                          ubw= "" , 
                          data = data.work,
                          ...
)
{ 
  # 30/01/2015 : built on the previous one
  spe.ivreg <- specialreg.fit(y = y,v=v, endo=endo,
                              exo =exo,  
                              iv=iv,
                              trimtype = trimtype,
                              trimlevel = trimlevel, 
                              hetero = hetero, 
                              hetv = hetv, 
                              udensmethod = udensmethod,
                              ubw= ubw , 
                              data = data,
                              ...)
  
  # I think I don't need this, but error if not !!
  y <- as.factor(y)
  y.num <- as.numeric(y==1)
  v <- as.numeric(v)  
  endo <- as.matrix(endo)
  exo <- as.matrix(exo)
  iv <- as.matrix(iv)
  
  
  ## Step 5: Computing marginal effects 
  # we keep close to stata notations 
  beta.hat <- as.matrix(spe.ivreg$coefficients)
  xbeta.hat <- spe.ivreg$fitted.values
  length(xbeta.hat) <- length(y)    # One has to have the same length  (Trimming + estimation)
  # Computing the bandwidth 
  
  aif.bw <- silverman.bw(xbeta.hat)
  
  #Computing Mi the average index function (equ. A2.3 p 25)
  
  aif.hat <- npreg(tydat=y, txdat=xbeta.hat, bws = aif.bw,
                   ckertype="epanechnikov") 
  
  
  print(" Summary of AIF hat")
  print(" ---------------")
  print(summary(aif.hat))
  
  # Estimated AIF in equation  A2.3
  Mi <- fitted(aif.hat)
  length(Mi) <- length(v)    # One has to have the same length  (Trimming + estimation)
  
  #computing equation A2.2
  DMi <- y.num - Mi
  num <- npksum(txdat=xbeta.hat,
                tydat= DMi,
                operator="derivative",
                bws = aif.bw,
                bandwidth.divide=TRUE,
                ckertype="epanechnikov",
                ckerorder=2)
  
  deno <-  npksum(txdat=xbeta.hat,
                  bws = aif.bw,
                  bandwidth.divide=FALSE,
                  ckertype="epanechnikov",
                  ckerorder=2)
  
  mfxi <- num$ksum /deno$ksum
  #   
  #   mfx <- mfxi  %*% t(beta.hat)  
  #   
  #   print(" Summary of mfx (with only Xbeta (and beta)) ")
  #   print(" ---------------")
  #   print(t(t(colMeans(mfx))))
  #   #return(mfx)
  
  ### NEW version :   Adding one to beta at the end  !!! 
  beta.hat.1 <- c(1, beta.hat)
  names(beta.hat.1) <- c("Special", "cons", colnames(endo), colnames(exo))
  
  mfx.1 <- mfxi  %*% t(beta.hat.1)  #  
  print(" Summary of mfx (with  Xbeta and Beta+1)")
  print(" ---------------")
  print(t(t(colMeans(mfx.1))))
  #return(mfx)
  
  
  ##### Correct ? version with Xbeta + V !!!!!! 
  
  #Computing Mi the average index function (equ. A2.3 p 25)
  
  # In here we use Xbeta +v for computation  
  xbetav.hat <- xbeta.hat +v 
  betav.hat <- c(1, beta.hat)
  length(xbetav.hat) <- length(y)    # One has to have the same length  (Trimming + estimation)
  
  #Compution the bandwidth
  aifv.bw <- silverman.bw(xbetav.hat)
  
  #Computing Mi the average index function (equ. A2.3 p 25)
  
  aifv.hat <- npreg(tydat=y, txdat=xbetav.hat, bws = aifv.bw,
                    ckertype="epanechnikov") 
  
  print(" Summary of AIF hat")
  print(" ---------------")
  print(summary(aifv.hat))
  
  # Estimated AIF in equation  A2.3
  Miv <- fitted(aifv.hat)
  length(Miv) <- length(v)    # One has to have the same length  (Trimming + estimation)
  
  #computing equation A2.2
  DMiv <- y.num - Miv
  numv <- npksum(txdat=xbetav.hat,
                 tydat= DMiv,
                 operator="derivative",
                 bws = aifv.bw,
                 bandwidth.divide=TRUE,
                 ckertype="epanechnikov",
                 ckerorder=2)
  
  denov <-  npksum(txdat=xbetav.hat,
                   bws = aifv.bw,
                   bandwidth.divide=FALSE,
                   ckertype="epanechnikov",
                   ckerorder=2)
  
  mfxiv <- numv$ksum /denov$ksum
  
  # recalling the names..
  names(betav.hat) <- c("Special", "cons", colnames(endo), colnames(exo))
  
  mfxv <- mfxiv  %*% t(betav.hat ) 
  # Estimating the average marginal effects   
  
  print(" Summary of mfx (computed with Xbeta+V (and beta+1)")
  print(" ---------------")
  print(t(t(colMeans(mfxv))))
  print(" ")
  
  # Comparaison des deux approches  
  print(" Comparing  mfx computed with Xbeta+V or Xbeta (and adding 1 to beta")
  print(" ---------------")
  print(cbind(colMeans(mfxv), colMeans(mfx.1)))
  
  return(mfxv)
  
  ## END
}


