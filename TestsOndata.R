# Test des différentes fonctions sur nos donnéees


rm(list=ls())


library(np)
library(foreign)
library(AER)
library(Formula)
library(texreg)

#setwd("D:/progs/Celine/Special")   
setwd("c:/Chris/progs/Celine/Special")   

##  Testing wintrim 
# generating the same numbers than in the real dataset
set.seed(1234)
x.u <- rnorm(2771)

trimtest <- 0.005

# Par defaut toutes lmes fonctions sont à 5%: trim = 0.05

trim1 <- wintrim(x.u, trimtype = "WINSOR" , trimlevel = trimtest)   # <- defaut
trim2 <- wintrim2(x.u , trimtype = "WINSOR" , trimlevel = trimtest)   # <- defaut
trim3 <- wintrim.stata(x.u , trimtype = "WINSOR" , trimlevel = trimtest)  # <- defaut
trim4 <- wintrim.test(x.u , trimtype = "WINSOR" , trimlevel = trimtest)  # <- defaut
trim5 <- wintrim.testN(x.u , trimtype = "WINSOR" , trimlevel = trimtest)  # <- defaut


# La borne quantile est presque identiques au max 
quantile(abs(x.u), probs=1-trimtest)

lapply(list(trim1, trim2,trim3,trim4, trim5), summary)


# Les longueurs diffèrent 
length(trim1[!is.na(trim1)])
length(which(!is.na(trim2)))
length(which(!is.na(trim3)))
length(which(!is.na(trim4)))

# lapply(list(trim1, trim2,trim3,trim4), length)
length(which(abs(x.u) < quantile(abs(x.u), probs=1-trimtest)))

# Les fonctions wintrim.stata  et wintrim.test donnent des réultats approchants
# 1 de diference pour  trimlevel = 0.05
# 0 de différence pour trimlevel = 0.025 et 0.005



#  test des fonctions SpecialReg sur nos données 

## Load the dat directly after an estimation in STATA to have same sample 2771 obs. 
data.all<-read.dta("../Water/data/FinalFile.dta")

# # Load data with results (homogenous) !!
# data.all<-read.dta("FinalFileResults1.dta")
# 
# # Load data with results (heterogenous) !!
 data.all<-read.dta("FinalFileResultsHet.dta")


# Loading the chosen Dataset 

nrow(data.all)
data.work <- data.all

attach(data.work)
# Here we prepare data with standard notations 
D <- cbind(i_tap)
endo <- cbind(isatis_health)
iv <- cbind(itap_2008 ,iconcernwatpol_2008)
exo <- cbind( a2_age,i_under18, log_income, i_town, i_car, b08_locenv_water, i_can, i_fra)
hetv <- cbind(i_under18, log_income, i_town ,b08_locenv_water)

#####  Variables def pour test pas a pas TO REMOVE !!!
# y <- i_tap
# v <- Special
# trimtype <-"TRIM"
# trimlevel <- 0.025
# hetero = "HETERO"
# udensmethod = "SILVERMAN"  

# nouvelle formulation avec exo, endo et iv explicites !!! 
retour <- specialreg.fitN(D,Special, endo,  exo, iv,  
                          trimtype="TRIM", trimlevel=0.025 )  

### Test avec hetero !!!

retour2 <- specialreg.fitN(D,Special, endo,  exo, iv,  
                          trimtype="TRIM", trimlevel=0.025,
                          hetero = "HETERO", hetv =hetv )  

# Stata bw for homo = 0.07154983, Hetero = .23179664 

retour3 <- specialreg.fitN(D,Special, endo,  exo, iv,  
                           trimtype="TRIM", trimlevel=0.025,
                           hetero = "HETERO", hetv =hetv, 
                           udensmethod = "SILVERMAN")  


retour3 <- specialreg.fitN(D,Special, endo,  exo, iv,  
                           trimtype="WINSOR", trimlevel=0.025,
                           hetero = "HETERO", hetv =hetv, 
                           udensmethod = "FIXED", 
                           ubw = .23179664   )   


# Old procedure splittted in two : 
retour <- specialreg.fitMarg(D,Special, endo,  exo, iv,  
                           trimtype="TRIM", trimlevel=0.025,
                           hetero = "HETERO", hetv =hetv, 
                           udensmethod = "SILVERMAN" )  

# Test on parameters used in paper (table2) 

retour.fit <- specialreg.fit(D,Special, endo,  exo, iv,  
                             trimtype="TRIM", trimlevel=0.025,
                             hetero = "HETERO", hetv =hetv, 
                             udensmethod = "SILVERMAN" )   


retour.mfx <- specialreg.mfx(D,Special, endo,  exo, iv,  
                             trimtype="TRIM", trimlevel=0.025,
                             hetero = "HETERO", hetv =hetv, 
                             udensmethod = "SILVERMAN" )   

## Attenpt to recover results from table 4 


m1 <-specialreg.fit(D,Special, endo,  exo, iv,  
                    trimtype="WINSOR", trimlevel=0.005,
                    hetero = "HETERO", hetv =hetv, 
                    udensmethod = "SILVERMAN" )   

m2 <-specialreg.fit(D,Special, endo,  exo, iv,  
                    trimtype="TRIM", trimlevel=0.005,
                    hetero = "HETERO", hetv =hetv, 
                    udensmethod = "SILVERMAN" )   


m3 <-specialreg.fit(D,Special, endo,  exo, iv,  
                    trimtype="WINSOR", trimlevel=0.025,
                    hetero = "HETERO", hetv =hetv, 
                    udensmethod = "SILVERMAN" )   

m4 <-specialreg.fit(D,Special, endo,  exo, iv,  
                    trimtype="TRIM", trimlevel=0.025,
                    hetero = "HETERO", hetv =hetv, 
                    udensmethod = "SILVERMAN" )   


m5<-specialreg.fit(D,Special, endo,  exo, iv,  
                    trimtype="WINSOR", trimlevel=0.05,
                    hetero = "HETERO", hetv =hetv, 
                    udensmethod = "SILVERMAN" )   

m6 <-specialreg.fit(D,Special, endo,  exo, iv,  
                    trimtype="TRIM", trimlevel=0.05,
                    hetero = "HETERO", hetv =hetv, 
                    udensmethod = "SILVERMAN" )   


screenreg( l= list(m1, m2, m3, m4, m5, m6), 
           stars = c(0.01, 0.05, 0.10), 
           digits = 3)



# Test avec les résidus estimés de STATA : Cas HETERO

m1.h <-specialreg.fit.stata(D,Special, endo,  exo, iv,  
                    trimtype="WINSOR", trimlevel=0.005,
                    hetero = "HETERO", hetv =hetv, 
                    udensmethod = "SILVERMAN" )   


m2.h <-specialreg.fit.stata(D,Special, endo,  exo, iv,  
                    trimtype="TRIM", trimlevel=0.005,
                    hetero = "HETERO", hetv =hetv, 
                    udensmethod = "SILVERMAN" )   


m3.h <-specialreg.fit.stata(D,Special, endo,  exo, iv,  
                    trimtype="WINSOR", trimlevel=0.025,
                    hetero = "HETERO", hetv =hetv, 
                    udensmethod = "SILVERMAN" )   

m4.h <-specialreg.fit.stata(D,Special, endo,  exo, iv,  
                    trimtype="TRIM", trimlevel=0.025,
                    hetero = "HETERO", hetv =hetv, 
                    udensmethod = "SILVERMAN" )   


m5.h<-specialreg.fit.stata(D,Special, endo,  exo, iv,  
                   trimtype="WINSOR", trimlevel=0.05,
                   hetero = "HETERO", hetv =hetv, 
                   udensmethod = "SILVERMAN" )   

m6.h <-specialreg.fit.stata(D,Special, endo,  exo, iv,  
                    trimtype="TRIM", trimlevel=0.05,
                    hetero = "HETERO", hetv =hetv, 
                    udensmethod = "SILVERMAN" )   

screenreg( l= list(m1.h, m2.h, m3.h, m4.h, m5.h, m6.h), 
           stars = c(0.01, 0.05, 0.10), 
           digits = 3)



retour.mfx.stata <- specialreg.mfx.stata(D,Special, endo,  exo, iv,  
                             trimtype="TRIM", trimlevel=0.025,
                             hetero = "HETERO", hetv =hetv, 
                             udensmethod = "SILVERMAN" )   



################################## OLD STUFF (Usefull though !) ###############)

# On ajoute les X exo aux instruments

NewX <- cbind(isatis_health, i_under18 , log_income ,i_town, i_car, b08_locenv_water, a2_age, i_can, i_fra)
Newiv <- cbind(iv, exo)


# Test de la fonction complète specialreg.fit 

retour <- specialreg.fit(D,Special, NewX, Newiv,  trimtype="TRIM", trimlevel=0.025 )  
# <- bizarre dès la prmeière regression

lm(formula = Special ~ NewX + Newiv)
