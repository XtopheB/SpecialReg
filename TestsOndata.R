# Test des différentes fonctions sur nos donnéees


##  Testing wintrim 
# generating the same numbers than in the real dataset
set.seed(1234)
x.u <- rnorm(2771)

trimtest <- 0.005

# Par defaut toutes lmes fonctions sont à 5%: trim = 0.05

trim1 <- wintrim(x.u, trimtype = "TRIM" , trimlevel = trimtest)   # <- defaut
trim2 <- wintrim2(x.u , trimtype = "TRIM" , trimlevel = trimtest)   # <- defaut
trim3 <- wintrim.stata(x.u , trimtype = "TRIM" , trimlevel = trimtest)  # <- defaut
trim4 <- wintrim.test(x.u , trimtype = "TRIM" , trimlevel = trimtest)  # <- defaut

# La borne quantile est presque identiques au max 
quantile(abs(x.u), probs=1-trimtest)

lapply(list(trim1, trim2,trim3,trim4), summary)


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
nrow(data.all)
data.work <- data.all

attach(data.work)
# Here we prepare data with standard notations 
D <- i_tap
endo <- isatis_health
iv <- cbind(itap_2008 ,iconcernwatpol_2008)
exo <- cbind(i_under18, log_income, i_town, i_car, b08_locenv_water, a2_age, i_can, i_fra)

# On ajoute les X exo aux instruments

NewX <- cbind(isatis_health, i_under18 , log_income ,i_town, i_car, b08_locenv_water, a2_age, i_can, i_fra)
Newiv <- cbind(iv, exo)


# Test de la fonction complète specialreg.fit 

retour <- specialreg.fit(D,Special, NewX, Newiv,  trimtype="TRIM", trimlevel=0.025 )  
    # <- bizarre dès la prmeière regression

lm(formula = Special ~ NewX + Newiv)


# nouvelle formulation avec exo, endo et iv explicites !!! 
retour <- specialreg.fitN(D,Special, endo,  exo, iv,  trimtype="TRIM", trimlevel=0.025 )  


