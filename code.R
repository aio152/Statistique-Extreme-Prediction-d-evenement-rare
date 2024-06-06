################################## Packages ####################################
library(evd)
library(glue)
library(tseries)
library(urca)
########################## Importation des données #############################
setwd("/home/alassaneissakaomar/Documents/Statistique Spatial")

data = read.table('BERLIN.txt')
data
data=data[data[,1]!='1945',]
data_ts = ts(data[,4], start = 1876, end=2023, frequency = 365)
plot(data_ts, xlab= "Années", ylab = "Temperatures", main  = "Representation de nos données par serie temporelle")
plot(data)


#Test d'egalité de moyenne:
data44 = data[data[,1]=='1944',]
data46 = data[data[,1]=='1946',]
t.test(data46[,4], data44[,4], paired=FALSE)

## Représentation des précipitations en fonction du temps (année)
plot(data[,1],data[,4],type='l')

## Histogramme des précipitations -> permet de voir l'allure de la distribution
hist(data[,4], freq = F, main = "Histogramme des temperatures")
# borné avec une queue très légère


#Test de dickey Foller

result1 <- ur.df(data_ts)
result2 <- adf.test(data_ts)
print(summary(result1))
print(result2)
############################## Approche GEV ####################################
####################### Modele avec temperature maximale #######################

# blocs par année à considérer => conserver les max par année
maxima = NULL
Annees = unique(data[,1])
Annees
for(k in 1:length(Annees)){
  maxima[k] = max(data[data[,1]==Annees[k],4])
}

maxima
plot(Annees,maxima, main = "Tracer des maximum en fonction des années")

## Estimation des paramètres du modèle
modele = fgev(maxima)
par(mfrow = c(2,2))
plot(modele)
dev.off()
## Vraisemblance profilée pour chaque paramètre
# LOC
modele$estimate[1]
modele$estimate[2]
modele$estimate[3]
plot(Vrais)

# intervalle de confiance à 95% = intersection avec droite en pointilles
# Ici environ entre 90 et 100

par(mfrow = c(2,2))
# Loc
Vrais = profile(modele,which = names(modele$estimate[1]),main="")
plot(Vrais)
# SCALE
Vrais = profile(modele,which = names(modele$estimate[2]),main="")
plot(Vrais)

# SHAPE
Vrais = profile(modele,which = names(modele$estimate[3]),main="")
plot(Vrais)

dev.off()

# plutot Weibull mais 0 compris

################################# Niveau de retour ####################################    
Mod10 = fgev(maxima,prob=1/10)
Mod100 = fgev(maxima,prob=1/100)
# Intervalle de confiance
# Mod10 [estm +- 1.96sqrt(var(estm))] = [35.23112 37.12748]; estm = 36.1793
# Mod100 [estm +- 1.96sqrt(var(estm))] = [36.95323 39.65337]; estm = 38.3033




####################### Modele avec temperature minimale #######################
plot(-data)

## Représentation des précipitations en fonction du temps (année)
plot(-data[,1],-data[,4],type='l')

## Histogramme des précipitations -> permet de voir l'allure de la distribution
hist(-data[,4])
# borné avec une queue très légère

# blocs par année à considérer => conserver les max par année
maxima_min = NULL
Annees = unique(data[,1])
Annees
for(k in 1:length(Annees)){
  maxima_min[k] = max(-data[data[,1]==Annees[k],4])
}

maxima_min
plot(Annees,-maxima_min)

## Estimation des paramètres du modèle
modele_2 = fgev(maxima_min)
par(mfrow = c(2,2))
plot(modele_2)
dev.off()
## Vraisemblance profilée pour chaque paramètre
# Estimation des parametres
-modele_2$estimate[1]
modele_2$estimate[2]
modele_2$estimate[3]

# intervalle de confiance à 95% = intersection avec droite en pointilles
# Ici environ entre 90 et 100

# Loc
Vrais = profile(modele_2,which = names(modele_2$estimate[1]),main="")
plot(Vrais)
# SCALE
Vrais = profile(modele,which = names(modele_2$estimate[2]),main="")
plot(Vrais)

# SHAPE
Vrais = profile(modele,which = names(modele_2$estimate[3]),main="")
plot(Vrais)
# plutot Weibull mais 0 compris

# Estimation des niveau de retour
Mod_2_10 = fgev(maxima_min,prob=1/10)
Mod_2_100 = fgev(maxima_min,prob=1/100)

# Intervalle de confiance
# Mod10 [estm ± 1.96sqrt(var(estm))] = [-13.64652 -11.03708]; estm = -12.3418
# Mod100 [estm ± 1.96sqrt(var(estm))] = [-18.3869 -14.5467]; estm = -16.4668 


######################### GPD #######################

######################### Temperature maximale #######################


#choix de seuil pour visualiser graphiquement les seuils
u_1 = 30
u_2 = 34


# Représentation des précipitations en fonction du temps (année)
plot(data_ts)
abline(u_1,0, col = 'red')
abline(u_2,0, col = 'blue')

##################### choix du seuil ##################

# Définir le vecteur de seuils entre 25 et 40
#Choix de r = 10 pour les clusters
seuils <- seq(25, 35, by = 0.5)

# Initialiser les vecteurs pour stocker les résultats
variance_values <- numeric(length(seuils))
shape_values <- numeric(length(seuils))

# Boucle pour chaque seuil
for (i in 1: length(seuils)) {
  u <- seuils[i]
  # Calculer le modèle GPD avec fpot
  model_GPD <- fpot(data[,4],threshold =  u, r = 10, cmax=TRUE)
  
  # Stocker les valeurs de variance et shape
  variance_values[i] <- model_GPD$estimate[1]
  
  shape_values[i] <- model_GPD$estimate[2]
}

# Tracer la variance en fonction de u
plot(seuils, variance_values, type = "o", col = "blue", xlab = "Seuil (u)", ylab = "Variance", main = "Tracer du seuil en fonction de la variance")

# Tracer le shape en fonction de u
plot(seuils, shape_values, type = "o", col = "red", xlab = "Seuil (u)", ylab = "Shape", main = "Tracer du seuil en fonction de la forme")

r1 = 10
r2 = 2

############## Representation des cluster ################################
par(mfrow = c(2,1))
clusters(data[,4],u_1,r1,plot=TRUE, col = 'white',xlab='Années', main = "u=30 , r =10")
clusters(data[,4],u_2,r2,plot=TRUE, col = 'white',xlab='Années', main = "u=34 , r =2")
dev.off()


######### Estimation des paramètres d'un modèle GPD pour chaque cas #########

model_GPD1 = fpot(data[,4],u_1, r = r1, cmax=TRUE)
par(mfrow = c(2,2))
plot(model_GPD1)

model_GPD2 = fpot(data[,4],u_2, r = r2, cmax=TRUE)
par(mfrow = c(2,2))
plot(model_GPD2)
dev.off()

########### # Estimation de la valeur de l'indice extremal theta et interprétation#################
exi(data[,4],u=u_1, r = r1) # estm 0.3129003
exi(data[,4],u=u_2, r = r2) # estm 0.7

################################### REtour ################################
# Retour pour 10 ans
retour_10 = fpot(data[,4],u_1, npp=length(data[,4]), mper = 10, cmax=TRUE,r = r1)
v1c =fpot(data[,4],u_2, npp=length(data[,4]), mper = 10, cmax=TRUE, r= r2)
plot(v1c)

# Retour pour 100 ans
retour_100 = fpot(data[,4],u_1, npp=length(data[,4]), mper = 100, cmax=TRUE, r = r1)
fpot(data[,4],u_2, npp=length(data[,4]), mper = 100, cmax=TRUE, r = r1)


########################### Intervalle de confiance ########################
#  retour_10    [37.69572, 40.61228] ; estm 39.154
#  retour_100     [37.91106, 41.11814]; estm 39.5146


######################### Temperature minimale #######################

#choix de seuil pour visualiser graphiquement les seuils
v_1 = 10
v_2 = 0


# Représentation des précipitations en fonction du temps (année)
plot(-data_ts)
abline(v_1,0, col = 'red')
plot(-data_ts)
abline(v_2,0, col = 'red')

##################### choix du seuil ##################

# Définir le vecteur de seuils entre 0 et -10
#Choix de r = 10 pour les clusters
seuils <- seq(0, 10, by = 0.5)

# Initialiser les vecteurs pour stocker les résultats
variance_values <- numeric(length(seuils))
shape_values <- numeric(length(seuils))

# Boucle pour chaque seuil
for (i in 1: length(seuils)) {
  u <- seuils[i]
  # Calculer le modèle GPD avec fpot
  model_GPD <- fpot(-data[,4],threshold =  u, r = 10, cmax=TRUE)
  
  # Stocker les valeurs de variance et shape
  variance_values[i] <- model_GPD$std.err[1]
  
  shape_values[i] <- model_GPD$std.err[2]
}

# Tracer la variance en fonction de u
plot(seuils, variance_values, type = "o", col = "blue", xlab = "Seuil (u)", ylab = "Variance", main = "Tracer du seuil en fonction de la variance")

# Tracer le shape en fonction de u
plot(seuils, shape_values, type = "o", col = "red", xlab = "Seuil (u)", ylab = "Shape", main = "Tracer du seuil en fonction de la forme")

r1 = 10

############## Representation des cluster ################################
clusters(data[,4],v_1,r1,plot=TRUE, col = 'white',xlab='Années')
dev.off()


######### Estimation des paramètres d'un modèle GPD pour chaque cas #########

model_GPD1_min = fpot(data[,4], 4, r = r1, cmax=TRUE)
par(mfrow = c(2,2))
plot(model_GPD1)
dev.off()

############ Estimation de la valeur de l'indice extremal theta et interprétation#################
exi(data[,4],u=4, r = r1) # estm 0.006536838

################################### REtour ################################
# Retour pour 10 ans
retour_10 = fpot(-data[,4],4, npp=length(data[,4]), mper = 10, cmax=TRUE,r = r1)


# Retour pour 100 ans
retour_100 = fpot(-data[,4],4, npp=length(data[,4]), mper = 100, cmax=TRUE, r = r1)

########################### Intervalle de confiance ########################
#  retour_10    [-16.58348, -22.10412] ; estm -19.3438
#  retour_100     [-24.23373, -20.31373]; estm -20.6535
