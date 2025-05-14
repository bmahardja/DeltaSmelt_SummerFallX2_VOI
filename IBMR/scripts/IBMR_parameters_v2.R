############################################################
### Delta smelt IBMR parameters ############################
### William Smith (USFWS; BDFWO); 27 June 2022 #############
############################################################

# v2: added pcr parameter to control predation to consumption effect ratio for turbidity
# - used sum of squares to fit simulated Feb lengths to observed Feb lengths, changed max turbidity effect (min1) to 0.63

### Model parameters ###
# 1. Model dimensions
# 2. Rose and Kimmerer parameters
# 3. New parameters for this version

### 1. Model dimensions ###
n.prey <- 12 # Number of zooplankton prey groups
n.strata <- 12 # Number of spatial strata
n.years <- length(first:length(yearz))

### 2. Rose and Kimmerer parameter values ###
# Rose parameters from Rose and Kimmerer (2013a) Table 1
Sp <- 0.15									# fraction of weight lost due to spawning

# Maximum consumption (C.max) parameters - larvae, postlarvae, juveniles, adults
a.c <- c(0.18,0.18,0.18,0.1,0.1)			# weight multiplier - 0.18 larvae & postlarvae, 0.1 juveniles & adults
b.c <- c(-0.275,-0.275,-0.275,-0.54,-0.54)	# weight exponent - -0.275 larvae & postlarvae, -0.54 juveniles & adults
CQ <- c(7,7,10,10,10)						# Temperature at CK1 of maximum (deg C)
T.0 <- c(17,17,20,20,20)					# Temperature at 0.98 of maximum (deg C)
T.M <- c(20,20,23,23,23)					# Temperature at 0.98 of maximum (deg C)
T.L <- c(28,28,27,27,27)					# Temperature at CK4 of maximum (deg C)
CK.1 <- c(0.4,0.4,0.4,0.4,0.4)				# effect at temperature CQ
CK.4 <- c(0.01,0.01,0.01,0.001,0.001)			# effect at temperature T.L

# Metabolism (R) parameters
a.r <- c(0.0027,0.0027,0.0027,0.0027,0.0027)	# weight multiplier
b.r <- c(-0.216,-0.216,-0.216,-0.216,-0.216)	# weight exponent
R.Q <- c(0.036,0.036,0.036,0.036,0.036)			# exponent for temperature effect
S.d <- c(0.175,0.175,0.175,0.175,0.175)			# Fraction of assimilated food lost to SDA (specific dynamic action)

# Egestion (F) and excretion (U) parameters
F.a <- c(0.16,0.16,0.16,0.16,0.16)				# Fraction of consumed food lost to egestion
U.a <- c(0.1,0.1,0.1,0.1,0.1)				# Fraction of assimilated food lost to excretion

stage.limit <- c(15,25)     # stage maximum lengths
e.s <- 4814									# J/g: convert g(prey)/g(delta smelt) to g(smelt)/g(smelt) - fixed
e.d <- c(1813,2590,2590,2590,2590,2590,2590,2590,2590,2590,2590,2590)		# energy density of prey items

V <- matrix(NA,n.prey,5)					# V.ij is vulnerability of prey type j to fish i. Set to 1 for all life stages eating all zooplankton types; except DS larvae values were 0 except for Limnoithona
V[1,] <- c(1,1,1,1,0) #limno
V[2,] <- c(1,1,1,1,1) #othcaljuv
V[3,] <- c(1,1,1,1,1) #pdiapjuv
V[4,] <- c(0,0,1,1,1) #othcalad
V[5,] <- c(0,0,1,1,1) #acartela
V[6,] <- c(0,0,1,1,1) #othclad
V[7,] <- c(1,1,1,0,0) #allcopnaup
V[8,] <- c(0,0,0,1,1) #daphnia
V[9,] <- c(0,1,1,1,1) #othcyc
V[10,] <- c(0,0,1,1,1) #other
V[11,] <- c(0,0,1,1,1) #eurytem
V[12,] <- c(0,0,1,1,1) #pdiapfor
K <- matrix(NA,n.prey,5)					# K.ik is half-saturation constant for fish i feeding on each prey type k - calculated outside model to obtain realistic diet and consumption rates
K[1,] <- c(NA,2.5,120,1.5,100) #limno				# From 2nd Rose model_apr2020
K[2,] <- c(NA,0.375,0.24,7.5,3) #othcaljuv
K[3,] <- c(NA,0.375,0.24,1.5,2) #pdiapjuv
K[4,] <- c(NA,250,6,0.75,0.6) #othcalad
K[5,] <- c(NA,250,36,0.75,0.25) #acartela
K[6,] <- c(NA,250,120,4.5,1) #othclad
K[7,] <- c(7.5,7.5,120,75,100) #allcopnaup
K[8,] <- c(NA,250,200,4.5,0.15) #daphnia
K[9,] <- c(NA,1.5,1.2,1.5,2) #othcyc
K[10,] <- c(NA,250,12,7.5,3) #other
K[11,] <- c(NA,250,6,0.375,0.25) #eurytem
K[12,] <- c(NA,250,2.4,0.375,0.25) #pdiapfor

### 3. New parameters ###
# kde parameters for super-adults
clean.cyc<-3 # number of loops to clean kde data
n.grid<-12

# Spawning parameters
egg.surv<-1
sp.temp.lo<-9
sp.temp.hi<-17.9 # no spawning beyond this temp

# length-weight parameters
ln.a<-c(0.000005,0.00000183)
ln.b<-c(3,3.38)

# Entrainment parameters
max.E <- c(0.15,0.15) # max daily log (entrainment risk), when OMR<=-5000, given occupancy of South Delta, make high enough that 1-exp(-31*max.E) is near one
entrain.strata<-c(0,0,1,0,0,0,0,0,0,0,0,0)# entrainment only occurs in SDelta strata
entrain.L <- 45
entrain.month<-c(rep(1,6),rep(0,5),1) # entrainment vulnerability by month
# Parameters to scale entrainment as f(OMR) (given occupancy of South Delta)
E.ad.x<-c(-5000,0) # model 'ramp' for -5000<OMR<-2500
E.ad.y<-c(max.E[1],0)
mod1<-glm(E.ad.y~E.ad.x)
E.pl.x<-c(-2500,0) # model 'ramp'
E.pl.y<-c(max.E[2],0)
mod2<-glm(E.pl.y~E.pl.x)
a.E<-vector()
b.E<-vector()
a.E[1] <- mod1$coefficients[1]
b.E[1] <- mod1$coefficients[2]
a.E[2] <- mod2$coefficients[1]
b.E[2] <- mod2$coefficients[2]

# starvation parameter; death if Wt_t<stv*Wt_t-1
stv <- 0.85

# Model turbidity effect on Cmax and natural mortality (M)
# based on Hasenbein et al. (2016), Fig. 2
max.Sec <- 248.51*5^-0.674 # min NTU measured=5, convert to Secchi
mid1.Sec <- 248.51*25^-0.674 # lower NTU at max feeding rate=25, convert to Secchi
min1 <- c(0.63,0.63) #c(8/25,0.85) # lowest Cmax scalar at low turbidity
turb.L.fx.mod <- glm(c(min1[1],min1[2])~c(20,45))
a.turb <- (max.Sec+mid1.Sec)/2 # glm logit regression parameters
b.turb <- c(-0.1,0.1)
pcr <- 1 # predation effect to consumption effect ratio

#temp_x1 <- c(max.Sec,mid1.Sec)
#temp_y1<-c(min1,1)
#mod1<-glm(temp_y1~temp_x1) # glm 'ramp' up from 5 to 35NTU
#a.turb<-vector()
#b.turb<-vector()
#a.turb[1] <- mod1$coefficients[1]
#b.turb[1] <- mod1$coefficients[2]

# Model turbidity effect on mortality
#min2 <- 0.66
#temp_y2<-c(1,min2) #
#mod2 <- glm(temp_y2~temp_x1) # glm 'ramp' down from 5 to 35NTU
#a.turb[2] <- mod2$coefficients[1]
#b.turb[2] <- mod2$coefficients[2]
