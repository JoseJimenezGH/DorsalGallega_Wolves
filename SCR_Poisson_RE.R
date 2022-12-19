#==============================================================================#
#                                                                              #
#                             CAPTURA-RECAPTURA SECR                           #
#                           LOBOS EN LA SIERRA DE SUIDO                        #
#                                - PONTEVEDRA -                                #    
#                                 José Jiménez                                 #
#                                 Daniel Cara                                  #
#                            03/02/2013 11:12:20 AM                            #
#                                                                              #
#==============================================================================#

# Se ha llevado a cabo un estudio experimental de estima de tamaño de la 
# población de lobo en la Sierra de Suido (Pontevedra) a partir de los datos
# de fototrampeo obtenidos a lo largo del año 2015. Para ello se ha tratado
# la información (historiales de captura, ubicación de las cámaras de
# fototrampeo y periodos de operatividad de las trampas) con una aproximación
# bayesiana de captura-recaptura espacialmente explícita (SCR).

library(coda)
library(scrbook)
library(secr)
library(lattice)
library(mcmcOutput)
library(raster)
library(rgdal)
library(fields)
library(ggplot2)

# PREPARACIÓN DE DATOS
#========================
# Preparamos la máscara de operatividad
setwd('C:/Users/Usuario/OneDrive/30 Proyecto Lobo Suido/01 data/')

# El número total de ocasiones de muestreo utilizados fue de 140 dias
MASK<-data.matrix(read.table("MASK.txt", header=FALSE))

source("Funciones_SCR.R")


library(secr)
wolf.ch <- read.capthist("captsecr.txt", "trapssecr.txt", detector='count', noccasions=119)
summary(wolf.ch)

traplocs<-as.matrix(traps(wolf.ch))
X<-data.matrix(traplocs)/1000
medX<- mean(X[,1])
medY<- mean(X[,2])
X[,1]<-(X[,1]-medX)
X[,2]<-(X[,2]-medY)

plot(traplocs, asp=TRUE)
plot(wolf.ch, border = 0, tracks = TRUE, varycol = TRUE, add=TRUE)
wolf<-as.array(wolf.ch)
wolf<- aperm(wolf,c(1,3,2))
str(wolf)

datYknown<- wolf

# Operatividad total (cámaras). Una novedad en este estudio fue la 
# movilidad de las cámaras. Se fueron desplazando para maximizar el número de
# fotocapturas
cat("Camera trapping days=", sum(MASK), "\n") # Cámaras-día
# Media de dias de captura por cámara y trampa
cat("Mean trapping days per camera=", mean(apply((MASK),1,sum)), "\n")

K<-dim(MASK)[2] ;K     # Ocasiones de muestreo
J<-dim(MASK)[1] ;J     # Número de ubicaciones de las cámaras  

# Ploteamos la operatividad de las cámara (amarillo: operativas; rojo: no
# operativas)
x<-as.matrix(MASK)  ## Tenemos que convertir esto en matriz
image(1:K,1:J,t(x), yaxt = "n", xlab="Occasion", ylab="", cex.lab=1.25)
mtext(side = 2, "Camera", line = 2.5, cex=1.25)
axis(2, rev(seq(1, J, by=2)))
box()

# Resumen de la operatividad por ocasiones de muestreo
KT<-apply(MASK,1,sum)
colnames(MASK)<-c(1:K)
wtraps<-cbind(1:J,X,MASK)


M<-100 # Aumentado de datos que usamos (aproximadamente 3 x máximo estimado)

## Aumentado de datos
nind<- dim(datYknown)[1]; nind     # Es el número de individuos identificados

SCRsmy(datYknown)
Yaug <- array(0, dim = c(M, J, K)) # Aumentamos con M-nind "virtuales"

Yaug[1:nind, , ] <- datYknown
y<-apply(Yaug, c(1,2), sum)

# Capturas por individuo. Se observa como el número de capturas reflejan muy
# bien las dos manadas estudiadas, con un gran número de individuos en la 
# primera (probablemente consituida por agregación de dos manadas) y/o la
# presencia recurrente de individuos no territoriales
tab<-table(1:nind)
pts<-barplot(apply(apply(datYknown, c(1,3), sum),1,sum))
axis(side=1, at=pts, labels=F, tick=T)
axis(1, at = pts, labels = names(tab))

# Capturas por trampa.
tab<-table(1:J)
pts<-barplot(apply(apply(datYknown, c(2:3), sum),1,sum))
axis(side=1, at=pts, labels=F, tick=T)
axis(1, at = pts, labels = names(tab))

# Capturas por fecha. Hay un incremento de las capturas por trampa a medida
# que se van desplazando, con el conocimiento adquirido
tab<-table(1:K)
pts<-barplot(apply(apply(datYknown, c(2,3), sum),2,sum))
axis(side=1, at=pts, labels=F, tick=T)
axis(1, at = pts, labels = names(tab))


# COMPROBACIÓN DE ERRORES
#==========================
yr<-apply(datYknown,c(2,3),sum)
sum(yr)
MASK2<-MASK
MASK2[MASK2>0]<-1
sum(yr*MASK2)

# Para localizar los valores negativos (errores en MASK) usamos:
which((yr*MASK2)-yr<0, arr.in=TRUE) # No hay errores en MASK

# Capturas por dispositivo/dia. Rate of capture success
# Camera-trap
cat("Captures/100 trap-days=", 100*sum(yr)/(sum(MASK)),"\n")

# Generación del espacio de estados (S)
trapShape<-readOGR(dsn="C:/Users/Usuario/OneDrive/30 Proyecto Lobo Suido/01 data/GIS",layer="traps")
buffTrap<-buffer(trapShape, width = 9783.9)
plot(buffTrap)
points(trapShape)
A<-area(buffTrap)/1e8  # Aqui sale en ud/100 km2

(buff<- 9783.9/1000)
xl<-min(X[,1])-buff
xu<-max(X[,1])+buff
yl<-min(X[,2])-buff
yu<-max(X[,2])+buff
xlim=c(xl, xu)
ylim=c(yl, yu)

# Ploteado de capturas. Circulos proporcionales al numero de capturas
plot(traplocs, pch="+", col="blue", 
  xlim=c(min(traplocs[,1])-2000,max(traplocs[,1])+2000), 
  ylim=c(min(traplocs[,2])-2000,max(traplocs[,2])+2000), type="n", asp=TRUE)
tot<-apply(datYknown, 2,sum)
symbols(traplocs, circles=tot*50, inches=F,bg="#228B2219", fg=NULL, add=T)
points(traplocs, pch="+", col="blue")
# Spiderplot
spiderplotJJ5(datYknown, traplocs, buffer=2000, lwd=1)



library(nimble)
## define the model
code <- nimbleCode({
  
  alpha1 ~ dnorm(0,.01)
  sigma<- sqrt(1/(2*alpha1))
  psi ~ dunif(0,1)  # Probabilidad de que un ind. esté en la población
  # Random effects hyperparameter
  sigma.p ~ dunif(0, 10)
  mu0 ~ dnorm(0, 0.1)
  # Baseline detection rate
  for(j in 1:J){
    lp[j] ~ dnorm(mu0, sd=sigma.p)
    log(p0[j]) <- lp[j]
  }
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    d2[i,1:J]<- (s[i,1] - X[1:J,1])^2 + (s[i,2] - X[1:J,2])^2
    lam[i,1:J]<- p0[1:J]*exp(- alpha1*d2[i,1:J])
    mu[i,1:J] <- lam[i,1:J]*z[i]
    mu2[i,1:J] <- lam[i,1:J]*z[i]*KT[1:J]
    
    for(j in 1:J){
      outj[i,j] <- sqrt(d2[i,j]) > buffer   
      Y[i,j] ~ dpois(mu2[i,j])
      # components for fit diagnostics
      Ysim[i,j] ~ dpois(mu[i,j]*KT[j])       # simulated
      Yexp[i,j] <- lam[i,j] * z[i] * KT[j]   # expected
      # components for T1
      err1obs[i,j] <- (sqrt(Y[i,j]) - sqrt(Yexp[i,j]))^2
      err1sim[i,j] <- (sqrt(Ysim[i,j]) - sqrt(Yexp[i,j]))^2
    }
    # components for T2
    err2obs[i] <- (sqrt(sum(Y[i,1:J])) - sqrt(sum(Yexp[i,1:J])))^2
    err2sim[i] <- (sqrt(sum(Ysim[i,1:J])) - sqrt(sum(Yexp[i,1:J])))^2

    out[i] <- equals(sum(outj[i,1:J]),J)
    zeros[i] ~ dbern(out[i])
  }
    # components for T3
  for(j in 1:J){
    err3obs[j] <- (sqrt(sum(Y[1:M,j])) - sqrt(sum(Yexp[1:M,j])))^2
    err3sim[j] <- (sqrt(sum(Ysim[1:M,j])) - sqrt(sum(Yexp[1:M,j])))^2
  }
  
  # Fit diagnostics totals
  T1obs <- sum(err1obs[1:M,1:J])
  T1sim <- sum(err1sim[1:M,1:J])
  T2obs <- sum(err2obs[1:M])
  T2sim <- sum(err2sim[1:M])
  T3obs <- sum(err3obs[1:J])
  T3sim <- sum(err3sim[1:J])  

  N<-sum(z[1:M])
  D<-N/A
})

str(constants<-list(M=M, 
                    J=J, 
                    A=A))

str(data   <-  list(Y=y, 
                    X=X, 
                    KT=KT,
                    zeros=rep(0,M), 
                    buffer=buff,  
                    xlim=xlim, 
                    ylim=ylim))


str(inits  <-  list(alpha1=runif(1,0.01,0.06),
                    psi=runif(1,0,1),
                    sigma.p=runif(1,0,3),
                    mu0=runif(1,-5,-1),
                    Ysim=y,
                    s=X[sample(1:J, size=M, replace=TRUE),],
                    z=c(rep(1, nind), rbinom((M-nind),1,0.2))))

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
Rmodel$initializeInfo()
simNodes <- 'lp'
simNodeScalar <- Rmodel$expandNodeNames(simNodes)
allNodes <- Rmodel$getNodeNames()
nodesSorted <- allNodes[allNodes %in% simNodeScalar]
set.seed(1) # This makes the simulations here reproducible
for(n in nodesSorted) {
  Rmodel$simulate(n)
  depNodes <- Rmodel$getDependencies(n)
  Rmodel$calculate(depNodes)
}
Rmodel$lp
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)

params<-c('N','D','sigma','mu0','sigma.p','s','z','psi','T1obs','T1sim','T2obs','T2sim','T3obs','T3sim')
mcmc<-configureMCMC(Rmodel,  monitors=params, enableWAIC = TRUE, thin = 5) 
# Rebuild and compile with new sampler
mcmc$removeSamplers("s")
ACnodes <- paste0("s[", 1:constants$M, ", 1:2]")
for(node in ACnodes) {
  mcmc$addSampler(target = node,
                  type = "AF_slice",
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
}
mcmc$removeSamplers('z')
for(node in Rmodel$expandNodeNames('z')) mcmc$addSampler(target = node, type = 'slice')

MCMC <- buildMCMC(mcmc) 

CompMCMC <- compileNimble(MCMC, project = Rmodel)


## MCMC settings
nb<-5000*5
ni<-25000*5 +nb*5
nc<-3


outNim <- runMCMC(CompMCMC, niter = ni , nburnin = nb , inits=inits, 
                  nchains = nc, WAIC=TRUE, setSeed = FALSE, 
                  progressBar = TRUE, samplesAsCodaMCMC = TRUE)

outNim$WAIC
## nimbleList object of type waicList
## Field "WAIC":
## [1] 931.1034
## Field "lppd":
## [1] -363.3052
## Field "pWAIC":
## [1] 102.2465

#save(outNim, file="Lobo_Poisson_RE_GoF.RData")

summary(window(outNim$samples[,c('N','D','sigma','mu0','sigma.p','psi','T1obs','T1sim','T2obs','T2sim','T3obs','T3sim')],1))



samples<-as.matrix(outNim$samples)

# Assessing fit of the observation model
#========================================
par(mfrow=c(2,2),
    mar=c(6,5,4,4), las=1)  
options(scipen=999)
# Fit statistic 1: individual x trap frequencies.
(pvalue<-round(mean(samples[,'T1sim'] > samples[,'T1obs']), 2))
plot(samples[,'T1sim'], samples[,'T1obs'], xlim=c(70,170), ylim=c(70,170),
  xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", 
  bty ="n", las=1, col="grey50", main="Individual x trap\n frequencies", type='n',
  sub=paste("p-value=", pvalue,""))
idx1 = samples[,'T1sim']-samples[,'T1obs']>0
idx2 = samples[,'T1sim']-samples[,'T1obs']<0   
points(samples[,'T1sim'][idx2],samples[,'T1obs'][idx2], col='#00000033', pch=1, cex=.5)
points(samples[,'T1sim'][idx1],samples[,'T1obs'][idx1], col='#FF000033', pch=1, cex=.5)
segments(80,80,350,350, lwd=2)

# Fit statistic 2: Individual encounter frequencies.
(pvalue<-round(mean(samples[,'T2sim'] > samples[,'T2obs']), 2))
plot(samples[,'T2sim'], samples[,'T2obs'], xlim=c(5,27), ylim=c(5,27),
  xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", 
  bty ="n", las=1, col="grey50", main="Individual encounter\n frequencies", type='n',
  sub=paste("p-value=", pvalue,""))
idx1 = samples[,'T2sim']-samples[,'T2obs']>0
idx2 = samples[,'T2sim']-samples[,'T2obs']<0   
points(samples[,'T2sim'][idx2],samples[,'T2obs'][idx2], col='#00000033', pch=1, cex=.5) 
points(samples[,'T2sim'][idx1],samples[,'T2obs'][idx1], col='#FF000033', pch=1, cex=.5)
segments(0,0,100,100, lwd=2)

# Fit statistic 3: Trap frequencies.
(pvalue<-round(mean(samples[,'T3sim'] > samples[,'T3obs']), 2))
plot(samples[,'T3sim'], samples[,'T3obs'], xlim=c(5,27), ylim=c(5,27), 
  xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", 
  bty ="n", las=1, col="grey50", main="Trap encounter\n frequencies", type='n',
  sub=paste("p-value=", pvalue,""))
idx1 = samples[,'T3sim']-samples[,'T3obs']>0
idx2 = samples[,'T3sim']-samples[,'T3obs']<0   
points(samples[,'T3sim'][idx2],samples[,'T3obs'][idx2], col='#00000033', pch=1, cex=.5)
points(samples[,'T3sim'][idx1],samples[,'T3obs'][idx1], col='#FF000033', pch=1, cex=.5)
segments(0,0,120,120, lwd=2)

