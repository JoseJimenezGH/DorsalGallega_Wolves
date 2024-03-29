#==============================================================================#
#                                                                              #
#                 VIDEO CAMERA TRAP AND SPATIAL CAPTURE-RECAPTURE              #
#                          FOR WOLF DENSITY ESTIMATE                           #
#    Jose Jimenez, Daniel Cara, Francisco Garcia-Dominguez & Jose Barasona     # 
#                            01/12/2022 18:39:54                               #
#                                                                              #
#==============================================================================#

## Load libraries 
library(coda)
library(scrbook)
library(secr)
library(lattice)
library(mcmcOutput)
library(raster)
library(terra)
library(fields)
library(ggplot2)
library(nimble)

## Preparing data
#=================
setwd('C:/...')
source("Spiderplot_SCR.R")

# Operation mask. 119 days (occasions)
Oper<-data.matrix(read.table("Oper2.txt", header=FALSE))

wolf.ch <- secr::read.capthist("capt2.txt", "traps2.txt", detector='count', noccasions=119)
summary(wolf.ch)

traplocs<-as.matrix(secr::traps(wolf.ch))
X<-data.matrix(traplocs)/1000  # Scale
medX<- mean(X[,1])
medY<- mean(X[,2])
X[,1]<-(X[,1]-medX)    # Center
X[,2]<-(X[,2]-medY)

wolf<-as.array(wolf.ch)
wolf<- aperm(wolf,c(1,3,2))
str(wolf)

datYknown<- wolf

cat("Camera trapping days=", sum(Oper), "\n") # camera trap-day
# Mean days of capture per camera
cat("Mean trapping days per camera=", mean(apply((Oper),1,sum)), "\n")

K<-dim(Oper)[2] ;K     # sampling occasions
J<-dim(Oper)[1] ;J     # number of camera traps locations  

# Operation mask plot
x<-as.matrix(Oper)
image(1:K,1:J,t(x), yaxt = "n", xlab="Occasion", ylab="", cex.lab=1.25)
mtext(side = 2, "Camera", line = 2.5, cex=1.25)
axis(2, rev(seq(1, J, by=2)))
box()

# Summarize operation mask by occasions
KT<-apply(Oper,1,sum)

M<-80 # data augmentation
nind<- dim(datYknown)[1]; nind     # observed individuals
Yaug <- array(0, dim = c(M, J, K)) # augmented data

Yaug[1:nind, , ] <- datYknown
y<-apply(Yaug, c(1,2), sum)  # Summarized by traps

# State space
trapShape<-vect("C:/.../dataGIS/traps.shp")
buff_trap<-buffer(trapShape, width = 10447.71)
buffTrap<-aggregate(buff_trap)
plot(buffTrap)
points(trapShape)
area<-expanse(buffTrap)/1e8  # to transform ind/100 km2

(buff<- 10447.71/1000)   # 3*sigma scaled by 1000, like X
xl<-min(X[,1])-buff
xu<-max(X[,1])+buff
yl<-min(X[,2])-buff
yu<-max(X[,2])+buff
xlim=c(xl, xu)
ylim=c(yl, yu)

# Capture plot
plot(traplocs, pch="+", col="blue", 
  xlim=c(min(traplocs[,1])-2000,max(traplocs[,1])+2000), 
  ylim=c(min(traplocs[,2])-2000,max(traplocs[,2])+2000), type="n", asp=TRUE,
  xlab="UTM X", ylab="UTM Y", font.lab=2)
tot<-apply(datYknown, 2,sum)
symbols(traplocs, circles=tot*75, inches=F,bg="#228B2219", fg=NULL, add=T)
points(traplocs, pch="+", col="blue")
# Spiderplot
spiderplotJJ(datYknown, traplocs, buffer=2000, lwd=1)



library(nimble)
## define the model
code <- nimbleCode({
  
  alpha1 ~ dnorm(0,.01)
  sigma<- sqrt(1/(2*alpha1))
  psi ~ dbeta(1,1)
  # Random effects hyperparameter
  # Hyperparameters
  sigma.p ~ dunif(0, 10)
  mu0 ~ dnorm(0, 0.01)
  # Baseline detection rate
  for(j in 1:M){
    lp[i] ~ dnorm(mu0, sd=sigma.p)
    log(lam0[i]) <- lp[i]
  }
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    d2[i,1:J]<- (s[i,1] - X[1:J,1])^2 + (s[i,2] - X[1:J,2])^2
    lam[i,1:J]<- lam0[i]*exp(- alpha1*d2[i,1:J])
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

## Bundle constants
str(constants<-list(M=M, 
                    J=J, 
                    A=area))

## Bundle data
str(data   <-  list(Y=y, 
                    X=X, 
                    KT=KT,
                    zeros=rep(0,M), 
                    buffer=buff,  
                    xlim=xlim, 
                    ylim=ylim))

## Initial values
str(inits  <-  list(alpha1=runif(1,0.01,0.06),
                    psi=runif(1,0,1),
                    sigma.p=runif(1,0,3),
                    mu0=runif(1,-5,-1),
                    Ysim=y,
                    s=X[sample(1:J, size=M, replace=TRUE),],
                    z=c(rep(1, nind), rbinom((M-nind),1,0.2))))

## Create NIMBLE Model
Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
Rmodel$initializeInfo()
# initial values for "complex" quantities (by O. Gimenez in
# https://gist.github.com/oliviergimenez/e41e9cb99174f2124f948308e19ca7ec)
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

params<-c('N','D','sigma','mu0','sigma.p','psi','T1obs','T1sim','T2obs','T2sim','T3obs','T3sim')
mcmc<-configureMCMC(Rmodel,  monitors=params, enableWAIC = TRUE, thin = 5) 
# Rebuild and compile with new sampler
mcmc$removeSamplers("s")
ACnodes <- paste0("s[", 1:constants$M, ", 1:2]")
for(node in ACnodes) {
  mcmc$addSampler(target = node,
                  type = "RW_block",
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
}

MCMC <- buildMCMC(mcmc) 

CompMCMC <- compileNimble(MCMC, project = Rmodel)

## MCMC settings
nb<-5000*5
ni<-25000*5 +nb*5
nc<-3

outNim <- runMCMC(CompMCMC, niter = ni , nburnin = nb , inits=inits, 
                  nchains = nc, WAIC=TRUE, setSeed = FALSE, 
                  progressBar = TRUE, samplesAsCodaMCMC = TRUE)

## Results
outNim$WAIC
summary(outNim$samples)
