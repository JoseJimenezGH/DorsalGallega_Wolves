#==============================================================================#
#                                                                              #
#                 VIDEO CAMERA TRAP AND SPATIAL CAPTURE-RECAPTURE              #
#                          FOR WOLF DENSITY ESTIMATE                           #
#           Negative binomial Model + random effects - Faster version          #
#      Jose Jimenez, Daniel Cara, Francisco Garcia­Dominguez & Jose Barasona   # 
#                            01/03/2023 17:13:21                               #
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
Oper<-data.matrix(read.table("Oper.txt", header=FALSE))

wolf.ch <- secr::read.capthist("capt.txt", "traps.txt", detector='count', noccasions=119)
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
trapShape<-vect(".../dataGIS/traps.shp")
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


## Define nimble model
code <- nimbleCode({
  
  ## Scale for the half normal detection function
  alpha1 ~ dnorm(0,.01)
  sigma<- sqrt(1/(2*alpha1))
  psi ~ dbeta(1,1)  # data augmentation parameter
  r ~ dgamma(.01,.01)
  # Hyperparameters
  sigma.p ~ dunif(0, 10)
  mu0 ~ dnorm(0, 0.01)
  # Baseline detection rate
  for(j in 1:J){
    lp[j] ~ dnorm(mu0, sd=sigma.p)
    log(p0[j]) <- lp[j]
  }
  
  for(i in 1:M){
    z[i] ~ dbern(psi)
    ## Activity centres of the observed individuals
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    d2[i,1:J]<- (s[i,1] - X[1:J,1])^2 + (s[i,2] - X[1:J,2])^2
    #lam[i,1:J]<- p0[1:J]*exp(- alpha1*d2[i,1:J])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], 
                                   X = X[1:J,1:2], 
                                   J=J,
                                   sigma=sigma, 
                                   p0=p0[1:J], 
                                   z=z[i])

    p[i,1:J] <- r/(r+lam[i,1:J])
    Y[i,1:J] ~ dNBVector(p=p[i,1:J],r=r*KT[1:J],z=z[i])
    
    for(j in 1:J){   
      outj[i,j] <- sqrt(d2[i,j]) > buffer  # zero-trick by R. Chandler in
      # https://groups.google.com/g/spatialcapturerecapture/c/NzqUovn8jF0/m/Plg2g6O6AgAJ
    }
    out[i] <- equals(sum(outj[i,1:J]),J)
    zeros[i] ~ dbern(out[i])    
  }

  ## Number of individuals in the population
  N<-sum(z[1:M])
  D<-N/A
})

GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), p0=double(1), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     ans <- p0*exp(-d2/(2*sigma^2))
     return(ans)
    }
  }
)

dNBVector <- nimbleFunction(
  run = function(x = double(1), p = double(1), r = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dnbinom(x, p = p, size = r, log = TRUE))
      return(logProb)
    }
  }
)

rNBVector <- nimbleFunction(
  run = function(n = integer(0),p = double(1),r = double(1), z = double(0)) {
    returnType(double(1))
    J <- length(p)
    ans<- numeric(J, value=0)
    return(ans)
  }
)

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
str(inits  <-  list(alpha1=runif(1,0,0.1),
                    psi=runif(1,0,1),
                    sigma.p=runif(1,0,3),
                    mu0=runif(1,-5,-1),
                    r=runif(1,0,4),
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
#Rmodel$lp
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
## Parameters monitored
params<-c('N','D','sigma','mu0','sigma.p','psi','r')
mcmc<-configureMCMC(Rmodel,  monitors=params, thin = 1, enableWAIC = TRUE)

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
nb<-10000
ni<-50000 +nb
nc<-3

## Run MCMC
outNim <- runMCMC(CompMCMC, niter = ni , nburnin = nb , inits=inits, 
                  nchains = nc, WAIC=TRUE, setSeed = FALSE, 
                  progressBar = TRUE, samplesAsCodaMCMC = TRUE)

## Results
outNim$WAIC
summary(outNim$samples)

# We examine traceplots
xyplot(outNim$samples[,c('N','D','mu0','sigma','sigma.p','psi','r')])

# N histogram
hist(as.matrix(outNim$samples[,'N']),xlim=c(25,50), xlab="Population size", main="")
mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
abline(v=mode(as.matrix(outNim$samples[,'N'])), lty=2, lwd=3, col="blue")

