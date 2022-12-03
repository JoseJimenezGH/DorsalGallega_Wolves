#==============================================================================#
#                                                                              #
#                 VIDEO CAMERA TRAP AND SPATIAL CAPTURE-RECAPTURE              #
#                          FOR WOLF DENSITY ESTIMATE                           #
#    José Jiménez, Daniel Cara, Francisco García-Dominguez & Jose Barasona     # 
#                            01/12/2022 18:39:54                               #
#                                                                              #
#==============================================================================#

## Load libraries 
library(coda)
library(scrbook)
library(lattice)
library(raster)
library(rgdal)
library(fields)
library(nimble)

## Preparing data
#=================
setwd('C:/Users/Usuario/OneDrive/30 Proyecto Lobo Suido/05 R code/Suido/DEF')
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
trapShape<-readOGR(dsn="C:/Users/Usuario/OneDrive/30 Proyecto Lobo Suido/01 data/GIS",layer="traps")
buffTrap<-buffer(trapShape, width = 9947.409)
plot(buffTrap)
points(trapShape)
A<-area(buffTrap)/1e8  # Aqui sale en ud/100 km2

(buff<- 9947.409/1000)   # 3*sigma scaled by 1000, like X
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
  alpha1 ~ dnorm(0,.1)
  sigma<- sqrt(1/(2*alpha1))
  psi ~ dbeta(1,1)  # data augmentation parameter
  r ~ dgamma(.1,.1)
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

    for(j in 1:J){   
      d2[i,j]<- (s[i,1] - X[j,1])^2 + (s[i,2] - X[j,2])^2
      lam[i,j]<- p0[j]*exp(- alpha1*d2[i,j])
      mu[i,j] <- lam[i,j]*z[i]
      mu2[i,j] <- lam[i,j]*z[i]*KT[j]
      p[i,j] <- r/(r+mu2[i,j])
      Y[i,j] ~ dnegbin(p[i,j],r)
      outj[i,j] <- sqrt(d2[i,j]) > buffer  # zero-trick by R. Chandler in
      # https://groups.google.com/g/spatialcapturerecapture/c/NzqUovn8jF0/m/Plg2g6O6AgAJ
      
      # Components for fit diagnostics
      psim[i,j] <- r/(r+mu[i,j]*KT[j])
      Ysim[i,j] ~ dnegbin(psim[i,j],r)       # simulated
      Yexp[i,j] <- lam[i,j] * z[i] * KT[j]   # expected
      
      # Components for T1
      err1obs[i,j] <- (sqrt(Y[i,j]) - sqrt(Yexp[i,j]))^2
      err1sim[i,j] <- (sqrt(Ysim[i,j]) - sqrt(Yexp[i,j]))^2
    }
    # Components for T2
    err2obs[i] <- (sqrt(sum(Y[i,1:J])) - sqrt(sum(Yexp[i,1:J])))^2
    err2sim[i] <- (sqrt(sum(Ysim[i,1:J])) - sqrt(sum(Yexp[i,1:J])))^2

    out[i] <- equals(sum(outj[i,1:J]),J)
    zeros[i] ~ dbern(out[i])    
  }

  # Components for T3
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
  ## Number of individuals in the population
  N<-sum(z[1:M])
  D<-N/A
})


## Bundle constants
str(constants<-list(M=M,
                    J=J, 
                    A=A))
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
                    Ysim=y,
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
Rmodel$lp
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)
## Parameters monitored
params<-c('N','D','sigma','mu0','sigma.p','psi','r',
          'T1obs','T1sim','T2obs','T2sim','T3obs','T3sim')
mcmc<-configureMCMC(Rmodel,  monitors=params, thin = 5, enableWAIC = TRUE)

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
nb<-1000     # 5000*5
ni<-5000+nb  # 25000*5 +nb*5
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


# Assessing fit of the model
#========================================  
samples<-as.matrix(outNim$samples[,c('T1obs','T1sim','T2obs','T2sim','T3obs','T3sim')])
par(mfrow=c(2,2),
    mar=c(6,5,4,4), las=1)  
options(scipen=999)
# Fit statistic 1: individual x trap frequencies.
pvalue<-round(mean(samples[,'T1sim'] > samples[,'T1obs']), 2)
plot(samples[,'T1sim'], samples[,'T1obs'], xlim=c(100,300), ylim=c(100,300),
  xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", 
  bty ="n", las=1, col="grey50", main="Individual x trap\n frequencies", type='n',
  sub=paste("p-value=", pvalue,""))
idx1 = samples[,'T1sim']-samples[,'T1obs']>0
idx2 = samples[,'T1sim']-samples[,'T1obs']<0   
points(samples[,'T1sim'][idx2],samples[,'T1obs'][idx2], col='#00000033', pch=1, cex=.15)
points(samples[,'T1sim'][idx1],samples[,'T1obs'][idx1], col='#FF000033', pch=1, cex=.15)
segments(100,100,300,300, lwd=2)

# Fit statistic 2: Individual encounter frequencies.
pvalue<-round(mean(samples[,'T2sim'] > samples[,'T2obs']), 2)
plot(samples[,'T2sim'], samples[,'T2obs'], xlim=c(0,60), ylim=c(0,60), 
  xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", 
  bty ="n", las=1, col="grey50", main="Individual encounter\n frequencies", type='n',
  sub=paste("p-value=", pvalue,""))
idx1 = samples[,'T2sim']-samples[,'T2obs']>0
idx2 = samples[,'T2sim']-samples[,'T2obs']<0   
points(samples[,'T2sim'][idx2],samples[,'T2obs'][idx2], col='#00000033', pch=1, cex=.15) 
points(samples[,'T2sim'][idx1],samples[,'T2obs'][idx1], col='#FF000033', pch=1, cex=.15)
segments(0,0,60,60, lwd=2)

# Fit statistic 3: Trap frequencies.
pvalue<-round(mean(samples[,'T3sim'] > samples[,'T3obs']), 2)
plot(samples[,'T3sim'], samples[,'T3obs'], xlim=c(5,60), ylim=c(5,60), 
  xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", 
  bty ="n", las=1, col="grey50", main="Trap encounter\n frequencies", type='n',
  sub=paste("p-value=", pvalue,""))
idx1 = samples[,'T3sim']-samples[,'T3obs']>0
idx2 = samples[,'T3sim']-samples[,'T3obs']<0   
points(samples[,'T3sim'][idx2],samples[,'T3obs'][idx2], col='#00000033', pch=1, cex=.15)
points(samples[,'T3sim'][idx1],samples[,'T3obs'][idx1], col='#FF000033', pch=1, cex=.15)
segments(0,0,60,60, lwd=2)

