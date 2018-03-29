#simulations for 'Covariate adjustment and estimation of mean response in randomised trials'
#Jonathan Bartlett

library(MASS)
library(sandwich)
library(blockrand)

expit <- function(x) {
  exp(x)/(1+exp(x))
}

#define function which runs simulations

runNegBinSim <- function(nSim, n, b0, bx, bz, bxz, gammadisp, nbwm,randtype) {
  muhat1 <- array(0, dim=c(nSim,2))
  varest1 <- array(0, dim=c(nSim,3))

  muhat2 <- array(0, dim=c(nSim))
  varest2 <- array(0, dim=c(nSim,2))

  muhat3 <- array(0, dim=c(nSim))
  varest3 <- array(0, dim=c(nSim))

  for (sim in 1:nSim) {
    print(sim)

    x <- 1*(runif(n)<0.5)
    if (randtype==1) {
      #simple randomisation
      z <- 1*(runif(n)<0.5)
    } else if (randtype==2) {
      #random permuted block randomisation
      rand <- blockrand(n=n, num.levels=2)
      z <- as.numeric(rand$treatment)[1:n]-1
    } else if (randtype==3) {
      #stratified random permuted block
      z <- x
      n0 <- sum(x==0)
      z[x==0] <- as.numeric(blockrand(n=n0, num.levels=2)$treatment)[1:n0]-1
      n1 <- sum(x==1)
      z[x==1] <- as.numeric(blockrand(n=n1, num.levels=2)$treatment)[1:n1]-1
    }

    xb <- b0+bx*(x-0.5)+bz*(z-0.5)+bxz*(x-0.5)*(z-0.5)
    #specify shape of gamma
    k <- 2
    if (gammadisp==T) {
      #shape k
      gamma <- rgamma(n, k, k)
    } else {
      #we use log normal which has mean 1 and variance matching the gamma variance
      gamma <- exp(rnorm(n,-0.5*log(1+1/k),(log(1+1/k))^0.5))
    }
    
    #simulate follow-up time
    t <- rep(1,n)
    shortFup <- 1*(runif(n)<0.25)
    t[shortFup] <- runif(sum(shortFup))
    
    y <- rpois(n, t*gamma*exp(xb))

    simdata <- data.frame(y,x,z)

    if (nbwm==T) {
      mod <- glm.nb(y~x+z+offset(log(t)))
    } else {
      mod <- glm(y~x+z+offset(log(t)), family="poisson")
    }

    #calculate predictions for z=1
    pihat <- mean(z)

    #muhat1
    muhat1[sim,1] <- sum(y[z==1])/sum(t[z==1])
    #standard SE for this
    varest1[sim,1] <- (pihat^(-2))*(mean(t[z==1]))^(-2)*mean((z*(y-muhat1[sim,1]))^2)/n
    
    #muhat2
    activePredxb <- cbind(rep(1,n),x,rep(1,n)) %*% coef(mod)
    activePred <- exp(activePredxb)
    muhat2[sim] <- mean(activePred)
    #calculate delta method variances
    outer <- colMeans(as.numeric(exp(activePredxb)) * cbind(1,x,rep(1,n)))
    outer <- array(outer, dim=c(3,1))
    #sandwich variance treating covariates as fixed
    varest2[sim,1] <- t(outer) %*% vcovHC(mod, type = "HC0") %*% outer
    #variance estimator accounting random X assuming model correct
    varest2[sim,2] <- varest2[sim,1] + sum((activePred-muhat2[sim])^2)/(n^2)

    #muhat3
    muhat3[sim] <- muhat1[sim,1] - mean(((z-pihat)/pihat)*activePred)
    #variance of muhat3
    tauz <- mean(t[z==1])
    varest3[sim] <- (pihat^(-2))*tauz^(-2)*mean((z*(y-muhat3[sim]) -tauz*(z-pihat)*(activePred-muhat2[sim]) )^2)/n

  }

  list(muhat1=muhat1, muhat2=muhat2, muhat3=muhat3, varest1=varest1, varest2=varest2, varest3=varest3)

}

set.seed(698123)

nSim <- 10000
n <- 400

simRuns <- vector("list", 4)

simRuns[[1]] <- runNegBinSim(nSim=nSim,n=n,b0=0,bx=3,bz=1,bxz=0,gammadisp=T,nbwm=T,randtype=2)
#log normal dispersion parameter, but conditional mean model correct
simRuns[[2]] <- runNegBinSim(nSim=nSim,n=n,b0=0,bx=3,bz=1,bxz=0,gammadisp=F,nbwm=T,randtype=2)
#now conditional mean model wrong, gamma distributed dispersion
simRuns[[3]] <- runNegBinSim(nSim=nSim,n=n,b0=0,bx=3,bz=1,bxz=-1.5,gammadisp=T,nbwm=T,randtype=2)
#now use a Poisson working model
simRuns[[4]] <- runNegBinSim(nSim=nSim,n=n,b0=0,bx=3,bz=1,bxz=-1.5,gammadisp=T,nbwm=F,randtype=2)

resultPrint <- function(resultSet) {
  c(mean(resultSet$muhat1[,1]), 
    100*mean((resultSet$muhat1[,1]*exp(-1.96*resultSet$varest1[,1]^0.5 / resultSet$muhat1[,1])<mean(resultSet$muhat1[,1])) &
      (resultSet$muhat1[,1]*exp(1.96*resultSet$varest1[,1]^0.5 / resultSet$muhat1[,1])>mean(resultSet$muhat1[,1]))),
    mean(resultSet$muhat2)-mean(resultSet$muhat1[,1]),
    var(resultSet$muhat1[,1])/var(resultSet$muhat2),
    100*mean((resultSet$muhat2*exp(-1.96*resultSet$varest2[,1]^0.5 / resultSet$muhat2)<mean(resultSet$muhat1[,1])) &
           (resultSet$muhat2*exp(1.96*resultSet$varest2[,1]^0.5 / resultSet$muhat2)>mean(resultSet$muhat1[,1]))),
    100*mean((resultSet$muhat2*exp(-1.96*resultSet$varest2[,2]^0.5 / resultSet$muhat2)<mean(resultSet$muhat1[,1])) &
           (resultSet$muhat2*exp(1.96*resultSet$varest2[,2]^0.5 / resultSet$muhat2)>mean(resultSet$muhat1[,1]))),
    mean(resultSet$muhat3)-mean(resultSet$muhat1[,1]), 
    var(resultSet$muhat1[,1])/var(resultSet$muhat3),
    100*mean((resultSet$muhat3*exp(-1.96*resultSet$varest3^0.5 / resultSet$muhat3)<mean(resultSet$muhat1[,1])) &
           (resultSet$muhat3*exp(1.96*resultSet$varest3^0.5 / resultSet$muhat3)>mean(resultSet$muhat1[,1]))))
}

resTable <- matrix(unlist(lapply(simRuns, resultPrint)), byrow=F, ncol=length(simRuns))
colnames(resTable) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4")
rownames(resTable) <- c("Mean $\\hat{\\mu}_{1}(1)$", "95\\% CI Cov.",
                         "Bias $\\hat{\\mu}_{2}(1)$", "Rel. eff.", 
                          "Fixed $X$ 95\\% CI Cov.","Random $X$ 95\\% CI Cov.",
                         "Bias $\\hat{\\mu}_{3}(1)$", "Rel. eff.", "95\\% CI Cov.")
 
library(xtable)
print(xtable(data.frame(row = rownames(resTable),data.frame(resTable)), digits=2),sanitize.text.function=function(x){x}, include.rownames = F)


#stratified randomisation

simRuns2 <- vector("list", 4)

simRuns2[[1]] <- runNegBinSim(nSim=nSim,n=n,b0=0,bx=3,bz=1,bxz=0,gammadisp=T,nbwm=T,randtype=3)
#log normal dispersion parameter, but conditional mean model correct
simRuns2[[2]] <- runNegBinSim(nSim=nSim,n=n,b0=0,bx=3,bz=1,bxz=0,gammadisp=F,nbwm=T,randtype=3)
#now conditional mean model wrong, gamma distributed dispersion
simRuns2[[3]] <- runNegBinSim(nSim=nSim,n=n,b0=0,bx=3,bz=1,bxz=-1.5,gammadisp=T,nbwm=T,randtype=3)
#now use a Poisson working model
simRuns2[[4]] <- runNegBinSim(nSim=nSim,n=n,b0=0,bx=3,bz=1,bxz=-1.5,gammadisp=T,nbwm=F,randtype=3)

resultPrint <- function(resultSet) {
  c(mean(resultSet$muhat1[,1]), 
    100*mean((resultSet$muhat1[,1]*exp(-1.96*resultSet$varest1[,1]^0.5 / resultSet$muhat1[,1])<mean(resultSet$muhat1[,1])) &
               (resultSet$muhat1[,1]*exp(1.96*resultSet$varest1[,1]^0.5 / resultSet$muhat1[,1])>mean(resultSet$muhat1[,1]))),
    mean(resultSet$muhat2)-mean(resultSet$muhat1[,1]),
    var(resultSet$muhat1[,1])/var(resultSet$muhat2),
    100*mean((resultSet$muhat2*exp(-1.96*resultSet$varest2[,1]^0.5 / resultSet$muhat2)<mean(resultSet$muhat1[,1])) &
               (resultSet$muhat2*exp(1.96*resultSet$varest2[,1]^0.5 / resultSet$muhat2)>mean(resultSet$muhat1[,1]))),
    100*mean((resultSet$muhat2*exp(-1.96*resultSet$varest2[,2]^0.5 / resultSet$muhat2)<mean(resultSet$muhat1[,1])) &
               (resultSet$muhat2*exp(1.96*resultSet$varest2[,2]^0.5 / resultSet$muhat2)>mean(resultSet$muhat1[,1]))),
    mean(resultSet$muhat3)-mean(resultSet$muhat1[,1]), 
    var(resultSet$muhat1[,1])/var(resultSet$muhat3),
    100*mean((resultSet$muhat3*exp(-1.96*resultSet$varest3^0.5 / resultSet$muhat3)<mean(resultSet$muhat1[,1])) &
               (resultSet$muhat3*exp(1.96*resultSet$varest3^0.5 / resultSet$muhat3)>mean(resultSet$muhat1[,1]))))
}

resTable <- matrix(unlist(lapply(simRuns2, resultPrint)), byrow=F, ncol=length(simRuns2))
colnames(resTable) <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4")
rownames(resTable) <- c("Mean $\\hat{\\mu}_{1}(1)$", "95\\% CI Cov.",
                        "Bias $\\hat{\\mu}_{2}(1)$", "Rel. eff.", 
                        "Fixed $X$ 95\\% CI Cov.","Random $X$ 95\\% CI Cov.",
                        "Bias $\\hat{\\mu}_{3}(1)$", "Rel. eff.", "95\\% CI Cov.")

print(xtable(data.frame(row = rownames(resTable),data.frame(resTable)), digits=2),sanitize.text.function=function(x){x}, include.rownames = F)
