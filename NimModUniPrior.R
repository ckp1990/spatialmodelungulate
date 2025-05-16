### R SCRIPT – PART 3 ###
###NimModUniPrior.R script###
###This is the model specification script for the model WITHOUT indicator variables###
### observed distance from transect and cluster size are detection covariates###
### alpha1 and alpha2 are site-level abundance covariates###
### beta1, beta2, beta3, beta4 are grid-cell level abundance covariates###
library(nimble)
NimModUniPrior<- nimbleCode({
  sigma0~dunif(-5,5)
  p~dunif(0,10)
  ## modeling sigma parameter to include cluster size effects##
  for(k in 1:ngs) {
    log(sigma[k])<-sigma0 + p*(grszMean[k]-1)
    sigma2[k]<-sigma[k]*sigma[k]
  }
  for (m in 1:(ngs*grszcat)) {
    summand[m] <- exp(m * log(lams) - logFactorial[m])
    summandForMean[m] <- m * summand[m]
  }
  ##modeling detection function using half-normal function##
  for(k in 1:ngs){
    mncell[k,1] <- (sqrt(2*3.1416) * sigma[k] / distBreaks[ndistcat]) *(phi(distBreaks[1]/sigma[k]) - 0.5)
    for(j in 2:ndistcat) {
      mncell[k,j] <- (sqrt(2*3.1416) * sigma[k] / distBreaks[ndistcat]) * (phi(distBreaks[j]/sigma[k]) - phi(distBreaks[j-1]/sigma[k]))
    }}
  ## modeling cluster size as a zero-truncated Poisson variable##
  gs[1] <- exp(-lams) / (1 - exp(-lams)) * sum(summand[1:grszcat])
  grszMean[1] <- sum(summandForMean[1:grszcat]) / sum(summand[1:grszcat])
  for (k in 2:ngs) {
    gs[k] <- exp(-lams) / (1 - exp(-lams)) * sum(summand[((k-1)*grszcat+1):(k*grszcat)])
    grszMean[k] <- sum(summandForMean[((k-1)*grszcat+1):(k*grszcat)]) / sum(summand[((k-1)*grszcat+1):(k*grszcat)])
  }
  lams~dunif(0,50)
  ## modeling site-level covariate effects, variable sampling efforts and spatial misalignment##
  for(i in 1:ntrans){
    lam[i]<- ndistwalk[i]*(bigM[i,7]*Z[bigM[i,2]] + bigM[i,8]*Z[bigM[i,3]] + bigM[i,9]*Z[bigM[i,4]]+bigM[i,10]*Z[bigM[i,5]] + bigM[i,11]*Z[bigM[i,6]]) * exp(alpha1*covsites[i,2] + alpha2*covsites[i,3])
    for (k in 1:ngs) {
      for(j in 1:ndistcat) {
        mu[k,j,i]<- lam[i]*mncell[k,j]*gs[k]
        newy[k,j,i]~dpois(mu[k,j,i])
      }}}
  ##modeling landscape-level covariate effects and random spatial effects##
  for(k in 1:Ngrid){
    log(Z[k])<-b[k] + beta1*cov1km[k,3] + beta2*cov1km[k,4] + beta3*cov1km[k,5] + beta4*cov1km[k,6]
  }
  ##specify Gaussian CAR prior##
  b[1:Ngrid] ~ car.normal(adj[1:13132], weights[1:13132], num[1:1792], tau)
  for(k in 1:sumNumNeigh) {
    weights[k] <- 1
  }
  ##specify prior distribution##
  alpha1~dunif(-5,5)
  alpha2~dunif(-5,5)
  beta1~dunif(-5,5)
  beta2~dunif(-5,5)
  beta3~dunif(-5,5)
  beta4~dunif(-5,5)
  sigs~dunif(0,5)
  tau<-1/(sigs*sigs)
})
###NimModUniPriorWITHindicators.R###
### observed distance from transect, cluster size are detection covariates###
### alpha1 and alpha2 are site-level abundance covariates###
### beta1, beta2, beta3, beta4 are grid-cell level abundance covariates###
### wa1, wa2, w1, w2, w3, w4 are indicator variables for the six covariates###
### UNIFORM PRIORS ###
library(nimble)
NimModUniPriorWITHindicators<- nimbleCode({
  sigma0~dunif(-5,5)
  p~dunif(0,10)
  ##modeling sigma parameter to include cluster size effects##
  for(k in 1:ngs) {
    log(sigma[k])<-sigma0 + p*(grszMean[k]-1)
    sigma2[k]<-sigma[k]*sigma[k]
  }
  for (m in 1:(ngs*grszcat)) {
    summand[m] <- exp(m * log(lams) - logFactorial[m])
    summandForMean[m] <- m * summand[m]
  }
  ##modeling detection function with half-normal function##
  for(k in 1:ngs){
    mncell[k,1] <- (sqrt(2*3.1416) * sigma[k] / distBreaks[ndistcat]) * (phi(distBreaks[1]/sigma[k]) - 0.5)
    for(j in 2:ndistcat) {
      mncell[k,j] <- (sqrt(2*3.1416) * sigma[k] / distBreaks[ndistcat]) * (phi(distBreaks[j]/sigma[k]) - phi(distBreaks[j-1]/sigma[k]))
    }}
  ##modeling cluster size as a zero-truncated Poisson variable##
  gs[1] <- exp(-lams) / (1 - exp(-lams)) * sum(summand[1:grszcat])
  grszMean[1] <- sum(summandForMean[1:grszcat]) / sum(summand[1:grszcat])
  for (k in 2:ngs) {
    gs[k] <- exp(-lams) / (1 - exp(-lams)) * sum(summand[((k-1)*grszcat+1):(k*grszcat)])
    grszMean[k] <- sum(summandForMean[((k-1)*grszcat+1):(k*grszcat)]) / sum(summand[((k-1)*grszcat+1):(k*grszcat)])
  }
  lams~dunif(0,50)
  ##modeling transect-level abundance, site-level covariate effects, indicator variables for each site-level covariate effect, variable sampling efforts and spatial misalignment##
  for(i in 1:ntrans){
    lam[i]<- ndistwalk[i]*(bigM[i,7]*Z[bigM[i,2]] + bigM[i,8]*Z[bigM[i,3]] + bigM[i,9]*Z[bigM[i,4]]+bigM[i,10]*Z[bigM[i,5]] + bigM[i,11]*Z[bigM[i,6]]) * exp(wa1*alpha1*covsites[i,2] + wa2*alpha2*covsites[i,3])
    for (k in 1:ngs) {
      for(j in 1:ndistcat) {
        mu[k,j,i]<- lam[i]*mncell[k,j]*gs[k]
        newy[k,j,i]~dpois(mu[k,j,i])
      }}}
  ##modeling grid-cell level covariate effects, indicator variables for each grid-cell level covariate effect and random spatial effects##
  for(k in 1:Ngrid){
    log(Z[k])<-b[k] + w1*beta1*cov1km[k,3] + w2*beta2*cov1km[k,4] + w3*beta3*cov1km[k,5] + w4*beta4*cov1km[k,6]
  }
  ##specify Gaussian CAR prior##
  b[1:Ngrid] ~ car.normal(adj[1:13132], weights[1:13132], num[1:1792], tau)
  for(k in 1:sumNumNeigh) {
    weights[k] <- 1
  }
  ##specify prior distributions##
  alpha1~dunif(-5,5)
  alpha2~dunif(-5,5)
  beta1~dunif(-5,5)
  beta2~dunif(-5,5)
  beta3~dunif(-5,5)
  beta4~dunif(-5,5)
  wa1~dbin(0.5,1)
  wa2~dbin(0.5,1)
  w1~dbin(0.5,1)
  w2~dbin(0.5,1)
  w3~dbin(0.5,1)
  w4~dbin(0.5,1)
  sigs~dunif(0,5)
  tau<-1/(sigs*sigs)
})
