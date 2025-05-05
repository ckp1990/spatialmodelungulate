### R SCRIPT – PART 4 ###
###posteriorSummariesWITHdetectionPlot.R###
### This script is for computing posterior summaries and their Monte Carlo standard errors, and #generating detection function plots ###
## Read MCMC samples and compute average group size ##
coda.samples <- as.mcmc(get(paste0(sps,"UP"))$samples)
mc= as.matrix(coda.samples)
mu = mc[,'lams']
aveGroupSize = mu / (1 - exp(-mu))
mc = cbind(mc, aveGroupSize)
## calculate Monte Carlo standard errors ##
post.stats.MCSE.mean<-mcse.mat(x = mc, method = "bm", g = NULL)
post.stats.MCSE.median<-mcse.mat(x = mc, method = "bm", g = median)
post.stats.MCSE.ql<-mcse.q.mat(x = mc, method = "bm", 0.025)
post.stats.MCSE.qu<-mcse.q.mat(x = mc, method = "bm", 0.975)
post.stats.hpd<-HPDinterval(as.mcmc(mc))
post.stats.SD<-as.matrix(apply(mc, 2, sd))
post.stats<-as.data.frame(cbind(post.stats.MCSE.mean,
                                cbind(post.stats.MCSE.median,
                                      cbind(post.stats.MCSE.ql,
                                            cbind(post.stats.MCSE.qu,
                                                  cbind(post.stats.hpd, post.stats.SD))))))
names(post.stats)<-c("Mean", "Mean.SE", "Median", "Median.SE", "2.5%", "2.5%.SE", "97.5%", "97.5%.SE", "95%_HPD_L", "95%_HPD_U", "SD")
## Export estimates and SE to file ##
write.csv(post.stats, paste0(sps,"_Monte_Carlo_estimates_SE_",
                             format(Sys.Date(), "%d%b%Y"),".csv"),row.names = T)
## Plot detection function for each group size category##
ngsCategories = length(grszMean)
grszMean.columns = paste0('grszMean', '[',1:ngsCategories,']')
grszM = mc[, grszMean.columns]
alpha0 = mc[,'sigma0']
alpha1 = mc[,'p']
sigm = exp(alpha0 + alpha1*(grszM -1))
B = distlimit+20
deltax = 5
x = seq(deltax, B, by=deltax)
probDetection = array(dim=c(nrow(sigm), ncol(sigm), length(x)))
for (k in 1:length(x)) {
  probDetection[,, k] = exp(-(x[k]*x[k])/(2*sigm*sigm))
}
grszColor = c('black','grey50','cyan','blue','blueviolet', 'green','palegreen4','khaki4','green4','red', 'firebrick4','rosybrown','goldenrod')
pdf(paste0(sps,'_DetectionPerGroup',format(Sys.Date(), "%d%b%Y"),'.pdf'), onefile=T)
for (j in grszSeq) {
  y = apply(probDetection[,j,], 2, mean)
  ylower = apply(probDetection[,j,], 2, quantile, prob=.025)
  yupper = apply(probDetection[,j,], 2, quantile, prob=.975)
  if (j==1) plot(x, y, type="l", cex.axis=1.4, cex.lab=1.6, las=1, lwd=2, xlab='Perpendicular distance from transect line (meters)', ylab='Probability of detection', col=grszColor[j])
  ind.rev = seq(length(x), 1, by=-1)
  lines(x, y, lwd=2, col=grszColor[j])
}
legend(x='topright', legend=as.character(round(apply(grszM, 2, mean)[grszSeq], digits=1)), fill=grszColor[seq(1:length((grszSeq)))]) ## fix this code for legend fill colour##
dev.off()
rm(i,j)
### Plotting posterior densities ###
csam<-as.mcmc(mc[,-(3:(Ngrid+2))])
png(filename = paste0(sps,'_densityPlots%03d',"_",format(Sys.Date(), "%d%b%Y"),'.png'), width = 8, height = 8, units = "cm", pointsize=6,res = 600)
densplot(csam, show.obs=F, ylab = "Density")
dev.off()
### Cross correlation plots ###
png(filename = paste0(sps,'_crosscorPlots',"_",format(Sys.Date(), "%d%b%Y"),'.png'), width = 8, height = 8, units = "cm", pointsize=6,res = 600)
crosscorr.plot(csam[,c(1:2,4:7)])
dev.off()
### uncomment the lines below, when using the ‘NimModUniWITHindicators.R’ script for assessing relative importance of covariates and generating posterior model weight summaries###
#pmw=paste(mc[,'wa1'],mc[,'wa2'],mc[,'w1'],mc[,'w2'],mc[,'w3'],mc[,'w4'])
#write.table(table(pmw), paste0(sps,'_pmw_',format(Sys.Date(), "%d%b%Y"),'.txt'), row.names = F)