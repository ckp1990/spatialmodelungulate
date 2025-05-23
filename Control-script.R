### R SCRIPT – PART 1 ###3
### ControlScript.R###
### Load required libraries and scripts ###
library(nimble)
library(coda)
library(spdep)
library(maptools)  ## not package called map tools
library(sf)
library(terra)
library(rgdal)
library(mcmcse)
library(matrixStats)
#setwd("..") 
### Choose species to be analysed ###
###The data from 5 species; CHT-chital, SBR-sambar, GAR-gaur, PIG-wild pig, MJK-muntjac were analyzed using these scripts###
sps='SBR'
### Specify input files ###
# line transect data file; see Appendix 1 for data structure #
dstf<-'Data/nhb.csv'
# transect details file; see Appendix 1 for data structure #
trf<-'Data/TRidentity.csv'
# Shape files with landscape level covariates; see Appendix 1 for data structure #
shp<-'Data/cov1km.shp'
# transect level covariates; see Appendix 1 for data structure #
trcov<-'Data/TRsbrcovariates1.csv'

### Specify values for estimating detection function ###
distcatsize=20; # distance category size; 15 to 20 categories are ideal#
distlimit=212; # maximum distance observed#
grszcat=9; # cluster category size; 10 to 15 categories are ideal#
grszlimit=9; # maximum cluster size observed#
### Specify sampling design ###
## Number of transects in the study design ##
uqid=83
# use the maximum number of grids intersected by any transect to arrive at this number; edit lines in the DataMatrix.R script based on sampling design details #
igrid=11
# Clean and remove transects which are not sampled, e.g., outside the park administrative boundary, and hence not walked in our case #
# List of transects outside the study area; these are specific to the data being analyzed #
subsetlist=c(40,54,61,67,72,78) ## looks like hard coding to me. 

## List of rows with wrong transect IDs ##
wrid<-c(58,59,60, 65,66,67)
## Assign correct ID to the above##
rrid=40
## Specify columns with transect number and transect level covariates; these depend on the data being analyzed ##
sb1<-c("tr.no","SBRpalplants.m2","thabitat.disturb.km")
## Specify columns with landscape level covariates; these depend on the data being analyzed ##
sb2<-c("VARSLOPE","WTRAVGDIST","ECODISTAVG","PCHSCR")
### Process data ###
source('DataMatrix.R')
### model definition ###
Ngrid= 1792 ## brut force
sumNumNeigh=length(unlist(adjmatrix))
adj=unlist(adjmatrix)
num=sapply(adjmatrix, length)
bigM=as.matrix(bigM2)
cov1km=as.matrix(cov1km)
covsites=as.matrix(covsites)
## define constants in the model ##
const <- list(sumNumNeigh=sumNumNeigh, bigM=bigM,
              adj=adj, num=num,
              cov1km=cov1km,
              covsites=covsites, ntrans=ntrans,
              ndistcat=ndistcat, ngs=ngs, newy=newy,
              Ngrid=Ngrid, ndistwalk=ndistwalk, grszcat=grszcat,
              distBreaks=distBreaks, grszMean=grszMean, logFactorial=logFactorial)
## read nimble model script for the analysis WITHOUT indicator variables##
### edit lines in the 'NimModUniPrior.R' script based on the number of covariates used ###
source('NimModUniPrior.R')
## comment the above line and use the line below for the analysis WITH indicator variables, after uncommenting ##
## edit lines in the 'NimModUniPriorWITHindicators.R' script based on the number of covariates used##
#source('/ungulate/NimModUniPriorWITHindicators.R')#
## define parameters for the model WITHOUT indicator variables ##
##parameters <- c("b", "p", "lams", "sigs", "beta1", "beta2", "beta3", "beta4", "alpha1", "alpha2", "gs", "grszMean","sigma","sigma0")
## comment the above line and use the line below to define parameters for the model WITH indicator variables after uncommenting##
parameters <- c("b", "p", "lams", "sigs", "beta1", "beta2", "beta3", "beta4", "alpha1", "alpha2", "gs", "grszMean","sigma","sigma0", "wa1", "wa2", "w1", "w2", "w3", "w4")##
## dataset used in the monograph has two site level covariates (alpha1, alpha2) and four grid-cell level covariates (beta1, beta2, beta3, beta4); correspondingly, there are six indicator variables (wa1, wa2, w1, w2, w3, w4); update the number of indicator variables depending on the number of covariates used in the analysis##
## specify either random or plausible initial values for the model WITHOUT indicator variables##
a=1
b=as.vector(rnorm(Ngrid, 0, .5))
#inits <- function(){list(b=b, p=0.2, lams=5, sigs=1, sigma0=3, beta1=a, beta2=a, beta3=a, beta4=a, alpha1=a, alpha2=a)}
## comment above lines and use the lines below to specify either random or plausible initial values for the model WITH indicator variables after uncommenting ##
##a=1##
##b=as.vector(rnorm(Ngrid, 0, .5))##
inits <- function(){list(b=b, p=0.2, lams=5, sigs=1, sigma0=3, beta1=a, beta2=a, beta3=a, beta4=a, alpha1=a, alpha2=a, wa1=a, wa2=a, w1=a, w2=a, w3=a, w4=a)}##
### running nimble model ###
# use the model specification script code ‘NimModUniPrior.R’ for the model WITHOUT indicator variables; and, ‘NimModUniPriorWITHindicators.R’ for the model WITH indicator variables#
t1=Sys.time()
assign(paste0(sps,"UP"),nimbleMCMC(
  code = NimModUniPriorWITHindicators,
  constants = const,
  inits = inits,
  monitors = parameters,
  niter = 220000,
  nburnin = 20000,
  summary=T,
  thin = 1))
print(Sys.time()-t1)
### save objects as a backup to avoid any loss ###
### comment this if you want to save read/write time ###
save.image()
### Compute summaries of posterior distribution ###
## check object ‘grszBreaks’ to determine which ‘groupsizes’ should be plotted ##
## uncomment lines in 'posteriorSummariesWITHdetectionPlot.R' if you are using the model WITH indicator variables ##
grszSeq = seq(1:length(grszBreaks))
source('D:/spatialmodelungulate-MTL/posteriorSummariesWITHdetectionPlot.R')
## Compute local, site-level and landscape-level abundances and generate density surface map ##
## Update lines 6-7 in the 'cell_abundance.R' script depending on the number of covariates used in the analysis ##
source('cell_abundance.R')
save.image()
