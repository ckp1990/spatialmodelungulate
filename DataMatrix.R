###R SCRIPT – PART 2 ###
###DataMatrix.R###
###This script is to organize data in a matrix form required for the analysis###
### Read all input files ###
x2=read.table(paste0(dstf),header=TRUE,sep=',') # line transect file
TRid=read.csv(paste0(trf)) # transect ID
modelgrid=maptools:::read.shape(paste0(shp)) # shape file with grids
modelpoly=sf::st_read(paste0(shp)) # shape file with grids
covsites=read.csv(paste0(trcov))
### Subset distance data for species that is being analyzed ###
x=subset(x2,x2$walk.no!='NA')
y=subset(x,x$species==paste0(sps))
l=length(x$tr.no)
ntrans=max(x$tr.no)
### save covariates for use with more complicated models ###
cov=rep(0,ntrans); tottrlen=cov; nwalks2=cov;
for (i in 1:l) { j=x$tr.no[i]
tottrlen[j]=x$total.tr.length[i]
if (x$walk.no[i]>nwalks2[j]) nwalks2[j]=x$walk.no[i]
}
### compute cluster size category for each distance category for each transect ###
### For each transect, you must have cluster size as rows, distance category as columns ###
l=length(y$tr.no); t5=floor(y$p.dist/distcatsize)+1; ndistcat=max(t5);
hold=floor(y$gr.sz/grszcat)+ 1; ngs=max(hold);
newy2=array(0,dim=c(ngs,ndistcat,ntrans));
for (i in 1:l) { k=y$tr.no[i]
grsz=y$gr.sz[i]
if (grsz>grszlimit) grsz=grszlimit
gs=floor(grsz/grszcat)+ 1
dist=y$p.dist[i]
if (dist>distlimit) dist=distlimit
distcat=floor(dist/distcatsize)+1
newy2[gs,distcat,k]=newy2[gs,distcat,k]+1
}
distBreaks = seq(distcatsize, ndistcat*distcatsize, by=distcatsize) # endpoints of distance categories
grszBreaks = seq(grszcat, ngs*grszcat, by=grszcat) # endpoints of cluster size categories
grszMean = seq(sum(1:grszcat)/grszcat, by=grszcat, length=ngs) # means of cluster size categories
logFactorial = lgamma((1:grszBreaks[ngs]) + 1)
### remove transects which are OUTSIDE the park boundary, hence not sampled ###
newy=newy2[,,-subsetlist]
## Now, dim(newy) has 13 group size classes, 19 distance classes, 77 transects##
nwalks=nwalks2[-subsetlist]
ntrans=length(nwalks)
ndistwalk=tottrlen[-subsetlist]
### Generate adjacency matrix ###
poly <- modelgrid$att.data#maptools:::Map2poly(modelgrid)
adjmatrix=poly2nb(modelpoly)
### Scaling and centering of landscape level covariates ###
cov1km=poly
cov1km <- cov1km %>% as.data.frame() ##new 
#cov1km$LineID=1:length(modelpoly)
## Brut force code to turn PCHSCR into NUM
#cov1km$PCHSCR<-as.numeric(cov1km$PCHSCR)
## Brut force code to scale the data

cov1km[names(cov1km) %in% sb2] <- lapply(cov1km[names(cov1km) %in% sb2], 
                                         function(x) c(scale(x)))
### Scaling and centering of transect level covariates ###
covsites1<-covsites[c(names(covsites) %in% sb1)]

covsites1[names(covsites1) %in% sb1] <- lapply(covsites1[,names(covsites1) %in% sb1], function(x) c(scale(x)))
covsites1$tr.no<-covsites$tr.no
covsites=covsites1
rm(covsites1)
### Set up the Data matrix with grid Id and proportion of transect in that grid ###
TRid$Tr.no[c(as.numeric(paste0(wrid)))]=rrid
TRid$uqtemp=paste(TRid$Tr.no, TRid$Nummer, sep="a")
a=unique(TRid$uqtemp)
newTR=matrix(0, ncol=4, nrow=length(a))
newTR[,1]=a
for(i in 1:length(a)) {
  b=which(TRid$uqtemp == a[i])
  newTR[i,2]=TRid$Nummer[b[1]]
  newTR[i,3]=sum(TRid$Length[b])
  newTR[i,4]=TRid$Tr.no[b[1]]
}
newTR=as.data.frame(newTR)
newTR[,3]=as.numeric(as.character(newTR[,3]))
newTR[,2]=as.numeric(as.character(newTR[,2]))
newTR[,4]=as.numeric(as.character(newTR[,4]))
newTR[,5]=pmatch(newTR[ ,2], cov1km$NUMMER, duplicates=TRUE)
b=newTR[,5]=='NA'
c=which(b=="FALSE")
newTR[-c,5]=1
uqid=1:uqid
bigM=matrix(0, ncol=as.numeric(igrid), nrow=length(uqid))
bigM=as.data.frame(bigM)
for (i in 1:length(uqid)){
  bigM[i,1]=uqid[i]
  a=which(newTR[,4] == uqid[i])
  b=length(a)
  for(j in 2:(b+1)) {
    bigM[i,j]=newTR[a[j-1],5]
  }
  for(j in 7:(b+6)) {
    bigM[i,j]=newTR[a[j-6],3]
  }
}
### create proportions from length and return to bigM ###
hold7=bigM[,7]/(bigM[,7]+bigM[,8]+bigM[,9]+bigM[,10]+bigM[,11])
hold8=bigM[,8]/(bigM[,7]+bigM[,8]+bigM[,9]+bigM[,10]+bigM[,11])
hold9=bigM[,9]/(bigM[,7]+bigM[,8]+bigM[,9]+bigM[,10]+bigM[,11])
hold10=bigM[,10]/(bigM[,7]+bigM[,8]+bigM[,9]+bigM[,10]+bigM[,11])
hold11=bigM[,11]/(bigM[,7]+bigM[,8]+bigM[,9]+bigM[,10]+bigM[,11])
bigM[,7]=hold7
bigM[,8]=hold8
bigM[,9]=hold9
bigM[,10]=hold10
bigM[,11]=hold11
##subset out those transects NOT walked ##
bigM2=bigM[-subsetlist,]
for(i in 2:6){
  a=which(bigM2[,i] == 0)
  bigM2[a,i]=1
}

