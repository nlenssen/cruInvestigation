source("Code/namelist.Rnl")

## Parameters

# Papua New Guinea
# targetLon <- 142.8
# targetLat <- -5.6
# title <- "FMA Precipitation for Papua New Guinea"
# plotDes <- 'PNG'
# sInd <- 3

# Ethiopia
# targetLon <- 40.7
# targetLat <- 8
# title <- "JAS Precipitation for Ethiopia"
# plotDes <- 'ETH'
# sInd <- 8

# DRC
targetLon <- 23.5
targetLat <- -1
title <- "OND Precipitation for DRC"
plotDes <- 'DRC'
sInd <- 11

# NYC
targetLon <- -74
targetLat <- 41
title <- "MAM Precipitation for NYC"
plotDes <- 'NYC'
sInd <- 4
###############################################################################
# Load in the relevant PPT Data from CRU
###############################################################################

# load in the ncdf metadata
handle <- nc_open(sprintf('%s/Raw/%s',ddir,pptFile))

lon <- ncvar_get(handle,'lon')
lat <- ncvar_get(handle,'lat')

xInd <- which.min(abs(lon-targetLon))
yInd <- which.min(abs(lat-targetLat))

# set up all of the necessary time objects to properly subset the data
time <- ncvar_get(handle,'time')

tFull <- seq(1901,endYear+1,length=length(time)+1)[1:length(time)]
year <- floor(tFull)
month <- rep(1:12,length(time)/12)

fullTimeMat <- cbind(year,month)

# subset based on starting year and make all of the necessary time objects
timeInds <- which(year >= startYear)

timeMat <- fullTimeMat[timeInds,]

# grab the ppt monthly time series for the single location
ppt <- ncvar_get(handle, 'pre', start=c(xInd,yInd,timeInds[1]),
								count=c(1,1,length(timeInds)) )

counts <- ncvar_get(handle, 'stn', start=c(xInd,yInd,timeInds[1]),
								count=c(1,1,length(timeInds)) )
nc_close(handle)

# get the ppt series arranged all pretty
pptMonthly  <- matrix(ppt,ncol=12,byrow=T)
pptSeasonal <- seasonalSum12(ppt,tYear)

countsMonthly <- matrix(counts,ncol=12,byrow=T)
countsSeasonal <- seasonalSum12(counts,tYear)


###############################################################################
# Load in the relevant PPT Data from CHIRPS
###############################################################################
tYear2 <- 1981:2018
# load in the ncdf metadata
handle <- nc_open(sprintf('%s/Raw/%s',ddir,'CHIRPS.mon.0.5x0.5.1981.2018.nc'))

lon2 <- ncvar_get(handle,'X')
lat2 <- ncvar_get(handle,'Y')

# set up all of the necessary time objects to properly subset the data
time2 <- ncvar_get(handle,'T')

nYears2 <- ceiling(length(time2)/12)

tFull2 <- seq(1981,by=1/12,length=length(time2))
year2 <- floor(tFull2)
month2 <- rep(1:12,ceiling(length(time2)/12))[1:length(time2)]

timeMat2 <- cbind(year2,month2)

# figure out which gridbox corresponds
xInd2 <- which(lon2 == lon[xInd])
yInd2 <- which(lat2 == lat[yInd])


# check some data things
testField <- ncvar_get(handle, 'precipitation', start=c(1,1,1), count=c(-1,-1,1))


# grab the ppt monthly time series for the single location
ppt2 <- ncvar_get(handle, 'precipitation', start=c(xInd2,yInd2,1),
								count=c(1,1,-1) )

#clean
ppt2[ppt2<0 | ppt2 > 10000] <- NA

nc_close(handle)

# add NAs as needed
ppt2Full <- c(ppt2,rep(NA,12 - length(time2) %%12))

# get the ppt series arranged all pretty
pptMonthly2  <- matrix(ppt2Full,ncol=12,byrow=T)
pptSeasonal2 <- seasonalSum12(ppt2Full,tYear2)

###############################################################################
# Plot up the two series and think about how to combine
###############################################################################

# get the series and time objs
cruSeries <- pptSeasonal[,sInd]
cruTime <- as.numeric(labels(cruSeries))
chirpsSeries <- pptSeasonal2[,sInd]
chirpsTime <- as.numeric(labels(chirpsSeries))

# build climatologies
cruClimYears <- 1961:1990
chirpsClimYears <- 1982:2011

cruClim    <- mean(cruSeries[which(cruTime %in% cruClimYears)],na.rm=T)
chirpsClim <- mean(chirpsSeries[which(chirpsTime %in% chirpsClimYears)],na.rm=T)
cruTerc    <- quantile(cruSeries[which(cruTime %in% cruClimYears)],probs=c(1/3,2/3))
chirpsTerc <- quantile(chirpsSeries[which(chirpsTime %in% chirpsClimYears)],
	probs=c(1/3,2/3))


xr <- range(cruTime,chirpsTime,na.rm=T)
yr <- range(cruSeries,chirpsSeries,na.rm=T)

# Plot the raw information
pdf(sprintf('%s/raw_%s.pdf',plotdir,plotDes),10,7)
plot(cruTime,cruSeries,type='l',xlim=xr,ylim=yr,main=title,lwd=2,
	xlab='Year', ylab='Precipitation (mm/season)')
points(cruTime,rep(cruClim,length(cruTime)),type='l',col='black',lwd=1.5)
points(cruTime,rep(cruTerc[1],length(cruTime)),type='l',col='black',lwd=1.5,lty=3)
points(cruTime,rep(cruTerc[2],length(cruTime)),type='l',col='black',lwd=1.5,lty=3)
points(chirpsTime,chirpsSeries,type='l',col='blue',lwd=2)
points(chirpsTime,rep(chirpsClim,length(chirpsTime)),type='l',col='blue',lwd=1.5)
points(chirpsTime,rep(chirpsTerc[1],length(chirpsTime)),
	type='l',col='blue',lwd=1.5,lty=3)
points(chirpsTime,rep(chirpsTerc[2],length(chirpsTime)),
	type='l',col='blue',lwd=1.5,lty=3)

par(new=TRUE)
plot(cruTime,countsSeasonal[,sInd],col='red',type='l',
	lty=2, xaxt="n", yaxt="n",ylab="", xlab="",ylim=c(0,35))

legend('topleft', c('CRU TS4.01', 'CHIRPS'), lwd=2, col=c('black', 'blue'))
dev.off()

# normalize the series individually and see how they line up
# using a gaussian standardization which is probably not good (should be gamma probably)

cruSd    <- sd(cruSeries[which(cruTime %in% cruClimYears)],na.rm=T)
chirpsSD <- sd(chirpsSeries[which(chirpsTime %in% chirpsClimYears)],na.rm=T)

cruNormal    <- (cruSeries- cruClim)/cruSd
chirpsNormal <- (chirpsSeries - chirpsClim)/chirpsSD

cruTercNormal    <- (cruTerc - cruClim)/cruSd
chirpsTercNormal <- (chirpsTerc - chirpsClim)/chirpsSD

xr <- range(cruTime,chirpsTime)
ymax <- max(abs(cruNormal),abs(chirpsNormal),na.rm=T)
yr <- c(-ymax, ymax)

pdf(sprintf('%s/standarized_%s.pdf',plotdir,plotDes),10,7)
plot(cruTime,cruNormal,type='l',xlim=xr,ylim=yr,main=title,lwd=2,
	xlab='Year', ylab='Normalized Precipitation (sigma/season)')
points(chirpsTime,chirpsNormal,type='l',col='blue',lwd=2)


points(cruTime,rep(cruTercNormal[1],length(cruTime)),
	type='l',col='black',lwd=1.5,lty=3)
points(cruTime,rep(cruTercNormal[2],length(cruTime)),
	type='l',col='black',lwd=1.5,lty=3)
points(chirpsTime,rep(chirpsTercNormal[1],length(chirpsTime)),
	type='l',col='blue',lwd=1.5,lty=3)
points(chirpsTime,rep(chirpsTercNormal[2],length(chirpsTime)),
	type='l',col='blue',lwd=1.5,lty=3)
abline(h=0)

countsNormal <- (countsSeasonal[,sInd])/sd(countsSeasonal[,sInd])
points(cruTime,countsNormal,col='red',type='l',lty=2)
legend('topleft', c('CRU TS4.01', 'CHIRPS'), lwd=2, col=c('black', 'blue'))
dev.off()

###############################################################################
# Explore error in tercile assessment
###############################################################################


compareExceedences <- function(s1,s2,t1,t2,terc1,terc2){
	overlapPeriod <- intersect(t1,t2)
	
	overlapTable  <- matrix(0,3,3)
	rownames(overlapTable) <- c('L1','M1','H1')
	colnames(overlapTable) <- c('L2','M2','H2')

	for(i in 1:length(overlapPeriod)){
		result1 <- whichTercile(s1[t1==overlapPeriod[i]],terc1)
		result2 <- whichTercile(s2[t2==overlapPeriod[i]],terc2)

		overlapTable[result1,result2] <- overlapTable[result1,result2] + 1
	}

	return(overlapTable)
}

compareExceedences(cruSeries,chirpsSeries,)

overlapPeriod <- intersect(cruTime,chirpsTime)
overlapTable  <- matrix(0,3,3)

rownames(overlapTable) <- c('L1','M1','H1')
colnames(overlapTable) <- c('L2','M2','H2')

for(i in 1:length(overlapPeriod)){
	cruResult <- whichTercile(cruSeries[cruTime==overlapPeriod[i]],cruTerc)
	chirpsResult <- whichTercile(chirpsSeries[chirpsTime==overlapPeriod[i]],chirpsTerc)

	overlapTable[cruResult,chirpsResult] <- overlapTable[cruResult,chirpsResult] + 1
}













