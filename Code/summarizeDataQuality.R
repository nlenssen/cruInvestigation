source("Code/namelist.Rnl")

###############################################################################
# Load in the relevant count Data
###############################################################################

# load in the ncdf metadata
handle <- nc_open(sprintf('%s/Raw/%s',ddir,pptFile))

lon <- ncvar_get(handle,'lon')
lat <- ncvar_get(handle,'lat')

# set up all of the necessary time objects to properly subset the data
time <- ncvar_get(handle,'time')

tFull <- seq(1901,endYear+1,length=length(time)+1)[1:length(time)]
year <- floor(tFull)
month <- rep(1:12,length(time)/12)

fullTimeMat <- cbind(year,month)

# subset based on starting year and make all of the necessary time objects
timeInds <- which(year >= startYear)

timeMat <- fullTimeMat[timeInds,]
continuousTime <- tFull[timeInds]
# grab the ppt monthly time series for the single location

counts <- ncvar_get(handle, 'stn', start=c(1,1,timeInds[1]),
								count=c(-1,-1,length(timeInds)) )

nc_close(handle)


# Fill array with NAs

counts[counts==-999] <- NA

###################
# mean/median plots
###################

totalMean   <- apply(counts,c(1,2),mean)
totalMedian <- apply(counts,c(1,2),median)
skew        <- totalMedian - totalMean

pdf(sprintf('%s/totalMean.pdf',plotdir),10,7)
image.plot(lon,lat,totalMean,ylim=c(-60,90),main="Mean Station Count")
world(add=T)
dev.off()

pdf(sprintf('%s/totalMedian.pdf',plotdir),10,7)
image.plot(lon,lat,totalMedian,ylim=c(-60,90),main="Median Station Count")
world(add=T)
dev.off()

pdf(sprintf('%s/totalSkewness.pdf',plotdir),10,7)
zmax <- max(skew,na.rm=T)
zr   <- c(-zmax,zmax)
image.plot(lon,lat,skew,ylim=c(-60,90),zlim=zr,main="Median - Mean",col=redBlue())
world(add=T)
dev.off()

##############################################
# Linear trend (over entire time period) plots
##############################################

regSlope <- function(ts,time = continuousTime){
	if(!all(is.na(ts))){
		fit <- lm(ts ~ time)
		return(coef(fit)[2])
	} else{
		return(NA)
	}
}

fullSlope <- apply(counts,c(1,2),regSlope)

# draw the plot
pdf(sprintf('%s/fullSlope.pdf',plotdir),10,7)
zmax <- max(abs(fullSlope),na.rm=T)
zr   <- c(-zmax,zmax)
image.plot(lon,lat,fullSlope,ylim=c(-60,90),zlim=zr,main="Linear Trend of Station Count (1950-2016)",col=redBlue())
world(add=T)
dev.off()

##############################################
# Make a plot for each decade
##############################################
decadeStart <- seq(1950,2010,by=10)

decadalMean <- decadalMedian <- decadalSlope <- 
	array(NA, dim=c(length(lon),length(lat),length(decadeStart)))

for(i in 1:(length(decadeStart))){
	timeStart <- which(continuousTime==decadeStart[i])
	if(decadeStart[i] < 2010){
		timeEnd <- which(continuousTime==decadeStart[i+1])-1
	} else{
		timeEnd <- length(continuousTime)
	}

	subInds <- timeStart:timeEnd

	# subset the data and compute the necessary functions
	subDat <- counts[,,subInds]

	decadalMean[,,i]   <- apply(subDat,c(1,2),mean)
	decadalMedian[,,i] <- apply(subDat,c(1,2),median)
	decadalSlope[,,i]  <- apply(subDat,c(1,2),regSlope,time=continuousTime[subInds])
}


# make all the plots

for(i in 1:length(decadeStart)){
	pdf(sprintf('%s/decadal/totalMean%1.0f.pdf',plotdir,i),10,7)
	image.plot(lon,lat,decadalMean[,,i],ylim=c(-60,90),
		main=sprintf("Mean Station Count (%4.0fs)",decadeStart[i]))
	world(add=T)
	dev.off()

	pdf(sprintf('%s/decadal/totalMedian%1.0f.pdf',plotdir,i),10,7)
	image.plot(lon,lat,decadalMedian[,,i] ,ylim=c(-60,90),
		main=sprintf("Median Station Count (%4.0fs)",decadeStart[i]))
	world(add=T)
	dev.off()

	pdf(sprintf('%s/decadal/totalSkewness%1.0f.pdf',plotdir,i),10,7)
	skew <- decadalMedian[,,i] - decadalMean[,,i]
	zmax <- max(skew,na.rm=T)
	zr   <- c(-zmax,zmax)
	image.plot(lon,lat,skew,ylim=c(-60,90),zlim=zr,
		main=sprintf("Median - Mean (%4.0fs)",decadeStart[i]),col=redBlue())
	world(add=T)
	dev.off()

	pdf(sprintf('%s/decadal/fullSlope%1.0f.pdf',plotdir,i),10,7)
	zmax <- max(abs(decadalSlope[,,i]),na.rm=T)
	zr   <- c(-zmax,zmax)
	image.plot(lon,lat,decadalSlope[,,i],ylim=c(-60,90),zlim=zr,
		main=sprintf("Linear Trend of Station Count (%4.0fs)",decadeStart[i]),col=redBlue())
	world(add=T)
	dev.off()
}

# Make the difference plot between the 1970s and the 1990s 2000s and 2010s

meanDiff90 <- decadalMean[,,5] - decadalMean[,,3]
meanDiff00 <- decadalMean[,,6] - decadalMean[,,3]
meanDiff10 <- decadalMean[,,7] - decadalMean[,,3]

zmax <- max(abs(meanDiff90),abs(meanDiff00),abs(meanDiff10),na.rm=T)
zr   <- c(-zmax,zmax)

pdf(sprintf('%s/decadalDifference1990.pdf',plotdir,i),10,7)
image.plot(lon,lat,meanDiff90,ylim=c(-60,90),zlim=zr,
	main="Change in Mean Coverage from 1970s to 1990s",col=redBlue())
world(add=T)
dev.off()

pdf(sprintf('%s/decadalDifference2000.pdf',plotdir,i),10,7)
image.plot(lon,lat,meanDiff00,ylim=c(-60,90),zlim=zr,
	main="Change in Mean Coverage from 1970s to 2000s",col=redBlue())
world(add=T)
dev.off()

pdf(sprintf('%s/decadalDifference2010.pdf',plotdir,i),10,7)
image.plot(lon,lat,meanDiff10,ylim=c(-60,90),zlim=zr,
	main="Change in Mean Coverage from 1970s to 2010s",col=redBlue())
world(add=T)
dev.off()





