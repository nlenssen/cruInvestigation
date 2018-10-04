source("Code/namelist.Rnl")

###############################################################################
# Load in the ENSO time series and identify enso years
###############################################################################
load(sprintf('%s/ninaMonthly.Rda',wddir))
load(sprintf('%s/ninaSeasonal.Rda',wddir))

###############################################################################
# Load in the relevant PPT Data
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




# just work with one season right now
sInd <- 12

# subset all the data to just the sInd seasonal time series
tempPptMonthly  <- pptMonthly[-1,sInd]
tempPptSeasonal <- pptSeasonal[,sInd]

tempCounts      <- countsMonthly[-1,sInd]
tempCountsSeasonal<- countsSeasonal[-1,sInd]

tYearFinal <- tYear[-1]

###############################################################################
# build the climatology for the gridbox at a monthly and seasonal time scale
###############################################################################
climMat <- matrix(c(1961,1971,1981,1990,2000,2010),ncol=2)

exceedMat <- matrix(NA, nrow=nrow(climMat),ncol=2)

for(i in 1:nrow(climMat)){
	pdf(sprintf('%s/clim%4.0f_%4.0f.pdf',plotdir,climMat[i,1],climMat[i,2]),10,7)
	climInds <- which(tYearFinal == climMat[i,1]):which(tYearFinal == climMat[i,2])
	seasonalClimatology <- mean(tempPptSeasonal[climInds])

	seasonalTerciles <- quantile(tempPptSeasonal[climInds],probs=c(1/3,2/3))

	# exploratory plot
	plot(tYearFinal,tempPptSeasonal,type='l',col='grey40',lwd=1.5,
		main=sprintf('NDJ precipitation for (%.2f,%.2f) with %4.0f-%4.0f Climatology',
			lon[xInd],lat[yInd],climMat[i,1],climMat[i,2]),
		xlab='Year', ylab='Total Seasonal Precipitation (mm)')
	abline(h=seasonalClimatology,col='red',lwd=1.5)
	abline(h=seasonalTerciles,lty=3,col='blue',lwd=1.5)

	points(tYearFinal,tempCounts*(50),col='black',type='l',lty=2)
	dev.off()
}

# think about how to calculate how this could effect exceedences?

# well, they are using the 1961-1990 as the climatology


###############################################################################
# view how this effects the variance of the series
###############################################################################

slideSize <- c(5,6,8,10,20)

for(i in 1:length(slideSize)){
	pdf(sprintf('%s/sdyear%1.0f.pdf',plotdir,slideSize[i]),10,7)
	maxInd <- length(tempPptSeasonal) - slideSize[i] + 1

	slidingSd <- rep(NA, maxInd)
	meanCount <- rep(NA, maxInd)
	# 

	for(j in 1:maxInd){
		slidingSd[j] <- sd(tempPptSeasonal[j:(j+slideSize[i])])
		meanCount[j]    <- mean(tempCountsSeasonal[j:(j+slideSize[i])])
	}


	plot(tYearFinal[1:maxInd]+slideSize[i]/2,slidingSd,type='l',
				xlim=range(tYearFinal),ylim=c(0,max(slidingSd,na.rm=T)),
				main=sprintf('NDJ %1.0f-year SD for (%.2f,%.2f)',
							slideSize[i],lon[xInd],lat[yInd]),
				xlab='Midpoint Year', ylab='Standard Deviation (mm)')

	points(tYearFinal[1:maxInd]+slideSize[i]/2,meanCount*(8),
		col='black',type='l',lty=2)
	dev.off()
}






















