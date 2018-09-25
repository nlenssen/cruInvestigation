source('Code/namelist.Rnl')


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

# grab the ppt monthly time series
ppt <- ncvar_get(handle, 'pre', start=c(1,1,timeInds[1]),
								count=c(-1,-1,length(timeInds)) )

counts <- ncvar_get(handle, 'stn', start=c(1,1,timeInds[1]),
								count=c(-1,-1,length(timeInds)) )
nc_close(handle)

# get the ppt series arranged all pretty
pptSeasonal <- seasonalField12(ppt,tYear,sum)
countsSeasonal <- seasonalField12(counts,tYear,min)

# get the correct time objects
timeMapWorking <- timeMat[12:(nrow(timeMat)-1),]
timeMap <- cbind(timeMapWorking,timeMapWorking[,1] + (timeMapWorking[,2]-1)/12)
save(lon,lat,timeMap,pptSeasonal,countsSeasonal,
	file=sprintf('%s/seasonalCRU.Rda',wddir))