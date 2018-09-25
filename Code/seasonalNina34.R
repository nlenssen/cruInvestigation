source('Code/namelist.Rnl')

nYear <- (endYear)-startYear+2
nMonth <- nYear * 12

tYear <- startYear:endYear

continuousTimeFull <- seq(startYear,endYear+2,length=nMonth+1)[1:nMonth]
tYearFull  <- rep(startYear:(endYear+2),each=12)[1:nMonth]
tMonthFull <- rep(1:12, nYear)[1:nMonth]

timeMapFull <- cbind(tYearFull,tMonthFull,continuousTimeFull)


# read in the raw data
nina34Raw <- as.matrix(read.table(sprintf('%s/Raw/%s',ddir,ensoFile),header=FALSE,skip=1))


# (0) create 3.4 anomalies for each month using 1981-2010 climatology
climStartInd <- which(tYear==1981)
climatologyLength <- 30

climatology <- apply(nina34Raw[climStartInd:(climStartInd+climatologyLength-1),-1],2,mean)

nina34Anomalies <- nina34Raw
for(i in 2:13){
	nina34Anomalies[,i] <- nina34Raw[,i] - climatology[i-1]
}

nina34 <- c(t(nina34Anomalies[3:(nrow(nina34Anomalies)),2:13]))

# (1) Smooth the montly 3.4 index with 3-pt binomial filter

coeffs <- c(1,4,6,4,1)
smoothNina <- filter(nina34,coeffs/sum(coeffs))

# (2) create seasonal mean time series
seasonalNina <- seasonalAverage12(smoothNina,tYear)

# (3) linearly detrend the time series
seasonalNinaVec <- c(t(seasonalNina))
tempTime <- 1:length(seasonalNinaVec)
dtdNina <- matrix(
	seasonalNinaVec - predict(lm(seasonalNinaVec ~ tempTime)),
	nrow(seasonalNina),ncol(seasonalNina),byrow=T)

# standardize each month?
meanSeasonal <- apply(dtdNina,2,mean)
sdSeasonal   <- apply(dtdNina,2,sd)

stdNina <- matrix(NA, nrow=nrow(seasonalNina),ncol=ncol(seasonalNina))

for(i in 1:ncol(stdNina)){
	stdNina[,i] <- (dtdNina[,i] - meanSeasonal[i])/sdSeasonal[i]
}

# check that we have approximate agreement with Deser in the NDJ years
ndjSeries <- cbind(tYear[-1],stdNina[,12])

# Deser gets 14 El Ninos vs our 10 over the overlapping time period
#	Missing 1963 (0.974), 1968 (0.836), 1987 (0.961), 2006 (0.736)
ninoYears <- ndjSeries[which(ndjSeries[,2] > 1),]

# We get the same amount of La Nina years at 9!
ninaYears <- ndjSeries[which(ndjSeries[,2] < -1),]


# save the series to work with for now. Chat with Lisa and Simon about the
# best way to process the two series
save(stdNina,file=sprintf('%s/ninaSeasonal.Rda',wddir))


# Following Deser et al 2017 (Section 2.c.1 p5)
# ENSO Compositing:

# We identify El Niño (EN) and La Niña (LN) events
# using a standard approach based on the Niño-3.4 SST
# index [SST anomalies averaged over 1208–1708W,
# 58S–58N (Barnston et al. 1997)]. Here we use ERSSTv3b
# to construct the Niño-3.4 index.

# (1) We first smooth the monthly Niño-3.4 SST index with 
#     a three-point binomial filter
# (2) average it over November–January (NDJ)
# (3) linearly detrend the time series
# (4) and identify EN (LN) events when the detrended index 
#     exceeds (falls below minus) one standard deviation 
#
# This procedure yields 18 EN and 14 LN events during 1920–
# 2013; these events are listed in Tables 1 and 2. Note that
# the Niño-3.4 index is based on NDJ, one month ahead
# of the DJF season used for the extratropical SLP composites
# to accommodate the time scale of the Rossby
# wave response to tropical heating [e.g., Alexander et al.
# (2002) and references therein]. We form ENSO composites
# by subtracting the average of the 14 LN events
# from the average of the 18 EN events. Similar ENSO
# composites are obtained using the Niño-4 SST index (SST
# anomalies averaged over 1608E–1508W, 58S–58N) in
# place of Niño-3.4 (not shown).