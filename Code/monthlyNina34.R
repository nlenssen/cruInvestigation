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

# remove the climatology for each month over 1981-2010
climStartInd <- which(tYear==1981)
climatologyLength <- 30

climatology <- apply(nina34Raw[climStartInd:(climStartInd+climatologyLength-1),-1],2,mean)

nina34Anomalies <- nina34Raw
for(i in 2:13){
	nina34Anomalies[,i] <- nina34Raw[,i] - climatology[i-1]
}

# rearrange to vector format
nina34 <- c(t(nina34Anomalies[3:(nrow(nina34Anomalies)),2:13]))

# smooth with 5-month moving average following Trenberth 
smoothNina <- filter(nina34,rep(1,5)/5)

# fix date stuff now that we lose the end points
naInds <- which(is.na(smoothNina))
continuousTime <- continuousTimeFull[-naInds]
timeMap <- timeMapFull[-naInds,]

# center, standardize, and detrend
smoothNinaFinal <- smoothNina[-naInds]
stdNina <- (smoothNinaFinal - mean(smoothNinaFinal))/sd(smoothNinaFinal)
dtdNina <- stdNina - predict(lm(stdNina~continuousTime))

nina34Final <- dtdNina

ninaMonthly <- list(nina34 = nina34Final, timeMap = timeMap)

save(ninaMonthly, file=sprintf('%s/ninaMonthly.Rda',wddir))

# plot(ninaMontly$timeMap[,3],ninaMontly$nina34,type='l')
# abline(h=0)
# abline(h=1,col='red')
# abline(h=-1,col='blue')



