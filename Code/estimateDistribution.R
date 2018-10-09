source('Code/namelist.Rnl')

# load in data, 
load(sprintf('%s/seasonalCRU.Rda',wddir))
load(sprintf('%s/cutoffInfo.Rda',wddir))

## Operate on each season seperately. Start with just DJF (season 1)
# Restrict to 50/50, pre-1990, and apply data poor mask
maxYear <- 1991

latInds <- which(lat <= 50 & lat >= -50)
lat2 <- lat[latInds]

sInd <- 1

timeInds <- which(timeMap[,2]==sInd & timeMap[,1] < maxYear)
subFieldTemp <- pptSeasonal[,latInds,timeInds]

# element-multiply each array slice by the mask
subField <- subFieldTemp * c(dataPoorMask[,latInds])


## get climatology (1961-1990) and convert to mult-anomalies
climStartInd <- which(timeMap[timeInds,1]==1961)
climInds <- climStartInd:(climStartInd+29)

climField <- apply(subField[,,climInds],c(1,2),mean)

# convert to mult-amonalies
multAnom    <- array(NA, dim=dim(subField))
logMultAnom <- array(NA, dim=dim(subField))

for(i in 1:dim(subField)[3]){
	multAnom[,,i] <- subField[,,i]/climField
	
	tempLog <- log(multAnom[,,i])
	tempLog[tempLog==-Inf] <- NaN
	logMultAnom[,,i] <- tempLog
}

## Get non-parametric terciles from entire period
npTercile <- aperm(apply(subField,c(1,2),quantile,
	 probs=c(1/3,2/3),na.rm=T),c(2,3,1))
npTercileAnom <- aperm(apply(multAnom,c(1,2),quantile,
	 probs=c(1/3,2/3),na.rm=T),c(2,3,1))

#######
# TODO: Visualize!!!
#######



testFit   <- array(NA, dim=c(nrow(subField),ncol(subField),2))
testFitSD <- array(NA, dim=c(nrow(subField),ncol(subField),2))
fitObjs <- list()

pb   <- txtProgressBar(1, nrow(subField), style=3)
for(i in 1:nrow(subField)){
	setTxtProgressBar(pb, i)
	for(j in 1:ncol(subField)){
		tempFit <- gammaFit(subField[i,j,])
		
		if(!is.na(tempFit)){
			testFit[i,j,] <- tempFit$estimate
			testFitSD[i,j,] <- tempFit$sd
		}
		fitObjs[[nrow(subField) * (i-1) + j]] <- tempFit
	}
}


## check rough spatial structure

set.panel(2,1)
image.plot(lon,lat2,testFit[,,1],
	main=expression(paste('Shape (', k, ') Parameter from Gamma MLE')))
world(add=T,col='black')

image.plot(lon,lat2,testFitSD[,,1],main='Standard Error')
world(add=T,col='black')



## check SE on gamma fit for the data
set.panel(2,1)
image.plot(lon,lat2,testFit[,,2],
	main=expression(paste('Scale (',theta, ') Parameter from Gamma MLE')))
world(add=T,col='black')

image.plot(lon,lat2,testFitSD[,,2],main='Standard Error')
world(add=T,col='black')

## Think about other goodness of fit measures? compare with gaussian?


