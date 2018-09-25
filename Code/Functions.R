seasonalAverage12 <- function(ts,tYear){

	nMonth <- length(ts)
	nYear <- length(tYear)

	outMat <- matrix(NA, nrow=length(tYear)-1, ncol=12)
	
	rownames(outMat) <- tYear[2:nYear]
	colnames(outMat) <- c('DJF', 'JFM', 'FMA', 
						  'MAM', 'AMJ', 'MJJ',
						  'JJA', 'JAS', 'ASO',
						  'SON', 'OND', 'NDJ')

	for(i in 2:nYear){
		tsInd <- 1 + 12*(i-1)

		for(j in 1:12){
			seasonInds <- ( (j-2) : j)
			outMat[i-1,j] <- mean(ts[tsInd+seasonInds])
		}
	}

	return(outMat)
}


seasonalSum12 <- function(ts,tYear){

	nMonth <- length(ts)
	nYear <- length(tYear)

	outMat <- matrix(NA, nrow=length(tYear)-1, ncol=12)
	
	rownames(outMat) <- tYear[2:nYear]
	colnames(outMat) <- c('DJF', 'JFM', 'FMA', 
						  'MAM', 'AMJ', 'MJJ',
						  'JJA', 'JAS', 'ASO',
						  'SON', 'OND', 'NDJ')

	for(i in 2:nYear){
		tsInd <- 1 + 12*(i-1)

		for(j in 1:12){
			seasonInds <- ( (j-2) : j)
			outMat[i-1,j] <- sum(ts[tsInd+seasonInds])
		}
	}

	return(outMat)
}

seasonalMin12 <- function(ts,tYear){

	nMonth <- length(ts)
	nYear <- length(tYear)

	outMat <- matrix(NA, nrow=length(tYear)-1, ncol=12)
	
	rownames(outMat) <- tYear[2:nYear]
	colnames(outMat) <- c('DJF', 'JFM', 'FMA', 
						  'MAM', 'AMJ', 'MJJ',
						  'JJA', 'JAS', 'ASO',
						  'SON', 'OND', 'NDJ')

	for(i in 2:nYear){
		tsInd <- 1 + 12*(i-1)

		for(j in 1:12){
			seasonInds <- ( (j-2) : j)
			outMat[i-1,j] <- min(ts[tsInd+seasonInds])
		}
	}

	return(outMat)
}

seasonalField12 <- function(arr,tYear,FUN){
	nlon <- dim(arr)[1]
	nlat <- dim(arr)[2]

	outTime <- (length(tYear)-1) * 12 # not sure about this dimensionality rn

	outArray <- array(NA, dim = c(nlon,nlat,outTime))

	pb   <- txtProgressBar(2, nYear, style=3)

	for(i in 2:nYear){
		setTxtProgressBar(pb, i)

		tsInd <- 1 + 12*(i-1)

		for(j in 1:12){
			seasonInds <- ( (j-2) : j)
			if(max(tsInd+seasonInds) <= dim(arr)[3]){
				outArray[,,((i-2)*12 + j)] <- apply(arr[,,tsInd+seasonInds],c(1,2),FUN)
			}
		}
	}

	return(outArray)
}


redBlue <- function(n=256){
	require(fields)
	require(RColorBrewer)
	rev(designer.colors(n, brewer.pal(11,'RdBu')))
}