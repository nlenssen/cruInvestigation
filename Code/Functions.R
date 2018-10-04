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


# Function to compare two series empirical tercile given a climatology
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

whichTercile <- function(val,terc){
	if(val<terc[1]){
		return(1)
	} else if(val < terc[2]){
		return(2)
	} else{
		return(3)
	}
}



redBlue <- function(n=256){
	require(fields)
	require(RColorBrewer)
	rev(designer.colors(n, brewer.pal(11,'RdBu')))
}

symPlot <- function(x,y,z,zmax=NULL,pal=redBlue(),...){
	if(is.null(zmax)){
		zmax <- max(abs(z),na.rm=T)
		zr   <- c(-zmax,zmax)
	}

	image.plot(x,y,z,zlim=zr,col=pal,...)
	world(add=T)
}


