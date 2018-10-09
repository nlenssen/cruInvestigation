source('Code/namelist.Rnl')

# load in the processed enso series
load(sprintf('%s/ninaSeasonal.Rda',wddir))

# we do the calculation for each season (column)
ninaIndicator <- matrix(0,nrow=nrow(stdNina),ncol=ncol(stdNina))


for(i in 1:ncol(stdNina)){
	# set the top M nina years for each season to -1
	ninaIndicator[which(stdNina[,i] %in% sort(stdNina[,i],decreasing=FALSE)[1:M]),i] <- -1

	# set the top M nino years for each season to 1
	ninaIndicator[which(stdNina[,i] %in% sort(stdNina[,i],decreasing=TRUE)[1:M]),i] <- 1
}


save(stdNina,ninaIndicator,file=sprintf('%s/ninaSeasonalInd.Rda',wddir))