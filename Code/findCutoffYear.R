# the goal of this script is to find the year in which each gridbox
# becomes 'data poor' in the CRU dataset. 

source('Code/namelist.Rnl')

# load in all the CRU data for now

load(sprintf('%s/seasonalCRU.Rda',wddir))

# Find the first occurance of a station less than 2 for each gridbox

cutoffInd <- function(x){
	min(which(x<2))
}

cutoffField <- apply(countsSeasonal,c(1,2),cutoffInd)

# handle the locations that stay data rich and are ocean 
cutoffField[cutoffField == Inf] <- nrow(timeMap)
cutoffField[is.na(pptSeasonal[,,1])]  <- NA


# restrict to \pm 50 degrees lat to be on chirps grid
latInds <- which(lat <= 50 & lat >= -50)
lat2 <- lat[latInds]

# convert to time
cutoffDate <- apply(cutoffField,c(1,2), function(x) timeMap[x,3])

pal <- designer.colors(256,brewer.pal(9,"OrRd"))

pdf(sprintf('%s/yearCutoff.pdf',plotdir),10,7)
image.plot(lon,lat2,cutoffDate[,latInds],col=pal,
	xlab='',ylab='',main='CRU TS 4.01 Precipitation Data Poor Year (1980 contour)')
contour(lon,lat2,cutoffDate[,latInds],levels=c(1990,1980),add=T,col=c('green','blue'))
world(add=T,col='grey40')
dev.off()

# make the mask to ignore pre-1960 for now
dataPoorMask <- ifelse(cutoffDate < 1965,NA,1)

# save relevant files for further work
save(cutoffField,cutoffDate,dataPoorMask,
	file=sprintf('%s/cutoffInfo.Rda',wddir))