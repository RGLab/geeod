library(devtools)
rvData <- readRDS("data-raw/rv144_cc.rds")
names(rvData) <- tolower(names(rvData))

#more convenient stimulous names\
newStimNames <- c("sebctrl","negctrl","env")
rvData$stim <- as.character(rvData$stim)
for(i in 1:length(unique(rvData$stim))) {
  rvData$stim[which(rvData$stim==unique(rvData$stim)[i])] <- newStimNames[[i]]
}

rvData <- rvData[which(rvData$population %in% leaves),]
rvData <- rvData[-which(rvData$parentcount==0),]
rv144 <- rvData
rm(rvData)

use_data(rv144,overwrite=TRUE)


