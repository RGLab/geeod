require(devtools)
statsData <- readRDS("data-raw/stats.rds")
pdata <- readRDS("data-raw/pdata.rds")

names(statsData) <- tolower(names(statsData))
names(pdata) <- tolower(names(pdata))
pdata$stim <- sapply(pdata$stim,as.character)

names(pdata)[1:2] <- c("experiment","order")
pdata$visitno[which(grepl("Pos",pdata$visitno))] <- "pos"
pdata$stim[which(pdata$stim %in% c("ctr SPZ","Cells Only"))] <- "control" #is this correct?

malaria <- merge(statsData,pdata,by.x="name",by.y="name")
malaria <- malaria[-which(malaria$order==99),] #removing the 99 sample order because I don't know what they mean

use_data(malaria,overwrite=TRUE)





