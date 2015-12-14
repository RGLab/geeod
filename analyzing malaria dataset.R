leaves <- unique(malaria$population)
leaves <- leaves[which(!(leaves %in% unique(malaria$parent)))]

#subsetting the data quite arbitrarily
leaves <- leaves[1:7]
data <- malaria[which(malaria$population %in% leaves),]

fittedList <- lapply(1:length(leaves),c)
for(i in 1:length(leaves)) {
  tempdata <- data[which(data$population==leaves[i]),]
  tempdata <- tempdata[-which(tempdata$visitno=="pos"),]
  fit <- geem_betabinomial(count~stim*visitno,N=parentcount,id=paste(tempdata$ptid,tempdata$visitno,sep=""),data=tempdata,waves=NULL,
                           corstr="exchangeable")
  fittedList[[i]] <- fit
  print(fit$rhos)
  coef <- coef(fit)[9:20]
  var <- fit$var[9:20,9:20]
  wald <- as.numeric(t(coef) %*% solve(var) %*% coef)
  pval <- 1-pchisq(wald,length(coef))
  print(pval)

  yhat <- predictBB(fit,tempdata$parentcount)
  residual <- tempdata$count-yhat
  residual <- matrix(residual,ncol=20)
}

