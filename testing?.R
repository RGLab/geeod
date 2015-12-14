leaves <- unique(rv144$population)
leaves <- leaves[!(leaves %in% unique(rv144$parent))]
#  Checking if workds for all leaves
# for(i in 1:length(leaves)) {
#   print(i)
#   data <- rv144[which(rv144$population==leaves[[i]]),]
  result <- geem_betabinomial(count ~ age + gender + stim*vaccine,
                              id=data$ptid,N=parentcount,waves=NULL,data=data,
                              rho.fixed=FALSE,corstr = "exchangeable")
#   print(result$rhos)
#   fit2 <- geeglm(cbind(count,parentcount-count) ~ age + gender + stim*vaccine,id=ptid,
#          corstr = "exchangeable",data=data,family="binomial")
# }

# Checking what happens with more than one leaf
data <- rv144[which(rv144$population %in% leaves[1:16]),]
result <- geem_betabinomial(count ~ population*(age + gender + stim*vaccine),
                            id=data$ptid,N=parentcount,waves=NULL,data=data,
                            rho.fixed=FALSE,corstr = "exchangeable")

data <- rv144[which(rv144$population==leaves[[1]]),]
data$ptid <- as.character(data$ptid)
nothreeobs <- which(table(data$ptid)!=3)
if(!is.null(nothreeobs)) data <- data[-which(data$ptid==names(nothreeobs)),]
data$ptid <- as.factor(data$ptid)

N <- data$parentcount
result <- geem_betabinomial(count ~ age + gender + stim*vaccine,
                  id=data$ptid,N=parentcount,waves=NULL,data=data)
rho <- result$rhos[length(result$rhos)]

#  checking how well estimates rho ------------------------
generate.temp.data <- function(data,n,rho) {
  require(VGAM)
  nsubject <- length(unique(data$ptid))
  boot_patient <- rmultinom(1,n,prob=rep(1,nsubject))

  construct.dataset <- function(index) {
    data[rep(which(data$ptid==unique(data$ptid)[index]),boot_patient[index]),]
  }
  newdata <- lapply(1:length(boot_patient),construct.dataset)
  newdata <- do.call("rbind",newdata)
  newdata$ptid <- as.vector(sapply(1:n,function(x) rep(x,3)))

  mu <- newdata$count
  #newdata$count <- rbetabinom(nrow(newdata),size=newdata$parentcount,prob=newdata$count/newdata$parentcount,rho=rho)
  #newdata$count[is.na(newdata$count)] <- 0
  return(list(data=newdata,mu=newdata$count))
}

M <- 50
estRhos <- numeric(M)
for(i in 1:M) {
  tempdata <- generate.temp.data(data,400,rho)
  mu <- tempdata$mu
  tempdata <- tempdata$data

  fitbetabinom <- geem_betabinomial(count ~ age + gender + stim*vaccine,
                                    id=data$ptid,N=parentcount,waves=NULL,data=tempdata,
                                    rho.fixed=FALSE,corstr = "exchangeable")


  estRhos[i] <- fitbetabinom$rhos[length(fitbetabinom$rhos)]
  print(estRhos[i])
}

mean(estRhos-rho)
sqrt(mean((estRhos-rho)^2))
#Looks good!

# Testing for false positive rate ---------------------------

data <- rv144[which(rv144$population==leaves[[3]]),]
data$ptid <- as.character(data$ptid)
nothreeobs <- which(table(data$ptid)!=3)
if(!is.null(nothreeobs)) data <- data[-which(data$ptid==names(nothreeobs)),]
data$ptid <- as.factor(data$ptid)

# stim has an actual effect, more than neg and less than seb...
M <- 100
pvalBinom <- numeric(M)
pvalBetaBinom <- numeric(M)
for(i in 1:M) {
  tempdata <- generate.temp.data(data,1000,rho)
  mu <- tempdata$mu
  tempdata <- tempdata$data

  vaccineProp <- mean(tempdata$vaccine=="VACCINE")
  IDs <- unique(tempdata$ptid)
  tempdata$vaccine <- as.numeric(tempdata$vaccine)
  for(j in 1:length(IDs)) {
    subject_rows <- which(tempdata$ptid==IDs[j])
    tempdata$vaccine[subject_rows] <- rbinom(1,1,0.5)
  }

  fitbinom <- geem_betabinomial(count ~ age + gender + stim*vaccine,
                                    id=ptid,N=parentcount,waves=NULL,data=tempdata,
                                    rho.fixed=TRUE,corstr = "exchangeable")
  fitbetabinom <- geem_betabinomial(count ~ age + gender + stim*vaccine,
                                    id=ptid,N=parentcount,waves=NULL,data=tempdata,
                                    rho.fixed=FALSE,corstr = "exchangeable")

  var <- fitbinom$var[8:9,8:9]
  coef <- coef(fitbinom)[8:9]
  pvalBinom[i] <- 1-pchisq(as.numeric(t(coef) %*% solve(var) %*% coef),2)
  var <- fitbetabinom$var[8:9,8:9]
  coef <- coef(fitbetabinom)[8:9]
  pvalBetaBinom[i] <- 1-pchisq(as.numeric(t(coef) %*% solve(var) %*% coef),2)

  print(c(pvalBinom[i],pvalBetaBinom[i]))
}

mean(pvalBinom<0.05)
mean(pvalBetaBinom<0.05)

mean(pvalBinom<0.1)
mean(pvalBetaBinom<0.1)


