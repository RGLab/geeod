bootstrap.patients.vaccine.permute <- function(data,n,vaccineProp) {
  ids <- unique(data$ptid)
  indexes <- lapply(ids,function(x) x=which(c(data$ptid==x)))

  bootsamp <- sample(length(ids),n,replace=TRUE)
  bootfunction <- function(i,vaccineProp) {
    lines <- data[indexes[[bootsamp[i]]],]
    lines$ptid <- i
    lines$vaccine <- rbinom(1,1,prob=vaccineProp)
    lines
  }

  newdata <- lapply(1:length(bootsamp),bootfunction,vaccineProp)
  newdata <- do.call("rbind",newdata)
}

bootstrap.proportionally <- function(data,n,vaccineProp) {
  ids <- unique(data$ptid)
  indexes <- lapply(ids,function(x) x=which(c(data$ptid==x)))
  vaccineStatus <- sapply(indexes,function(x) data$vaccine[x[[1]]])

  sampVaccine <- sample(which(vaccineStatus=="VACCINE"),round(n*vaccineProp),replace=TRUE)
  sampPlacebo <- sample(which(vaccineStatus=="PLACEBO"),n-round(n*vaccineProp),replace=TRUE)
  bootsamp <- c(sampVaccine,sampPlacebo)
  randomVaccine <- sample(c(rep(1,round(n*vaccineProp)),rep(0,n-round(n*vaccineProp))),n,replace=FALSE)
  bootfunction <- function(i) {
    lines <- data[indexes[[bootsamp[i]]],]
    lines$ptid <- i
    lines$vaccine <- randomVaccine[i]
    lines
  }

  newdata <- lapply(1:length(bootsamp),bootfunction)
  newdata <- do.call("rbind",newdata)
}

leaves <- unique(rv144$population)
leaves <- leaves[!(leaves %in% unique(rv144$parent))]

run.typeIsim <- function(config,leaves,replications) {

  checkInverse <- function(m) class(try(solve(t(m)%*%m),silent=T))=="matrix"

  n <- config[[1]]
  vaccineProp <- config[[2]]

  pvalBinom <- matrix(nrow=replications,ncol=length(leaves))
  pvalBetaBinom <- matrix(nrow=replications,ncol=length(leaves))
  allLeaves <- rv144[-which(rv144$stim=="env"),]
  for(i in 1:length(leaves)) {
    data <- rv144[which(rv144$population %in% leaves[1:4]),]

    for(m in 1:replications) {
      modelmat <- NULL
      while(FALSE==checkInverse(modelmat)) {
        countPlacebo <- 0
        tempdata <- bootstrap.proportionally(data,n,vaccineProp)
        modelmat <- model.matrix(count ~ age + gender + stim*vaccine,data=tempdata)
      }

      fitbinom <- try(geem_betabinomial(count ~ age + gender + stim*vaccine,
                                    id=ptid,N=parentcount,waves=NULL,data=tempdata,
                                    rho.fixed=TRUE,corstr = "exchangeable"))
      fitbetabinom <- try(geem_betabinomial(count ~population/(age + gender + stim*vaccine),
                                        id=ptid,N=parentcount,waves=NULL,data=data,
                                        rho.fixed=FALSE,corstr = "exchangeable"))
      if(class(fitbinom)=="geem") {
        pvalBinom[m,i] <- summary(fitbinom)$p[7]
      } else {
        pvalBinom[m,i] <- NA
      }

      if(class(fitbetabinom)=="geem") {
        pvalBetaBinom[m,i] <- summary(fitbetabinom)$p[7]
      } else {
        pvalBetaBinom[m,i] <- NA
      }
      print(c(m,leaves[[i]],binom=round(mean(pvalBinom[1:m,i]<0.05,na.rm=TRUE),3),bb=round(mean(pvalBetaBinom[1:m,i]<0.05,na.rm=TRUE),3)))
    }
  }

  pvalBinom <- data.frame(pvalBinom)
  names(pvalBinom) <- leaves

  pvalBetaBinom <- data.frame(pvalBetaBinom)
  names(pvalBetaBinom) <- leaves

  config <- t(replicate(replications,config,simplify=TRUE))
  pvalBinom <- data.frame(config,reg="binom",pvalBinom)
  #print(pvalBinom)
  pvalBetaBinom <- data.frame(config,reg="betabinom",pvalBetaBinom)
  #print(pvalBetaBinom)
  result <- rbind(pvalBinom,pvalBetaBinom)
  return(result)
}

configurations <- expand.grid(n=c(32,64,128,256,512),vaccineProp=seq(from=0.5,to=0.85,length.out=3))

simResults <- apply(configurations,1,run.typeIsim,leaves[1:7],replications=100)
#save(simResults,file="type I error experiment 120715.Robj")
require(ggplot2)
require(reshape2)
results <- do.call("rbind",simResults)
results <- melt(results,id.vars=c("n","vaccineProp","reg"))
names(results)[4:5] <- c("cytokine","pval")
ggplot(results,aes(color=as.factor(n), x=pval,linetype=reg)) + stat_ecdf() +
  facet_grid(.~vaccineProp,labeller="label_both") + theme_bw() +
  xlim(0,1) + ylim(0,1) + geom_abline(intercept=0,slope=1)

