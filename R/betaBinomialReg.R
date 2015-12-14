#' Prediction for Beta-Binomial regression.
#'
#' @param fit A fitted geem_betabinomial object.
#' @param M A vector of total cell counts.
#' @param X An optional matrix of explanatory variables.
#'
#' @return Returns a vector of expected counts
#'
predictBB <- function(fit,N,X=NULL) {
  eta <- predict(fit)
  expectation <- N*exp(eta)/(1+exp(eta))
  return(expectation)
}

#' Fit a beta-binomial generalized estimating equation model.
#'
#' Calculate coefficients and nuisance parameters for the beta-binomial
#' model using generalized estimating equations. Implementation is based
#' on the geeM package.
#'
#' @param formula See corresponding documentation to \code{geem}.
#' @param N A vector of total counts of length \code{reponse}.
#' @param id See corresponding documentation to \code{geem}.
#' @param data See corresponding documentation to \code{geem}.
#' @param corstr See corresponding documentation to \code{geem}.
#' @param Mv See corresponding documentation to \code{geem}.
#' @param weights See corresponding documentation to \code{geem}.
#' @param corr.mat See corresponding documentation to \code{geem}.
#' @param init.beta See corresponding documentation to \code{geem}.
#' @param init.alpha See corresponding documentation to \code{geem}.
#' @param init.phi See corresponding documentation to \code{geem}.
#' @param scale.fix See corresponding documentation to \code{geem}.
#' @param nodummy See corresponding documentation to \code{geem}.
#' @param sandwich See corresponding documentation to \code{geem}.
#' @param maxit See corresponding documentation to \code{geem}.
#' @param tol See corresponding documentation to \code{geem}.
#'
#'
#' @return An object of class "geem_betabinomial
geem_betabinomial <- function(formula,N,id,waves,data=parent.frame(),
                              corstr="independence",weights=NULL,corr.mat=NULL,
                              init.beta=NULL,init.alpha=NULL,init.phi=1,init.rho=0,
                              rho.fixed=FALSE,
                              scale.fix=FALSE,nodummy=FALSE,sandwich=TRUE,
                              maxit=20,tol=1e-05) {

  #Family functions
  LinkFun <- function(mu) {
    p <- mu/N
    log(p/(1-p))
  }

  InvLink <- function(eta) {
    N/(1+exp(-eta))
  }

  InvLinkDeriv <- function(eta) {
    N*exp(-eta)/(1+exp(-eta))^2
  }

  VarFun <- function(mu) {
    p <- mu/N
    N*p*(1-p)*(1+(N-1)*rho)
  }

  #Preparing call
  call <- as.list(match.call())
  call <- call[-1]
  rho <- init.rho
  rhos <- rho
  FunList <- list(LinkFun=LinkFun,VarFun=VarFun,InvLink=InvLink,InvLinkDeriv=InvLinkDeriv)

  #getting N
  if(typeof(data) == "environment"){
    N <- N
  } else {
    if(length(call$N) == 1){
      N.col <- which(colnames(data) == call$N)
      if(length(N.col) > 0) {
        N <- data[,N.col]
      } else {
        N <- eval(call$N, envir=parent.frame())
      }
    }else if(is.null(call$N)){
      stop("N must be specified!")
    }
  }

  varArgs <- list(N=N,rho=rho)
  call$family <- FunList

  #initializing beta with geeglm
  require(geepack)
  env <- parent.frame()
  modelMatCall <- list(object=call$formula,data=call$data)
  modelMat <- do.call("model.matrix",modelMatCall,envir=env)
  solve(t(modelMat)%*%modelMat)

  modelFrameCall <- list(formula=formula,data=data)
  response <- model.response(do.call("model.frame",modelFrameCall))

  callGeeglm <- call
  callGeeglm$family <- "binomial"
  callGeeglm$formula <- cbind(response,N-response) ~ modelMat-1
  if(any(names(callGeeglm)=="rho.fixed")) {
    callGeeglm <- callGeeglm[-which(names(callGeeglm)=="rho.fixed")]
  }
  callGeeglm <- callGeeglm[-which(names(callGeeglm)=="N")]
  fit <- do.call("geeglm",callGeeglm,envir=env)

  betas <- matrix(coef(fit),ncol=1)

  for(i in 2:10) {
    #If rho doesn't need to be update then break
    if(rho.fixed==TRUE & i > 2) {
      break
    }

    rhoOld <- rho
    #update rho
    if(rho.fixed==FALSE) {
      yhat <- predictBB(fit,N)
      residuals <- data$count - yhat

      rho <- lm(residuals^2~ offset(yhat*(1-yhat/N)) + I(yhat*(1-yhat/N)*(N-1)))$coefficients[2]
      if(rho < 0 | rho >=1) {
        warning("estimate rho not in [0,1)")
        rho <- min(max(rho,0),.9995)
        if(rho==0) break
      }
      rhos <- c(rhos,rho)
    }

    call$init.beta <- as.vector(betas[,i-1])
    try(fit <- do.call("geem",call,envir=env))
    betas <- cbind(betas,coef(fit))

    varArgs$rho <- rho
    if(abs(rhoOld-rho)<10^-6 | rho==0) break
  }

  fit$rhos <- as.vector(rhos)
  betas <- data.frame(betas)
  names(betas) <- paste("rho",round(rhos,4),sep="")
  fit$betaMat <- betas
  return(fit)
}
