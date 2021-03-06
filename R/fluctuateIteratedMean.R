#' fluctuateIteratedMean
#' 
#' This function performs a fluctuation of an initial estimate of the G-computation regression at 
#' a specified time \code{t} using a call to \code{glm} (i.e., a logistic submodel) or a call to
#' \code{optim} (if bounds are specified). The structure of the function is specific to how it is 
#' called within \code{mean.tmle}.
#' In particular, \code{wideDataList} must have a very specific structure for this 
#' function to run properly. The list should consist of \code{data.frame} objects. The first should
#' have all rows set to their observed value of \code{trt}. The remaining should 
#' in turn have all rows set to each value of \code{trtOfInterest} in the
#' \code{survtmle} call. The latter will be used to obtain predictions that are then mapped
#' into the estimates of the cumulative incidence funciton at \code{t0}. 
#' Currently the code requires each \code{data.frame} to have
#' named columns for each name in \code{names(adjustVars)}, as well as a column named
#' \code{trt}. It must also have a columns named \code{Nj.Y} where j corresponds with the
#' numeric values input in \code{allJ}. These are the indicators of failure due to the various
#' causes before time \code{t} and are necessary for determining who to include in the 
#' fluctuation regression. Similarly, each \code{data.frame} should have a column call \code{C.Y} 
#' where Y is again \code{t-1}, so that right censored observations are not included in
#' the regressions. The function will fit a logistic regression with \code{Qj.star.t+1} 
#' as outcome (also needed as a column in \code{wideDataList}) with offset \code{qlogis(Qj.star.t)}
#' and number of additional covariates given by \code{length(trtOfInterest)}. These additional 
#' covariates should be columns in the each \code{data.frame} in \code{wideDataList} called 
#' \code{H.z.t} where \code{z} corresponds to a each unique value of \code{trtOfInterest}.
#' The function returns the same \code{wideDataList}, but with a column called 
#' \code{Qj.star.t} added to it, which is the fluctuated initial regression estimate evaluated
#' at the observed data points. 
#' 
#' @param wideDataList A list of \code{data.frame} objects. 
#' @param t The timepoint at which to compute the iterated mean. 
#' @param uniqtrt The values of \code{trtOfInterest} passed to \code{mean.tmle}.
#' @param whichJ Numeric value indicating the cause of failure for which regression should be 
#' computed.
#' @param allJ Numeric vector indicating the labels of all causes of failure. 
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate. Needed only because
#' the naming convention for the regression if \code{t==t0} is different than if \code{t!=t0}.
#' @param adjustVars Object of class \code{data.frame} that contains the variables to adjust 
#' for in the regression. 
#' @param bounds A list of bounds... XXX NEED MORE DESCRIPTION HERE XXX
#' @param ... Other arguments. Not currently used. 
#' 
#' @export 
#' 
#' @return The function then returns a list that is exactly the same as the input \code{wideDataList}, 
#' but with a column named \code{Qj.star.t} added to it, which is the fluctuated conditional mean of 
#' \code{Qj.star.t+1} evaluated at the each of the rows of each \code{data.frame} in 
#' \code{wideDataList}. 


fluctuateIteratedMean <- function(wideDataList, t, uniqtrt, whichJ, allJ, t0, bounds=NULL,...){
  outcomeName <- ifelse(t==t0, paste("N",whichJ,".",t0,sep=""), paste("Q",whichJ,"star.",t+1,sep=""))
  
  ## determine who to include in estimation
  include <- rep(T, nrow(wideDataList[[1]]))
  if(t!=1){
    for(j in allJ){
      # exclude previously failed subjects
      eval(parse(text=paste("include[wideDataList[[1]]$N",j,".",t-1,"==1] <- F",sep="")))
    }
    # exclude previously censored subjects
    eval(parse(text=paste("include[wideDataList[[1]]$C.",t-1,"==1] <- F",sep="")))
  }
  if(is.null(bounds)){
    wideDataList <- lapply(wideDataList, function(x,t){
      # check for 0's and 1's
      eval(parse(text=paste("x$Q",whichJ,".",t,"[x$Q",whichJ,".",t,"<.Machine$double.neg.eps] <- .Machine$double.neg.eps",sep="")))
      eval(parse(text=paste("x$Q",whichJ,".",t,"[x$Q",whichJ,".",t,">1-.Machine$double.neg.eps] <- 1-.Machine$double.neg.eps",sep="")))
      x
    },t=t)
    
    flucForm <- paste(outcomeName, "~ -1 + offset(qlogis(Q",whichJ,".",t,")) +", paste0("H",uniqtrt,".",t, collapse="+"),sep="")

    # fluctuation model
    suppressWarnings(
      flucMod <- glm(as.formula(flucForm), family="binomial",data=wideDataList[[1]][include,], start=rep(0, length(uniqtrt)))
    )
    # get predictions back
    wideDataList <- lapply(wideDataList, function(x,t){
      eval(parse(text=paste("x$Q",whichJ,"star.",t,"<- predict(flucMod, newdata=x,type='response')",sep="")))
      x
    },t=t)
    
  }else{
    cleverCovariates <- paste0("H",uniqtrt,".",t)
    
    # calculate offset term and outcome
    wideDataList <- lapply(wideDataList, function(x){
      eval(parse(text=paste("x$thisOutcome <- (x[,outcomeName] - x$l",whichJ,".",t,")",
                            "/(x$u",whichJ,".",t," - x$l",whichJ,".",t,")",sep="")))
      eval(parse(text=paste("x$thisScale <- x$u",whichJ,".",t," - wideDataList[[1]]$l",whichJ,".",t,sep="")))
      
      eval(parse(text=paste("x$Qtilde",whichJ,".",t,
                            " <- x$N",whichJ,".",t-1," + (1-(x$NnotJ.",t-1,"+ x$N",whichJ,".",t-1,
                            ")) * (x$Q",whichJ,".",t," - x$l",whichJ,".",t,")/x$thisScale", 
                            sep="")))
      
      eval(parse(text=paste("x$Qtilde",whichJ,".",t,"[x$Qtilde",whichJ,".",t,"==0] <- .Machine$double.neg.eps",sep="")))
      eval(parse(text=paste("x$Qtilde",whichJ,".",t,"[x$Qtilde",whichJ,".",t,"==1] <- 1-.Machine$double.neg.eps",sep="")))
      
      x$thisOffset <- 0
      eval(parse(text=paste("x$thisOffset[x$NnotJ.",t-1," + x$N",whichJ,".",t-1,"==0] <- ",
                            "qlogis(x$Qtilde",whichJ,".",t,"[x$NnotJ.",t-1," + x$N",whichJ,".",t-1,"==0])",
                            sep="")))
      x
    })
    
    if(length(cleverCovariates)>1){
      #           fluc.mod <- optim(par=rep(0,length(cleverCovariates)), 
      #                             fn=LogLikelihood.offset, 
      #                             Y=wideDataList[[1]]$thisOutcome[include], 
      #                             H=as.matrix(Diagonal(x=wideDataList[[1]]$thisScale[include])%*%
      #                                           as.matrix(wideDataList[[1]][include,cleverCovariates])),
      #                             offset=wideDataList[[1]]$thisOffset[include],
      #                             method="BFGS",gr=grad.offset,
      #                             control=list(reltol=1e-12, maxit=50000))
      #           
      fluc.mod <- optim(par=rep(0,length(cleverCovariates)), 
                        fn=LogLikelihood.offset, 
                        Y=wideDataList[[1]]$thisOutcome[include], 
                        H=as.matrix(wideDataList[[1]][include,cleverCovariates]),
                        offset=wideDataList[[1]]$thisOffset[include],
                        method="BFGS",gr=grad.offset,
                        control=list(reltol=1e-7, maxit=50000))
    }else{
      fluc.mod <- optim(par=rep(0,length(cleverCovariates)), 
                        fn=LogLikelihood.offset, 
                        Y=wideDataList[[1]]$thisOutcome[include], 
                        H=as.matrix(Diagonal(x=wideDataList[[1]]$thisScale[include])%*%
                                      as.matrix(wideDataList[[1]][include,cleverCovariates])),
                        offset=wideDataList[[1]]$thisOffset[include],
                        method="Brent",lower=-1000,upper=1000,
                        control=list(reltol=1e-7, maxit=50000))
    }
    
    if(fluc.mod$convergence!=0){
      stop("fluctuation convergence failure")
    }else{
      beta <- fluc.mod$par
      
      #           wideDataList <- lapply(wideDataList, function(x){
      #             eval(parse(text=paste("x$Q",whichJ,"star.",t,
      #                                   " <- x$N",whichJ,".",t-1," + (1-(x$NnotJ.",t-1,"+ x$N",whichJ,".",t-1,
      #                                   ")) * (plogis(x$thisOffset + ",
      #                                   "as.matrix(Diagonal(x=x$thisScale)%*%as.matrix(x[,cleverCovariates]))%*%as.matrix(beta))",
      #                                   "*x$thisScale + x$l",whichJ,".",t,")",sep="")))
      wideDataList <- lapply(wideDataList, function(x){
        eval(parse(text=paste("x$Q",whichJ,"star.",t,
                              " <- x$N",whichJ,".",t-1," + (1-(x$NnotJ.",t-1,"+ x$N",whichJ,".",t-1,
                              ")) * (plogis(x$thisOffset + ",
                              "as.matrix(x[,cleverCovariates])%*%as.matrix(beta))",
                              "*x$thisScale + x$l",whichJ,".",t,")",sep="")))
        x 
      })
    }
  }
  wideDataList
}
