#' getHazardInfluenceCurve
#' 
#' This function computes the hazard-based efficient influence curve at the final
#' estimate of the fluctuated cause-specific hazard functions and evaluates 
#' it on the observed data. The influence-function is computed on the long-format data
#' but is subsequently summed over all timepoints for each observation and the function
#' returns a new short form data set with columns added corresponding to the sum over
#' all time points of the estimated efficient influence function evaluated at that observation.
#' 
#' @param dataList A list of \code{data.frame} objects. See \code{?makeDataList} for more information.
#' @param dat A \code{data.frame} in short form. See \code{?makeDataList} for more information.
#' @param allJ Numeric vector indicating the labels of all causes of failure. 
#' @param ofInterestJ Numeric vector indicating \code{ftypeOfInterest} that was passed to 
#' \code{hazard.tmle}. 
#' @param nJ The number of unique failure types. 
#' @param uniqtrt The values of \code{trtOfInterest} passed to \code{mean.tmle}.
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate. 
#' @param verbose A boolean indicating whether the function should print messages to indicate progress.
#' @param ... Other arguments. Not currently used. 
#' 
#' @return An object of class \code{data.frame} with columns \code{D.jX.zZ} added for each value
#' of X in \code{ofInterestJ} and each value of Z in \code{uniqtrt}. These are the sum over
#' all time points of the estimated efficient influence function evaluated at that observation.
#' 
#' @export


getHazardInfluenceCurve <- function(
  dataList, dat, allJ, ofInterestJ, nJ, uniqtrt, t0, verbose, ...
){
  for(z in uniqtrt){
    for(j in ofInterestJ){
      eval(parse(text=paste("dat$margF",j,".z",z,".t0 <- mean(dataList[[1]]$F",j,".z",z,".t0[dataList[[1]]$t==1])",sep="")))
      eval(parse(text=paste("dat$F",j,".z",z,".t0 <- dataList[[1]]$F",j,".z",z,".t0[dataList[[1]]$t==1]",sep="")))      
      thisD <- NULL
      for(jTild in allJ){
        thisD <- eval(parse(text=paste("cbind(thisD, dataList[[1]]$H",j,".j",ifelse(jTild==j,"Self","NotSelf"),".z",z,"/(1-dataList[[1]]$hazNot",j,")*(dataList[[1]]$N",jTild," - dataList[[1]]$Q",jTild,"Haz))",sep="")))
      }
      eval(parse(text=paste("dat$D.j",j,".z",z," <- unlist(by(rowSums(thisD), dataList[[1]]$id, FUN=sum)) + ",
                            "dat$F",j,".z",z,".t0 - dat$margF",j,".z",z,".t0",sep="")))
    }
  }
  
  dat  
}
