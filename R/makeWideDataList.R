#' makeWideDataList
#' 
#' The function takes a \code{data.frame} and \code{list} consisting of
#' short and long format right-censored failure times. The function reshapes the long format into 
#' the wide format needed for calls to \code{mean.tmle}. The list returned by the function
#' will have number of entries equal to \code{length(trtOfInterest) + 1}. The first will 
#' contain the observed \code{trt} columns and will set \code{C.t} (the censoring counting process)
#' equal to the observed value of censoring. The subsequent entries will set \code{trt} equal to 
#' each level of \code{trtOfInterest} and set \code{C.t} to zero for everyone. 
#' 
#' @param dat The short form \code{data.frame}
#' @param allJ Numeric vector indicating the labels of all causes of failure. 
#' @param uniqtrt The values of \code{trtOfInterest} passed to \code{mean.tmle}.
#' @param adjustVars A data.frame of adjustment variables that will be used in estimating the 
#' conditional treatment, censoring, and failure (hazard or conditional mean) probabilities. 
#' @param dataList A list of long format \code{data.frame} objects. See \code{?makeDataList} for 
#' more details on formatting.
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate. 
#' @param ... Other arguments. Not currently used. 
#' 
#' @export 
#' 
#' @return A list of \code{data.frame} objects as described above. 


makeWideDataList <- function(dat, allJ, uniqtrt, adjustVars, dataList,t0,...){
  wideDataList <- vector(mode="list",length=length(dataList))
  wideDataList[[1]] <- data.frame(dat$trt, dat[,names(adjustVars)], 
                                  reshape(dataList[[2]][,!(names(dataList[[2]]) %in% 
                                                             c("trt",names(adjustVars),"ftime","ftype"))],
                                          direction="wide",timevar="t",idvar="id"))
  colnames(wideDataList[[1]])[1] <- c("trt")
  colnames(wideDataList[[1]])[2:(1+ncol(adjustVars))] <- names(adjustVars)
  # set Nj0=0 for all j -- makes things easier to run in a loop later
  eval(parse(text=paste(paste0("wideDataList[[1]]$N",allJ,".0", collapse="<-"), "<- wideDataList[[1]]$C.0 <- 0",sep="")))
  
  
  wideDataList[2:length(dataList)] <- lapply(dataList[2:length(dataList)], function(x){
    out <- data.frame(dat[,names(adjustVars)], 
                      reshape(x[,!(names(x) %in% c("trt",names(adjustVars),"ftime","ftype"))],
                              direction="wide",timevar="t",idvar="id")
                      ,row.names=NULL)
    out[,paste0("C.",1:t0)] <- 0
    names(out)[1:(ncol(adjustVars))] <- names(adjustVars)
    # set Nj0=0 for all j -- makes things easier to run in a loop later
    eval(parse(text=paste(paste0("out$N",allJ,".0", collapse="<-"), "<- out$C.0 <- 0",sep="")))
    out
  })
  names(wideDataList) <- c("obs",uniqtrt)
  
  for(z in uniqtrt) eval(parse(text=paste("wideDataList$'",z,"'$trt <-", z,sep="")))
  
  wideDataList <- lapply(wideDataList, function(x){
    # make clever covariates
    for(z in uniqtrt){
      for(t in 1:t0){
        eval(parse(text=paste("x$H",z,".",t," <- (x$trt==",z," & x$C.",t-1,"==0)/(x$G_dC.",t,"*x$g_",z,".",t,")",sep="")))
      }
      eval(parse(text=paste("x$H",z,".",0," <- (x$trt==",z,")/(x$g_",z,".",1,")",sep="")))
    }
    x
  })
  
  wideDataList
}
