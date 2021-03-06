#' makeDataList
#' 
#' The function takes a \code{data.frame} of 
#' short format right-censored failure times and reshapes the long format into 
#' the wide format needed for calls to both \code{mean.tmle} and \code{hazard.tmle}. 
#' The list returned by the function will have number of entries equal to 
#' \code{length(trtOfInterest) + 1}. The first will have the number of rows for each observation
#' equal to the \code{ftime} corresponding to that observation. The subsequent entries will
#' have \code{t0} rows for each observation and will set \code{trt} column equal to each value
#' of \code{trtOfInterest} in turn.
#' 
#' @param dat The short form \code{data.frame}
#' @param J The unique values of \code{ftype} passed to \code{survtmle}. 
#' @param ntrt The number of \code{trt} values of interest. 
#' @param uniqtrt The unique values of \code{trtOfInterest} passed to \code{mean.tmle}.
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate. 
#' @param bounds XXX NEED MORE DESCRIPTION HERE XXX
#' @param ... Other arguments. Not currently used. 
#' 
#' @export 
#' 
#' @return A list of \code{data.frame} objects as described above. 


makeDataList <- function(dat, J, ntrt, uniqtrt, t0, bounds=NULL,...){
  n <- nrow(dat)
  dataList <- vector(mode="list",length=ntrt+1)
  
  # first element used for estimation
  dataList[[1]] <- dat[rep(1:nrow(dat),dat$ftime),]
  for(j in J){
    eval(parse(text=paste("dataList[[1]]$N",j," <- 0",sep="")))
    eval(parse(text=paste("dataList[[1]]$N",j,"[cumsum(dat$ftime)] <- as.numeric(dat$ftype==j)",sep="")))  
  }
  dataList[[1]]$C <- 0
  dataList[[1]]$C[cumsum(dat$ftime)] <- as.numeric(dat$ftype==0)
  
  n.row.ii <- nrow(dataList[[1]])
  row.names(dataList[[1]])[row.names(dataList[[1]]) %in% paste(row.names(dat))] <- paste(row.names(dat),".0",sep="")
  dataList[[1]]$t <- as.numeric(paste(unlist(strsplit(row.names(dataList[[1]]),".",fixed=T))[seq(2,n.row.ii*2,2)]))+1
  
  if(!is.null(bounds)){
    boundFormat <- data.frame(t=bounds[[1]]$t)
    for(j in J){
      if(paste("l",j,sep="") %in% names(bounds[[j]])){
        eval(parse(text=paste("boundFormat$l",j," <- bounds[[j]]$l",j,sep="")))
      }else{
        eval(parse(text=paste("boundFormat$l",j," <- 0",sep="")))
      }
      if(paste("u",j,sep="") %in% names(bounds[[j]])){
        eval(parse(text=paste("boundFormat$u",j," <- bounds[[j]]$u",j,sep="")))
      }else{
        eval(parse(text=paste("boundFormat$u",j," <- 1",sep="")))
      }
    }
    suppressMessages(
      dataList[[1]] <- join(x=dataList[[1]],y=boundFormat,type="left")
    )  
  }else{
    for(j in J){
      eval(parse(text=paste("dataList[[1]]$l",j," <- 0",sep="")))
      eval(parse(text=paste("dataList[[1]]$u",j," <- 1",sep="")))
    }
  }
  
  # subsequent elements used for prediction
  for(i in 1:ntrt){
    dataList[[i+1]] <- dat[sort(rep(1:nrow(dat),t0)),]
    dataList[[i+1]]$t <- rep(1:t0,n)
    for(j in J){
      typejEvents <- dat$id[which(dat$ftype==j)]
      eval(parse(text=paste("dataList[[i+1]]$N",j," <- 0",sep="")))
      eval(parse(text=paste("dataList[[i+1]]$N",j,"[dataList[[i+1]]$id %in% typejEvents &  dataList[[i+1]]$t >= dataList[[i+1]]$ftime] <- 1",sep="")))
    }
    censEvents <- dat$id[which(dat$ftype==0)]
    dataList[[i+1]]$C <- 0
    dataList[[i+1]]$C[dataList[[i+1]]$id %in% censEvents & dataList[[i+1]]$t >= dataList[[i+1]]$ftime] <- 1
    dataList[[i+1]]$trt <- uniqtrt[i]
    dataList[[i+1]]$ftime <- t0
    
    if(!is.null(bounds)){
      suppressMessages(
        dataList[[i+1]] <- join(x=dataList[[i+1]],y=boundFormat,type="left")
      ) 
    }else{
      for(j in J){
        eval(parse(text=paste("dataList[[",i,"+1]]$l",j," <- 0",sep="")))
        eval(parse(text=paste("dataList[[",i,"+1]]$u",j," <- 1",sep="")))
      }
    }
  }
  names(dataList) <- c("obs",uniqtrt)
  return(dataList)
}
