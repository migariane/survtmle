% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mean.tmle.R
\name{mean.tmle}
\alias{mean.tmle}
\title{mean.tmle}
\usage{
\method{mean}{tmle}(ftime, ftype, trt, t0 = max(ftime), incidence = TRUE,
  adjustVars = NULL, SL.ftime = NULL, SL.ctime = NULL, SL.trt = NULL,
  glm.ftime = NULL, glm.ctime = NULL, glm.trt = "1", returnIC = TRUE,
  returnModels = FALSE, ftypeOfInterest = "all", trtOfInterest = "all",
  t0.bin = TRUE, bounds = NULL, verbose, ...)
}
\arguments{
\item{ftime}{A numeric vector of failure times. Missing values are not supported. Right 
censored observations should have a non-missing value for \code{ftime} while 
\code{ftype} should be set to 0.}

\item{ftype}{A numeric vector indicating the type of failure with type 0 reserved for 
right censored observations. Each unique value will be treated as an (unordered) separate
type of failure.}

\item{trt}{A numeric vector indicating observed treatment assignment. Each unique value will 
be treated as an (unordered) separate type of treatment. Currently, only two unique values of 
\code{trt} are supported.}

\item{t0}{The time at which to return cumulative incidence estimates. By default this is set
to \code{max(ftime)}.}

\item{incidence}{If \code{TRUE} the function return estimates of cumulative incidence. If 
\code{FALSE} the function returns esimtates of survival.}

\item{adjustVars}{A data.frame of adjustment variables that will be used in estimating the 
conditional treatment, censoring, and failure (hazard or conditional mean) probabilities.}

\item{SL.ftime}{A character vector or list specification to be passed to the \code{SL.library} argument 
in the call to \code{SuperLearner} for the outcome regression (either cause-specific hazards or 
condtional mean). See \code{?SuperLearner} for more information on how to specify valid 
\code{SuperLearner} libraries. It is expected that the wrappers used in the library will play nicely
with the input variables, which will be called \code{"trt"} and \code{names(adjustVars)}.}

\item{SL.ctime}{A character vector or list specification to be passed to the \code{SL.library} argument 
in the call to \code{SuperLearner} for the estimate of the conditional hazard for censoring. It is expected 
that the wrappers used in the library will play nicely with the input variables, which will be 
called \code{"trt"} and \code{names(adjustVars)}.}

\item{SL.trt}{A character vector or list specification to be passed to the \code{SL.library} argument 
in the call to \code{SuperLearner} for the estimate of the conditional probabiltiy of treatment.
It is expected that the wrappers used in the library will play nicely
with the input variables, which will be \code{names(adjustVars)}.}

\item{glm.ftime}{A character specification of the right-hand side of the equation passed to the
\code{formula} option of a call to \code{glm} for the outcome regression (either cause-specific hazards or 
condtional mean). Ignored if \code{SL.ftime != NULL}. Use \code{"trt"} to specify the treatment 
in this formula (see examples). The formula can additionally include any variables found in 
\code{names(adjustVars)}.}

\item{glm.ctime}{A character specification of the right-hand side of the equation passed to the
\code{formula} option of a call to \code{glm} for the estimate of the conditional hazard for 
censoring. Ignored if \code{SL.ctime != NULL}. Use \code{"trt"} to specify the treatment 
in this formula (see examples). The formula can additionally include any variables found in 
\code{names(adjustVars)}.}

\item{glm.trt}{A character specification of the right-hand side of the equation passed to the
\code{formula} option of a call to \code{glm} for the estimate of the conditional probabiltiy of 
treatment. Ignored if \code{SL.trt != NULL}. By default set to "1", corresponding to using empirical
estimates of each value of \code{trt}. The formula can include any variables found in 
\code{names(adjustVars)}.}

\item{returnIC}{A boolean indicating whether to return vectors of influence curve estimates. These are
needed for some post-hoc comarisons, so it is recommended to leave as \code{TRUE} (the default)
unless the user is sure these estimates will not be needed later.}

\item{returnModels}{A boolean indicating whether to return the \code{SuperLearner} or \code{glm} 
objects used to estimate the nuisance parameters. Must be set to \code{TRUE} if the user plans to 
use calls to \code{timepoints} to obtain estimates at times other than \code{t0}. See \code{?timepoints}
for more information.}

\item{ftypeOfInterest}{An input specifying what failure types to compute estimates of incidence for. 
The default value is \code{"all"}, which computes estimates for values \code{unique(ftype)}. Can alternatively
be set to a vector of values found in \code{ftype}.}

\item{trtOfInterest}{An input specifying which levels of \code{trt} are of interest. The default value
is \code{"all"}, which computes estimates for values \code{unique(trt)}. Can alternatively be set to a 
vector of values found in \code{trt}.}

\item{bounds}{A list of bounds... XXX NEED MORE DESCRIPTION HERE XXX}

\item{verbose}{A boolean indicating whether the function should print messages to indicate progress.}

\item{...}{Other arguments. Not currently used.}
}
\value{
call The call to \code{survtmle}.

est A numeric vector of point estimates -- one for each combination of \code{ftypeOfInterest}
and \code{trtOfInterest}.

var A covariance matrix of all the point estimates

meanIC The empirical mean of the efficient influence function at the estimated, fluctuated
nuisance parameters. If all goes well, each value should be small. This should be confirmed, 
particularly if \code{bounds} were specified.

ic The efficient influence function at the estimated, fluctuated nuisance parameters,
evaluated on each of the observations (summed over all times). These may be used to perform
post-hoc comparisons. More details coming.

ftimeMod If \code{returnModels=TRUE} the fit object(s) for the call to \code{glm} or 
\code{SuperLearner} for the outcome regression models. If \code{method="mean"} this will be a list
of length \code{length(ftypeOfInterest)} each of length \code{t0} (one regression for each
failure type and for each time point). If \code{method="hazard"} this will be a list of length
\code{length(ftypeOfInterest)} with one model corresponding with the (pseudo-)hazard for each
cause of failure. If \code{returnModels=FALSE}, this will equal \code{NULL}.

ctimeMod If \code{returnModels=TRUE} the fit object for the call to \code{glm} or 
\code{SuperLearner} for the censoring hazard regression model.  
If \code{returnModels=FALSE}, this will equal \code{NULL}.

trtMod If \code{returnModels=TRUE} the fit object for the call to \code{glm} or 
\code{SuperLearner} for the conditioanl probability of \code{trt} regression model. 
If \code{returnModels=FALSE}, this will equal \code{NULL}.
}
\description{
This function estimates the marginal cumulative incidence for failures of
specified types using targeted minimum loss-based estimation based on the 
G-computation representation of cumulative incidence. The function is called
by \code{survtmle} whenever \code{method="mean"} is specified. However, power 
users could, in theory, make calls directly to this function.
}
\examples{

## Single failure type examples
# simulate data
set.seed(1234)
n <- 100
trt <- rbinom(n,1,0.5)
adjustVars <- data.frame(W1 = round(runif(n)), W2 = round(runif(n,0,2)))

ftime <- round(1 + runif(n,1,4) - trt + adjustVars$W1 + adjustVars$W2)
ftype <- round(runif(n,0,1))

#' # Fit 1
# fit mean.tmle object with glm estimators for treatment, censoring, and failure
fit1 <- mean.tmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
glm.trt = "W1 + W2", 
glm.ftime = "trt + W1 + W2", glm.ctime = "trt + W1 + W2", 
method="mean", returnModels = TRUE)

}

