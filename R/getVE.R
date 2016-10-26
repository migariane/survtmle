#' getVE
#'
#' A function that uses output from a survtmle object to compute
#' vaccine efficacy and vaccine sieve effect estimates
#'
#' @param out A \code{survtmle} object with two types of failure
#' and two treatment types.
#' @param tol Level of incidence below which not to compute VE.
#' If incidence is very small in the placebo arm the asymmetric scale of VE can
#' cause very small values of VE.
#' @return A list of data.frames whose first row is the vaccine efficacy against
#' type 1 endpoints, second row is vaccine efficacy against type 2 endpoints
#' and third row is vaccine sieve effect. In each row the first entry corresponds
#' to the effect estimate, the second and third to the lower and upper confidence
#' limits, and the fourth to the p-value for the Wald test of the
#' relevant hypothesis test (VE = 0 or VSE = 1).

getVE <- function(out, tol = 1e-10){
    fnOut <- lapply(out, function(rslt){
        # check for very small values of cumulative incidence
        # that lead to very large/small values of VE
        if(all(rslt$est[1:2] > tol)){
            # ve against type 1 endpoints
            ve1 <- 1 - rslt$est[2]/rslt$est[1]
            # relative risk on the log scale
            log.rr1 <- log(rslt$est[2]/rslt$est[1])
            # gradient used for se computation on log-scale
            a <- matrix(c(-1/rslt$est[1], 1/rslt$est[2]),nrow=2)
            # standard error on the log scale for relative risk
            se.log.rr1 <- sqrt(t(a)%*%rslt$var[1:2,1:2]%*%a)
            # confidence interval
            ve1.low <- 1 - exp(log.rr1 + 1.96*se.log.rr1)
            ve1.high <- 1 - exp(log.rr1 - 1.96*se.log.rr1)
            # two sided p-value
            p1 <- 2*pnorm(-abs(log.rr1/se.log.rr1))
        }else{
            ve1 <- ve1.low <- ve1.high <- p1 <- NA
        }

        if(all(rslt$est[3:4] > tol)){
           # type 2 proceed exactly as above
            ve2 <- 1 - rslt$est[4]/rslt$est[3]
            log.rr2 <- log(rslt$est[4]/rslt$est[3])
            a <- matrix(c(-1/rslt$est[3], 1/rslt$est[4]),nrow=2)
            se.log.rr2 <- sqrt(t(a)%*%rslt$var[(3:4),(3:4)]%*%a)
            ve2.low <- 1 - exp(log.rr2 + 1.96*se.log.rr2)
            ve2.high <- 1 - exp(log.rr2 - 1.96*se.log.rr2)
            p2 <- 2*pnorm(-abs(log.rr2/se.log.rr2))
        }else{
            ve2 <- ve2.low <- ve2.high <- p2 <- NA
        }

        if(all(rslt$est > tol)){
           # vaccine sieve effect
            sieve <- exp(log.rr2 - log.rr1)
            # on the log scale
            log.sieve <- log.rr2 - log.rr1
            # gradient used to compute standard error on the log scale
            a <- matrix(c(1/rslt$est[1], -1/rslt$est[2],
                          -1/rslt$est[3], 1/rslt$est[4]),nrow=4)
            # standard error
            se.log.sieve <- sqrt(t(a)%*%rslt$var%*%a)
            # confidence interval
            sieve.low <- exp(log.sieve -1.96*se.log.sieve)
            sieve.high <- exp(log.sieve +1.96*se.log.sieve)
            # p-value for hypothesis test
            sieve.p <- 2*pnorm(-abs(log.sieve/se.log.sieve))
        }else{
            sieve <- sieve.low <- sieve.high <- sieve.p <- NA
        }

        # return a matrix
        tmp <- data.frame(
            rbind(c(ve1, ve1.low, ve1.high, p1),
                  c(ve2, ve2.low, ve2.high, p2),
                  c(sieve, sieve.low, sieve.high, sieve.p)
        ))
        row.names(tmp) <- c("ve1","ve2","vse")
        colnames(tmp) <- c("est","ci.low","ci.high","p")
        return(tmp)
    })
    fnOut
}
