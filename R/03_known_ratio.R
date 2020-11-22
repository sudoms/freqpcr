# 03 known ratio

#' Log-likelihood of getting Cq values when exact allele content is known.
#'
#' Internal function to return the log-likelihood getting the four Cq measurements under true allele frequency for each bulk sample is known. As the largest difference from \code{\link{freqpcr_loglike}()}, this function does not take account for the uncertainty in individual DNA yield.
#' @param X A numeric vector that stores the current parameter sizes of \code{meanDNA}, \code{targetScale}, \code{baseChange}, \code{sdMeasure}, \code{zeroAmount}, and \code{EPCR} in log scale.
#' @param target0,target1,housek0,housek1 Measured Cq values. Numeric vectors having the same length as \code{A} and \code{trueY}.
#' @inheritParams knownqpcr
#' @return A scalar of the log likelihood.
#' @family likelihood definitions
knownqpcr_loglike <- function(X, A, trueY, target0, target1, housek0, housek1) {
    meanDNA <- exp(X[1])
    targetScale <- exp(X[2])
    baseChange <- exp(X[3])
    sdMeasure <- exp(X[4])
    zeroAmount <- exp(X[5])
    EPCR <- exp(X[6])
    ans <- c(LogLike=0.0) # (scalar) returning value of log-likelihood.

    for(j in c(1:length(A))) {
        xR <- meanDNA*A[j]*trueY[j] # as frequency is 'true', no stochasticity in DNA quantity.
        xS <- meanDNA*A[j]*(1-trueY[j])
        # housek: housekeeping gene, target: gene of interest. w: R+S, d: R+zS
        mar <-  ( log(dnorm(housek0[j], mean=-log(xS+xR)/log(1+EPCR), sd=sdMeasure))
                + log(dnorm(target0[j], mean=-log((xS+xR)*targetScale)/log(1+EPCR), sd=sdMeasure))
                + log(dnorm(housek1[j], mean=-log((xS+xR)*baseChange)/log(1+EPCR), sd=sdMeasure))
                + log(dnorm(target1[j], mean=-log((zeroAmount*xS+xR)*targetScale*baseChange)/log(1+EPCR),
                        sd=sdMeasure))
                )
        ans <- ans + mar
    }
    return(ans)
}



#' Estimate auxiliary parameters using the bulk samples with known allele ratios.
#'
#' Function to estimate the auxiliary experimental parameters using the sample solutions where the exact allele ratios (the argument \code{trueY}) are known.
#' @param A Sample sizes as a numeric vector. It is the counterpart of the \code{N} argument in \code{\link{freqpcr}()}, but the elements of \code{A} in \code{\link{knownqpcr}()} are not restricted to integer. Including the case you have arranged the artificial sample solution with arbitrary DNA concentration, it is convenient to determine \code{A} in order to ensure that the equivalent amount of DNA extracted from an individual is approximately 1, though it is not a prerequisite and you virtually can put any positive numeric.
#' @param trueY A numeric vector of the same length as \code{A}. The size of \code{trueY[i]} signifies the true frequency of the mutant allele in the \emph{i}th bulk sample. The values must be between 0 and 1.
#' @param target0,target1,housek0,housek1 Measured Cq values. Numeric vectors having the same length as \code{A} and \code{trueY}.
#' @param XInit A named vector specifying the initial sizes of the auxiliary parameters in the optimization. Defined in the natural log scale e.g. \code{zeroAmount = -5} corresponds to the residue rate of \code{exp(-5)} = 0.007. It is highly recommended to keep the default value.
#' @param method A string specifying the optimization algorithm used in \code{\link[stats]{optim}()}. The default is \code{BFGS}, which is plausible in most situation.
#' @param trace Non-negative integer. If positive, \code{\link[stats]{optim}()} outputs trace information. The default is 0 (no information).
#' @param report The frequency of reports if \code{trace} is positive. Defaults to every 10 iterations.
#' @param quiet Suppress the output of the function sent to stdout? The default is \code{FALSE}.
#' @inheritParams freqpcr
#' @return A table containing the estimated values for the following parameters:
#' \cr
#' \code{meanDNA} is the template DNA concentration.
#' \cr
#' \code{targetScale} (\eqn{\delta_{T}}) is the relative template DNA amount of the target to the houskeeping loci.
#' \cr
#' \code{baseChange} (\eqn{\delta_{B}}) is the change rate in the amount of the internal reference after digestion (in the RED-delta delta Cq method).
#' \cr
#' \code{sdMeasure} (\eqn{\sigma_{c}}) is the measurement error on each of the four Cq values.
#' \cr
#' \code{EPCR} (\eqn{\eta}) is the amplification efficiency per PCR cycle.
#' @examples
#' # First, make a dummy Cq dataset of six bulk DNA samples,
#' # each of which comprises eight haploid individuals.
#' ntrap <- 6
#' npertrap <- 8
#' dmy_cq <- makeCqList(   rand.seed=1, P=0.75, K=2, ntrap=ntrap, npertrap=npertrap,
#'                         scaleDNA=1e-07, targetScale=1.5, baseChange=0.3, EPCR=0.95,
#'                         zeroAmount=1e-3, sdMeasure=0.3, diploid=FALSE   )
#' print(dmy_cq)
#'
#' # Then, estimate the auxiliary parameters.
#' # Note that "dmy_cq" does not include information on the parameters.
#' knownqpcr(  dmy_cq@N, rep(0.75, ntrap),
#'             dmy_cq@target0, dmy_cq@target1, dmy_cq@housek0, dmy_cq@housek1,
#'             quiet=TRUE  )
#' @export
#' @family estimation procedures
knownqpcr <- function(  A, trueY, target0, target1, housek0, housek1,
                        XInit=c(meanDNA=-10, targetScale=0, baseChange=0, sdMeasure=1, zeroAmount=-5, EPCR=0),
                        method="BFGS", pvalue=0.05, trace=0, report=10, quiet=FALSE  ) {

    Z <- optim( par=XInit, fn=knownqpcr_loglike, gr=NULL, A, trueY, target0, target1, housek0, housek1,
                method=method, control=list(fnscale=-1, trace=trace, REPORT=report), hessian=TRUE )
    SE <- sqrt(-diag(solve(Z$hessian))) # When maximized, diag(inv._of_hessian) will be negative.
    qvalue <- -qnorm(p=pvalue/2, mean=0, sd=1) # if pvalue=0.05, it returns 2.5 percentile of N(0, 1)

    result <- matrix(rep(Z$par[1:6], each=3), ncol=3, byrow=TRUE)
    result[, 2] <- result[, 2]-qvalue*SE
    result[, 3] <- result[, 3]+qvalue*SE
    result <- exp(result)
    rownames(result) <- c(  "meanDNA (DNA content per unit)",
                            "targetScale (rel. target content)",
                            "baseChange (after digestion)",
                            "SD (Cq measurement error)",
                            "zeroAmount (S-target after digestion)",
                            "EPCR (PCR multiplification EPCRicacy)"  )
    colnames(result) <- c(  "Estimate",
                            paste(deparse(100*   pvalue/2) , "%", sep=""),
                            paste(deparse(100*(1-pvalue/2)), "%", sep="")  )
    if (quiet==FALSE) {
        title <- paste( "Maximum-likelihood estimates with the two-sided ",
                        deparse(100*(1-pvalue)), "% CIs", sep="" )
        cat("\n", title, "\n", sep="")
        print(result)
    }
    return( list(report=result, obj=Z) )
}
