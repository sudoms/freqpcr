# 03 known ratio


#' @title Estimate auxiliary parameters using the bulk samples with known allele ratios.
#'
#' @description Function to estimate the auxiliary experimental parameters using DNA solutions, where the dataset contains the samples with multiple allele mixing ratios (the argument \code{trueY}) and the exact mixing ratio are known for each sample. This function is used when all the replicates comprise of the observations on the \eqn{2\times2} combinations of the qPCR conditions: (target- versus housekeeping-gene) and (undigested or digested with endonuclease). If the combination is not complete, use another function \code{\link{knownqpcr_unpaired}()}, where you can deal with the situation e.g., only a part of the undigested samples was quantified.
#' @param A Sample sizes as a numeric vector. It is the counterpart of the \code{N} argument in \code{\link{freqpcr}()}, but an element of \code{A} is not restricted to integer. Considering the case you have arranged the artificial sample solution with arbitrary DNA concentration, it is convenient to determine \code{A} to ensure that the equivalent amount of DNA extracted from an individual is approximately 1. As the concentration in total is also adjusted with the parameter \code{targetScale} (see Value section), this variable should be used exclusively to reflect the differences in relative concentration between samples.
#' @param trueY A numeric vector of the same length as \code{A}. \code{trueY[i]} signifies the true frequency of the mutant allele in the \emph{i}th bulk sample. The values must be between 0 and 1.
#' @param housek0,target0,housek1,target1 Measured Cq values. Numeric vectors having the same length as \code{A} and \code{trueY}. Any of the values should not be duplicated (if so, use \code{\link{knownqpcr_unpaired}()} instead).
#' @param XInit A named vector specifying the initial sizes of the auxiliary parameters in the optimization. Defined in the natural log scale; e.g. \code{zeroAmount = -5} corresponds to the residue rate of \code{exp(-5)} = 0.007. It is highly recommended to keep the default values.
#' @param method A string specifying the optimization algorithm used in \code{\link[stats]{optim}()}. The default is \code{BFGS}, which is plausible in most situation.
#' @param trace Non-negative integer. If positive, \code{\link[stats]{optim}()} outputs trace information. The default is 0 (no information).
#' @param report The frequency of reports if \code{trace} is positive. Defaults to every 10 iterations.
#' @param quiet Suppress the output of the function sent to stdout? The default is \code{FALSE}.
#' @inheritParams freqpcr
#' @return A table containing the estimated values for the following parameters:
#' \enumerate{
#'     \item\code{meanDNA} is the template DNA concentration.
#'     \item\code{targetScale} (\eqn{\delta_{T}}) is the relative template DNA amount of the target to the houskeeping loci.
#'     \item\code{baseChange} (\eqn{\delta_{B}}) is the change rate in the amount of the internal reference after digestion (in the RED-delta delta Cq method).
#'     \item\code{sdMeasure} (\eqn{\sigma_{c}}) is the measurement error on each of the four Cq values.
#'     \item\code{EPCR} (\eqn{\eta}) is the amplification efficiency per PCR cycle.
#' }
#' @examples
#' # A dummy Cq dataset: four mixing ratios (0.1, 0.25, 0.5, 1) with four replicates.
#' trueY <- c(rep(0.1, 4), rep(0.25, 4), rep(0.5, 4), rep(1, 4))
#' housek0 <- c( 19.39, 19.78, 19.28, 19.58,  18.95, 19.91, 19.66, 19.96,
#'               20.05, 19.86, 19.55, 19.61,  19.86, 19.27, 19.59, 20.21 )
#' target0 <- c( 19.16, 19.08, 19.28, 19.03,  19.17, 19.67, 18.68, 19.52,
#'               18.92, 18.79, 18.8, 19.28,   19.57, 19.21, 19.05, 19.15 )
#' housek1 <- c( 21.61, 21.78, 21.25, 21.07,  22.04, 21.45, 20.72, 21.6,
#'               21.51, 21.27, 21.08, 21.7,   21.44, 21.46, 21.5, 21.8 )
#' target1 <- c( 24.3, 24.22, 24.13, 24.13,   22.74, 23.14, 23.02, 23.14,
#'               21.65, 22.62, 22.28, 21.65,  20.83, 20.82, 20.76, 21.3 )
#' d.cmp <- data.frame(A=rep(1, 16), trueY, housek0, target0, housek1, target1)
#' print(d.cmp)
#' p.cmp <- knownqpcr( d.cmp$A, d.cmp$trueY,
#'                     d.cmp$housek0, d.cmp$target0, d.cmp$housek1, d.cmp$target1 )
#' @export
#' @family estimation procedures
knownqpcr <- function(  A, trueY, housek0, target0, housek1, target1,
                        XInit=c(meanDNA=-10, targetScale=0, baseChange=0, sdMeasure=1, zeroAmount=-5, EPCR=0),
                        method="BFGS", pvalue=0.05, trace=0, report=10, quiet=FALSE  ) {

    Z <- optim( par=XInit, fn=.knownqpcr_loglike_paired, gr=NULL, A, trueY, housek0, target0, housek1, target1,
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


#' @title Log-likelihood of getting Cq values when exact allele content is known.
#'
#' @description Internal function to calculate the log-likelihood getting the sets of four Cq measurements under true allele frequency for each bulk sample is known. As the largest difference from \code{\link{.freqpcr_loglike}()}, this function does not take account for the uncertainty in individual DNA yield.
#' @param X A numeric vector that stores the current parameter sizes of \code{meanDNA}, \code{targetScale}, \code{baseChange}, \code{sdMeasure}, \code{zeroAmount}, and \code{EPCR} in log scale.
#' @param target0,target1,housek0,housek1 Measured Cq values. Numeric vectors having the same length as \code{A} and \code{trueY}.
#' @inheritParams knownqpcr
#' @return A scalar of the log likelihood.
#' @keywords internal
.knownqpcr_loglike_paired <- function(X, A, trueY, housek0, target0, housek1, target1) {
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




#' @title Estimate auxiliary parameters when the sample pairs are incomplete.
#'
#' @description A variant of \code{\link{knownqpcr}()} to estimate the auxiliary parameters when the exact mixing ratios are known for each the sample, but the qPCR analyses are conducted for only a part of the sample sets. To deal with the cases that the combination of the four Cq measurements is incomplete (thus the data lengths of \code{housek0}, \code{target0}, \code{housek1}, and \code{target1} are different), the function accepts the observed Cq values concatenated into a single vector (the argument \code{Cq}) accompanied with the experimental conditions (the arguments \code{Gene} and \code{Digest}).
#' @param A Sample sizes as a numeric vector.
#' @param Gene Numeric vector that specify the gene is housekeeping (= 0) or target (= 1), having the same length as \code{A} and \code{trueY}.
#' @param Digest Numeric vector having the same length as \code{A} and \code{trueY}. It specify the gene is intact (= 0) or digested with endonuclease (= 1) in the RED-\eqn{\Delta\Delta}Cq method. In general \eqn{\Delta\Delta}Cq analyses, the argument means all the alleles on the target locs will be amplified (= 0) (it has ) or only the R portion will be amplified (using an R-specific primer set).
#' @param Cq Measured Cq values. Numeric vectors having the same length as \code{A} and \code{trueY}.
#' @inheritParams knownqpcr
#' @return A table containing the estimated parameter values, which is same as \code{\link{knownqpcr}()}.
#' @examples
#' # A dummy Cq dataset: four mixing ratios (0.1, 0.25, 0.5, 1) with four replicates.
#' # Settings: targetScale=1.5, baseChange=0.3, EPCR=0.95, zeroAmount=1e-3, sdMeasure=0.3
#'
#' trueY <- c(rep(0.1, 4), rep(0.25, 4), rep(0.5, 4), rep(1, 4))
#' housek0 <- c( 19.39, 19.78, 19.28, 19.58,  18.95, 19.91, 19.66, 19.96,
#'               20.05, 19.86, 19.55, 19.61,  19.86, 19.27, 19.59, 20.21 )
#' target0 <- c( 19.16, 19.08, 19.28, 19.03,  19.17, 19.67, 18.68, 19.52,
#'               18.92, 18.79, 18.8, 19.28,   19.57, 19.21, 19.05, 19.15 )
#' housek1 <- c( 21.61, 21.78, 21.25, 21.07,  22.04, 21.45, 20.72, 21.6,
#'               21.51, 21.27, 21.08, 21.7,   21.44, 21.46, 21.5, 21.8 )
#' target1 <- c( 24.3, 24.22, 24.13, 24.13,   22.74, 23.14, 23.02, 23.14,
#'               21.65, 22.62, 22.28, 21.65,  20.83, 20.82, 20.76, 21.3 )
#'
#' # Complete observation: "wide" format.
#' d.cmp <- data.frame(A=rep(1, 16), trueY, housek0, target0, housek1, target1)
#' print(d.cmp)
#' p.cmp <- knownqpcr( d.cmp$A, d.cmp$trueY,
#'                     d.cmp$housek0, d.cmp$target0, d.cmp$housek1, d.cmp$target1 )
#' print(p.cmp)
#'
#' # Incomplete observation: prepared in the "long" format.
#' # Undegested (housek0 and target0) samples are analyzed only for trueY == 1.
#' d.long <- data.frame(
#'     A=rep(1, 4+4+16+16),
#'     Gene=c(rep(0, 4), rep(1, 4), rep(0, 16), rep(1, 16)),
#'     Digest=c(rep(0, 4+4), rep(1, 16+16)),
#'     trueY=c(rep(1, 8), rep(c(rep(0.1, 4), rep(0.25, 4), rep(0.5, 4), rep(1, 4)), 2)),
#'     Repl=rep(1:4, 1+1+4+4),
#'     Cq=c(housek0[13:16], target0[13:16], housek1, target1))
#' print(d.long)
#' p.long <- knownqpcr_unpaired( d.long$A, d.long$trueY,
#'                               d.long$Gene, d.long$Digest, d.long$Cq )
#' print(p.long)
#' @export
#' @family estimation procedures
knownqpcr_unpaired <- function(  A, trueY, Gene, Digest, Cq,
                        XInit=c(meanDNA=-10, targetScale=0, baseChange=0, sdMeasure=1, zeroAmount=-5, EPCR=0),
                        method="BFGS", pvalue=0.05, trace=0, report=10, quiet=FALSE  ) {

    Z <- optim( par=XInit, fn=.knownqpcr_loglike_unpaired, gr=NULL, A, trueY, Gene, Digest, Cq,
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


#' @title Log-likelihood of getting Cq values when exact allele content is known.
#'
#' @description Internal function to return the log-likelihood getting the four Cq measurements under true allele frequency for each bulk sample is known.
#' @param X A numeric vector that stores the current parameter sizes of \code{meanDNA}, \code{targetScale}, \code{baseChange}, \code{sdMeasure}, \code{zeroAmount}, and \code{EPCR} in log scale.
#' @inheritParams knownqpcr_unpaired
#' @return A scalar of the log likelihood.
#' @keywords internal
.knownqpcr_loglike_unpaired <- function(X, A, trueY, Gene, Digest, Cq) {
    meanDNA <- exp(X[1])
    targetScale <- exp(X[2])
    baseChange <- exp(X[3])
    sdMeasure <- exp(X[4])
    zeroAmount <- exp(X[5])
    EPCR <- exp(X[6])

    ans <- c(LogLike=0.0) # (scalar) returning value of log-likelihood.

    # vector of length(A == trueY)
    xR <- meanDNA*A*trueY
    xS <- meanDNA*A*(1-trueY)
    # housek: housekeeping gene, target: gene of interest. w: R+S, d: R+zS
    mar.housek0 <- -log(xS+xR)/log(1+EPCR)
    mar.target0 <- -log((xS+xR)*targetScale)/log(1+EPCR)
    mar.housek1 <- -log((xS+xR)*baseChange)/log(1+EPCR)
    mar.target1 <- -log((zeroAmount*xS+xR)*targetScale*baseChange)/log(1+EPCR)
    mar.matrix <- cbind(mar.housek0, mar.housek1, mar.target0, mar.target1)

    for(j in c(1:length(A))) {
        mar <-  ( log(dnorm(Cq[j], mean=mar.matrix[j, 2*(Gene[j])+(Digest[j]+1)], sd=sdMeasure)) )
        ans <- ans + mar
    }
    return(ans)
}
