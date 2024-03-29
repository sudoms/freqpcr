# 03 known ratio


# todo: introduce the random error structure for the Cq quartet (in which scale?).


#' @title Estimate auxiliary parameters using samples with known allele ratios.
#'
#' @description The function to estimate the auxiliary experimental parameters using DNA solutions, provided the dataset contains samples with multiple allele mixing ratios and the exact mixing ratio are known for each sample. This function is used when all replicates in the  dataset comprise the complete observations on the \eqn{2 \times 2} combinations of the qPCR conditions in a RED-\eqn{\Delta\Delta}Cq analysis: (loci for target or housekeeping genes) and (the target gene is undigested or digested with endonuclease). The quartet of the four Cq data, \code{housek0}, \code{target0} (these two are undigested samples amplified with housekeeping and target genes, respectively), \code{housek1}, and \code{target1} (digested samples) should be prepared as four numeric vectors having the same length.
#' \cr
#' One more variable, \code{trueY} is needed to run the estimation process; it a numeric vector having the same length with the four Cq data. It holds the exact allele-mixing ratio for each quartet (also see the code example). Optionally, you can adjust the relative DNA concentration between the replicates with a parameter vector \code{A}.
#' \cr
#' Since version 0.3.2, the \code{\link{knownqpcr}()} function can also deal with general \eqn{\Delta\Delta}Cq analyses. In such cases, samples with any mixing ratios are generally marked as `digested samples' i.e., either of  \code{housek1} or \code{target1}, depending on the loci to be amplified. The arguments of the corresponding undigested samples, \code{housek0} and \code{target0}, must not be specified by the user. Then, the parameter \code{baseChange} (\eqn{\delta_{B}}: the change rate of DNA contents before/after the endonuclease digestion) is not included in the estimation result.
#' @param housek0,target0,housek1,target1 Measured Cq values. Numeric vectors having the same length as \code{trueY}. Values must not be duplicated (any single Cq measure must not be recycled); if the dataset has missing Cq values, there are two ways. 1) Fill the missing values explicitly with NA and use \code{\link{knownqpcr}()}, or 2) use another function \code{\link{knownqpcr_unpaired}()}, which accepts a `long' format dataset.
#' \cr
#' In RED-\eqn{\Delta\Delta}Cq method, \code{housek0} and \code{target0} corresponds to the intact test samples (not digested with endonuclease) amplified with the housekeeping- and target-loci, respectively. In general \eqn{\Delta\Delta}Cq analyses, \code{housek0} and \code{target0} are absent, and only test samples (\code{housek1} and \code{target1}) are input by the user. At a first glance, the test samples seem to be unaffected by endonuclease when \code{trueY[i] == 1}, but they must also be input as \code{housek1} or \code{target1} because in fact their Cq values are affected by \code{baseChange}.
#' @param trueY A numeric vector having the same length as the Cq data. \code{trueY[i]} signifies the exact allele frequency in the \emph{i}th sample. The values must be between 0 and 1, and \code{NA} is not allowed. To improve the estimation accuracy, it is better to include the settings y == 0 (pure S solution) and y == 1 (pure R) in your dataset.
#' @param A Optionally, you can specify relative DNA content between the samples, as a numeric vector having the same length as the Cq data. If present, \code{A} must not include missing values. It is the counterpart of the \code{N} argument in \code{\link{freqpcr}()}, whereas an element of \code{A} is not restricted to integer. Because the concentration as a whole is also adjusted with the parameter \code{meanDNA} (see Value section), \code{A} is used exclusively to reflect the relative contents between the sample solutions. Otherwise, \code{A} should be left unspecified (the default is 1 for all replicates).
#' @param XInit Optionally, the named vector specifies the initial parameter values to be optimized. Defined in the natural log scale; e.g. \code{zeroAmount = -5} corresponds to the residue rate \code{exp(-5)} = 0.007. Keeping the default is highly recommended.
#' @param method A string specifying the optimization algorithm used in \code{\link[stats]{optim}()}. The default is \code{BFGS}, which is plausible in most situation.
#' @param trace Non-negative integer. If positive, \code{\link[stats]{optim}()} outputs trace information. The default is 0 (no information).
#' @param report The frequency of reports if \code{trace} is positive. Defaults to every 10 iterations.
#' @param verbose Send messages to stdout? Default is FALSE.
#' @inheritParams freqpcr
#' @return A table containing the estimated values for the following parameters:
#' \enumerate{
#'     \item\code{meanDNA} is the template DNA concentration (of housekeeping gene, per unit volume of sample solution) compared to the threshold line of PCR.
#'     \item\code{targetScale} (\eqn{\delta_{T}}) is the relative template DNA amount of the target to the houskeeping loci.
#'     \item\code{baseChange} (\eqn{\delta_{B}}) is the change rate in the DNA amount after endonuclease digestion in RED-\eqn{\Delta\Delta}Cq method. In general \eqn{\Delta\Delta}Cq analyses (neither \code{housek0} nor \code{target0} is input), this parameter is not returned. In both cases, \code{baseChange} is not required to run \code{\link{freqpcr}()}.
#'     \item\code{sdMeasure} (\eqn{\sigma_{c}}) is the measurement error (standard deviation) at each Cq value.
#'     \item\code{zeroAmount} (\eqn{z}) is the ratio of non-target allele amplified in qPCR (see the argument list of \code{\link{freqpcr}()}).
#'     \item\code{EPCR} (\eqn{\eta}) is the amplification efficiency per PCR cycle.
#' }
#' @examples
#' # A dummy Cq dataset: four mixing ratios with four replicates.
#' # K:2, scaleDNA:1e-11, targetScale:1.5, baseChange:0.3, zeroAmount:1e-3,
#' # sdMeasure:0.3, and EPCR:0.95. Assuming a RED-DeltaDeltaCq analyses.
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
#'
#' # In RED-DeltaDeltaCq analyses, four observations stem from one sample solution.
#' # Each argument must be specified with its name (name=source).
#' knownqpcr(  housek0=d.cmp$housek0, target0=d.cmp$target0,
#'             housek1=d.cmp$housek1, target1=d.cmp$target1,
#'             trueY=d.cmp$trueY, A=d.cmp$A, verbose=FALSE  )
#'
#' # In general DeltaDeltaCq analyses, the experimental design will not include
#' # dedicated control samples. The function then runs without 'housek0' and 'target0'.
#' knownqpcr(  housek1=d.cmp$housek1, target1=d.cmp$target1,
#'             trueY=d.cmp$trueY, A=d.cmp$A, verbose=TRUE  )
#' @export
#' @family estimation procedures
knownqpcr <- function(  housek0, target0, housek1, target1, trueY, A=rep(1, length(trueY)),
                        XInit=c(meanDNA=-10, targetScale=0, baseChange=0, sdMeasure=1, zeroAmount=-5, EPCR=0),
                        method="BFGS", pvalue=0.05, trace=0, report=10, verbose=FALSE  ) {

    arg.len <- length(unique(c(length(A), length(trueY), length(housek1), length(target1))))
    if (arg.len != 1) {
        stop("The arguments 'A', 'trueY', ('housek0', 'target0'), 'housek1', and 'target1' all must have the same length")
    }
    quartet <- TRUE
    if (missing(housek0) || missing(target0)) {
        # If either housek0 or target0 is absent, then the dataset is treated as "duo."
        quartet <- FALSE
        # The parameter 'baseChange' is not defined for duo. It is the case in general DeltaDeltaCq analyses.
        XInit <- XInit[names(XInit) != "baseChange"]
        warning("Either 'housek0' or 'target0' was not specified. 'baseChange' is not estimated", 
                immediate.=as.logical(verbose))
    } else {
        if (length(unique(c(length(housek0), length(target0), length(housek1), length(target1)))) != 1) {
            stop(paste( "As for the 'paired' format, ",
                        "'housek0', 'target0', 'housek1', and 'target1' all must have the same length.\n",
                        "  Fill the missing values with NA, or use 'knownqpcr_unpaired' function", sep="" ))
        }
    }

    dam0 <- cbind(trueY, A)
    is.na.YA <- is.na(rowSums(dam0))
    if (verbose & sum(is.na.YA) > 0) {
        warning("Either 'A' or 'trueY' contains missing value. The following samples are omitted",
                immediate.=as.logical(verbose))
        print(which(is.na.YA))
    }

    if (quartet) {
        # Homogeneity of the data length was already verified. Then trim the missing data element.
        ID <- c(1:length(trueY))[!is.na.YA]
        housek0 <- housek0[!is.na.YA]
        target0 <- target0[!is.na.YA]
        housek1 <- housek1[!is.na.YA]
        target1 <- target1[!is.na.YA]
        trueY <- trueY[!is.na.YA]
        A <- A[!is.na.YA]

        Cq.long <- c(housek0, target0, housek1, target1)
        ID.long <- rep(ID, 4)
        A.long <- rep(A, 4)
        Y.long <- rep(trueY, 4)
        Gene <- c(  rep(0, length(housek0)), rep(1, length(target0)), 
                    rep(0, length(housek1)), rep(1, length(target1))  )
        Digest <- c(rep(0, length(housek0) + length(target0)), 
                    rep(1, length(housek1) + length(target1)))
    } else {
        ID <- c(1:length(trueY))[!is.na.YA]
        housek1 <- housek1[!is.na.YA]
        target1 <- target1[!is.na.YA]
        trueY <- trueY[!is.na.YA]
        A <- A[!is.na.YA]

        Cq.long <- c(housek1, target1)
        ID.long <- rep(ID, 2)
        A.long <- rep(A, 2)
        Y.long <- rep(trueY, 2)
        Gene <- c(  rep(0, length(housek1)), rep(1, length(target1))  )
        Digest <-   rep(1, length(housek1) + length(target1))
    }

    dam <- data.frame(ID=ID.long, trueY=Y.long, Digest=Digest, Gene=Gene, A=A.long, Cq=Cq.long)
    # As Cq values allow missing (NA) data, NAs are removed here again, which is independent of 'is.na.YA'
    is.na.data <- is.na(rowSums(dam))
    dam <- dam[!is.na.data, ]
    if (verbose) {
        cat("  Internal dataset actually used for model fitting:\n")
        print(dam)
    }

    if (quartet) {
        Z <- optim( par=XInit, fn=.knownqpcr_loglike, gr=NULL,
                    dam$A, dam$trueY, dam$Digest, dam$Gene, dam$Cq,
                    method=method, control=list(fnscale=-1, trace=trace, REPORT=report), hessian=TRUE )
        SE <- sqrt(-diag(solve(Z$hessian))) # When maximized, diag(inv._of_hessian) will be negative.
        qvalue <- -qnorm(p=pvalue/2, mean=0, sd=1) # if pvalue=0.05, it returns 2.5 percentile of N(0, 1)

        result <- matrix(rep(Z$par[1:6], each=3), ncol=3, byrow=TRUE)
        result[, 2] <- result[, 2]-qvalue*SE
        result[, 3] <- result[, 3]+qvalue*SE
        result <- exp(result)
        rownames(result) <- c(  "meanDNA (DNA content per unit)",
                                "targetScale (rel. DNA quantity)",
                                "baseChange (after digestion)", # baseChange is only estimated in quartet
                                "sdMeasure (Cq measurement error)",
                                "zeroAmount (residue rate of digestion)",
                                "EPCR (amplification per cycle)"  )
    } else {
        Z <- optim( par=XInit, fn=.knownqpcr_loglike_duo, gr=NULL,
                    dam$A, dam$trueY, dam$Digest, dam$Gene, dam$Cq,
                    method=method, control=list(fnscale=-1, trace=trace, REPORT=report), hessian=TRUE )
        SE <- sqrt(-diag(solve(Z$hessian))) # When maximized, diag(inv._of_hessian) will be negative.
        qvalue <- -qnorm(p=pvalue/2, mean=0, sd=1) # if pvalue=0.05, it returns 2.5 percentile of N(0, 1)

        result <- matrix(rep(Z$par[1:5], each=3), ncol=3, byrow=TRUE)
        result[, 2] <- result[, 2]-qvalue*SE
        result[, 3] <- result[, 3]+qvalue*SE
        result <- exp(result)
        rownames(result) <- c(  "meanDNA (DNA content per unit)",
                                "targetScale (rel. DNA quantity)",
                                "sdMeasure (Cq measurement error)",
                                "zeroAmount (residue rate of digestion)",
                                "EPCR (amplification per cycle)"  )
    }
    colnames(result) <- c(  "Estimate",
                            paste(deparse(100*   pvalue/2) , "%", sep=""),
                            paste(deparse(100*(1-pvalue/2)), "%", sep="")  )

    if (verbose) {
        title <- paste( "Maximum-likelihood estimates with the two-sided ",
                        deparse(100*(1-pvalue)), "% CIs", sep="" )
        cat("\n", title, "\n", sep="")
        print(result)
        cat("\n")
    }
    return( list(report=result, obj=Z) )
}



#' @title Estimate auxiliary parameters when the sample pairs are incomplete.
#'
#' @description A variant of \code{\link{knownqpcr}()} that accepts the Cq values concatenated into a vector (the argument \code{Cq}) accompanied with the experimental conditions (the arguments \code{Digest} and \code{Gene}). Their exact allele mixing ratios are known as \code{trueY}.
#' @param Digest Numeric vector having the same length as \code{Gene}, \code{trueY}, and \code{Cq}. \code{NA} is not allowed. In the RED-\eqn{\Delta\Delta}Cq method, it specify the sample is intact (= 0) or digested with endonuclease (= 1). In other \eqn{\Delta\Delta}Cq-based analyses, there will be no control sample and all observations must be marked \code{Digest = 1}. 
#' @param Gene Numeric vector that specify each Cq measure (element of \code{Cq}) was taken with housekeeping (= 0) or target (= 1) locus. \code{NA} is not allowed.
#' @param trueY A numeric vector. \code{trueY[i]} signifies the exact frequency of the mutant allele in the \emph{i}th sample. The values must be between 0 and 1, and \code{NA} is not allowed. To improve the estimation accuracy, it is better to include the settings y == 0 (pure S solution) and y == 1 (pure R) in your dataset.
#' @param Cq Measured Cq values. This argument is a numeric vector and can contain \code{NA}s. The vector length must be the same as \code{Digest}, \code{Gene}, and \code{trueY} (i.e., missing values must be filled with \code{NA}).
#' @param A Optionally, you can specify relative DNA content between the samples, as a numeric vector having the same length as \code{Cq}. If present, \code{A} must not include missing values. It is the counterpart of the \code{N} argument in \code{\link{freqpcr}()}, whereas an element of \code{A} is not restricted to integer. Because the concentration as a whole is also adjusted with the parameter \code{meanDNA} (see Value section), \code{A} is used exclusively to reflect the relative contents between the sample solutions. Otherwise, \code{A} should be left unspecified (the default is 1 for all replicates).
#' @param XInit Optionally, the named vector specifies the initial parameter values to be optimized. Defined in the natural log scale; e.g. \code{zeroAmount = -5} corresponds to the residue rate \code{exp(-5)} = 0.007. Keeping the default is highly recommended.
#' @param method A string specifying the optimization algorithm used in \code{\link[stats]{optim}()}. The default is \code{BFGS}, which is plausible in most situation.
#' @param trace Non-negative integer. If positive, \code{\link[stats]{optim}()} outputs trace information. The default is 0 (no information).
#' @param report The frequency of reports if \code{trace} is positive. Defaults to every 10 iterations.
#' @param verbose Send messages to stdout? Default is FALSE.
#' @inheritParams freqpcr
#' @return A table containing the estimated parameter values. The format is same as \code{\link{knownqpcr}()}.
#' @examples
#' # A dummy Cq dataset: four mixing ratios with four replicates.
#' # K:2, scaleDNA:1e-11, targetScale:1.5, baseChange:0.3, zeroAmount:1e-3,
#' # sdMeasure:0.3, and EPCR:0.95. Assuming a RED-DeltaDeltaCq analyses.
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
#' # Incomplete observation dataset, prepared as the "long" format.
#' # If the undegested (Digest == 0) samples were only analyzed when trueY == 1.
#' d.long.all <- data.frame(
#'     trueY=rep(trueY, 4), Digest=c(rep(0, 16 + 16), rep(1, 16 + 16)),
#'     Gene=c(rep(0, 16), rep(1, 16), rep(0, 16), rep(1, 16)),
#'     A=rep(1, 16*4), Cq=c(housek0, target0, housek1, target1)  )
#' d.long <- d.long.all[d.long.all$Digest == 1 | d.long.all$trueY == 1, ]
#' print(d.long)
#'
#' knownqpcr_unpaired( Digest=d.long$Digest, Gene=d.long$Gene,
#'                     trueY=d.long$trueY, Cq=d.long$Cq, A=d.long$A )
#'
#' # In general DeltaDeltaCq analyses, the experimental design will not include
#' # dedicated control samples (Digest == 0).
#' d.long <- d.long.all[d.long.all$Digest == 1, ]
#' knownqpcr_unpaired( Digest=d.long$Digest, Gene=d.long$Gene,
#'                     trueY=d.long$trueY, Cq=d.long$Cq, A=d.long$A )
#' @export
#' @family estimation procedures
knownqpcr_unpaired <- function( Digest, Gene, trueY, Cq, A=rep(1, length(Cq)),
                                XInit=c(meanDNA=-10, targetScale=0, baseChange=0, sdMeasure=1, zeroAmount=-5, EPCR=0),
                                method="BFGS", pvalue=0.05, trace=0, report=10, verbose=FALSE ) {

    arg.len <- length(unique(c(length(Digest), length(Gene), length(trueY), length(Cq), length(A))))
    if (!arg.len) {
        stop("The arguments 'Digest', 'Gene', 'trueY, 'Cq', and 'A' all must have the same length")
    }
    quartet <- TRUE
    if (min(Digest, na.rm=TRUE) > 0) {
        # If all the samples are marked as "digested," then the dataset is treated as "duo."
        quartet <- FALSE
        # The parameter 'baseChange' is not defined for duo. It is the case in general DeltaDeltaCq analyses.
        XInit <- XInit[names(XInit) != "baseChange"]
        warning("As all samples are marked `digested', 'baseChange' is not estimated", 
                immediate.=as.logical(verbose))
    }

    # Trim the missing data element.
    dam <- data.frame(trueY, Digest, Gene, A, Cq)
    is.na.data <- is.na(rowSums(dam))
    dam <- dam[!is.na.data, ]
    if (verbose & sum(is.na.data) > 0) {
        warning("The data contain missing values. The following data points are omitted",
                immediate.=as.logical(verbose))
        print(which(is.na.data))
    }

    if (verbose) {
        dam2 <- dam
        colnames(dam2) <- c("Allele_ratio", "Sample_(control=0,_test=1)", "Gene_(housek=0,_target=1)",
                            "Rel_DNA_content", "Cq")
        cat("  Internal dataset actually used for model fitting:\n", sep="")
        print(dam2, max=ncol(dam2)*200)
    }

    if (quartet) {
        Z <- optim( par=XInit, fn=.knownqpcr_loglike, gr=NULL, dam$A, dam$trueY, dam$Digest, dam$Gene, dam$Cq,
                    method=method, control=list(fnscale=-1, trace=trace, REPORT=report), hessian=TRUE )
        SE <- sqrt(-diag(solve(Z$hessian))) # When maximized, diag(inv._of_hessian) will be negative.
        qvalue <- -qnorm(p=pvalue/2, mean=0, sd=1) # if pvalue=0.05, it returns 2.5 percentile of N(0, 1)

        result <- matrix(rep(Z$par[1:6], each=3), ncol=3, byrow=TRUE)
        rownames(result) <- c(  "meanDNA (DNA content per unit)",
                                "targetScale (rel. DNA quantity)",
                                "baseChange (after digestion)", # baseChange is only estimated in quartet
                                "sdMeasure (Cq measurement error)",
                                "zeroAmount (residue rate of digestion)",
                                "EPCR (amplification per cycle)"  )
        result[, 2] <- result[, 2]-qvalue*SE
        result[, 3] <- result[, 3]+qvalue*SE
        result <- exp(result)
    } else {
        Z <- optim( par=XInit, fn=.knownqpcr_loglike_duo, gr=NULL, dam$A, dam$trueY, dam$Digest, dam$Gene, dam$Cq,
                    method=method, control=list(fnscale=-1, trace=trace, REPORT=report), hessian=TRUE )
        SE <- sqrt(-diag(solve(Z$hessian))) # When maximized, diag(inv._of_hessian) will be negative.
        qvalue <- -qnorm(p=pvalue/2, mean=0, sd=1) # if pvalue=0.05, it returns 2.5 percentile of N(0, 1)

        result <- matrix(rep(Z$par[1:5], each=3), ncol=3, byrow=TRUE)
        rownames(result) <- c(  "meanDNA (DNA content per unit)",
                                "targetScale (rel. DNA quantity)",
                                "sdMeasure (Cq measurement error)",
                                "zeroAmount (residue rate of digestion)",
                                "EPCR (amplification per cycle)"  )
        result[, 2] <- result[, 2]-qvalue*SE
        result[, 3] <- result[, 3]+qvalue*SE
        result <- exp(result)
    }
    colnames(result) <- c(  "Estimate",
                            paste(deparse(100*   pvalue/2) , "%", sep=""),
                            paste(deparse(100*(1-pvalue/2)), "%", sep="")  )

    if (verbose) {
        title <- paste( "Maximum-likelihood estimates with the two-sided ",
                        deparse(100*(1-pvalue)), "% CIs", sep="" )
        cat("\n", title, "\n", sep="")
        print(result)
    }
    return( list(report=result, obj=Z) )
}



#' @title Log-likelihood of getting Cq values when exact allele mixing ratios are known.
#'
#' @description Internal function to return the log-likelihood getting the four Cq measurements under true allele frequency for each bulk sample is known.
#' @param X A numeric vector that stores the current parameter sizes of \code{meanDNA}, \code{targetScale}, \code{baseChange}, \code{sdMeasure}, \code{zeroAmount}, and \code{EPCR} in log scale.
#' @inheritParams knownqpcr
#' @return A scalar of the log likelihood.
#' @keywords internal
.knownqpcr_loglike <- function(X, A, trueY, Digest, Gene, Cq) {
    meanDNA <- exp(X[1])
    targetScale <- exp(X[2]) # relative DNA content of target gene / housekeeping gene"
    baseChange <- exp(X[3])
    sdMeasure <- exp(X[4])
    zeroAmount <- exp(X[5])
    EPCR <- exp(X[6])
    ans <- c(LogLike=0.0) # (scalar) returning value of log-likelihood.

    xR <- meanDNA*A*trueY
    xS <- meanDNA*A*(1-trueY)
    # housek: housekeeping gene, target: gene of interest. w: R+S, d: R+zS
    obs.housek0 <- -log(xS+xR)/log(1+EPCR)
    obs.target0 <- -log((xS+xR)*targetScale)/log(1+EPCR)
    obs.housek1 <- -log((xS+xR)*baseChange)/log(1+EPCR)
    obs.target1 <- -log((zeroAmount*xS+xR)*targetScale*baseChange)/log(1+EPCR)
    obs.matrix <- cbind(obs.housek0, obs.target0, obs.housek1, obs.target1)

    for(j in c(1:length(A))) {
        mar <- dnorm(Cq[j], mean=obs.matrix[j, 2*Digest[j] + Gene[j] + 1], sd=sdMeasure, log=TRUE)
        ans <- ans + mar
    }
    return(ans)
}


#' @title Log-likelihood of getting Cq values when exact allele mixing ratios are known, lacking the quartet structure.
#'
#' @description Internal function to return the log-likelihood getting the four Cq measurements under true allele frequency for each bulk sample is known. This is a variant for the 'duo' structure dataset and \code{baseChange} is absent. \code{obs.housek0} and \code{obs.target0} are defined for the parity of source code, but they actually accept no data.
#' @param X A numeric vector that stores the current parameter sizes of \code{meanDNA}, \code{targetScale}, \code{sdMeasure}, \code{zeroAmount}, and \code{EPCR} in log scale.
#' @inheritParams knownqpcr
#' @return A scalar of the log likelihood.
#' @keywords internal
.knownqpcr_loglike_duo <- function(X, A, trueY, Digest, Gene, Cq) {
    meanDNA <- exp(X[1])
    targetScale <- exp(X[2]) # relative DNA content of target gene / housekeeping gene
    sdMeasure <- exp(X[3])
    zeroAmount <- exp(X[4])
    EPCR <- exp(X[5])
    ans <- c(LogLike=0.0) # (scalar) returning value of log-likelihood.

    xR <- meanDNA*A*trueY
    xS <- meanDNA*A*(1-trueY)
    # housek: housekeeping gene, target: gene of interest. w: R+S, d: R+zS
    obs.housek0 <- -log(  xS+xR  )/log(1+EPCR)
    obs.target0 <- -log( (xS+xR)*targetScale )/log(1+EPCR)
    obs.housek1 <- -log(  xS+xR  )/log(1+EPCR) # baseChange is absent.
    obs.target1 <- -log( (zeroAmount*xS+xR)*targetScale )/log(1+EPCR)
    obs.matrix <- cbind(obs.housek0, obs.target0, obs.housek1, obs.target1)

    for(j in c(1:length(A))) {
        mar <- dnorm(Cq[j], mean=obs.matrix[j, 2*Digest[j] + Gene[j] + 1], sd=sdMeasure, log=TRUE)
        ans <- ans + mar
    }
    return(ans)
}
