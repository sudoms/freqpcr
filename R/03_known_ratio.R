# 03 known ratio


#' @title Estimate auxiliary parameters using samples with known allele ratios.
#'
#' @description The function to estimate the auxiliary experimental parameters using DNA solutions, provided the dataset contains samples with multiple allele mixing ratios and the exact mixing ratio are known for each sample. This function is used when all replicates in the  dataset comprise the complete observations on the \eqn{2 \times 2} combinations of the qPCR conditions in a RED-\eqn{\Delta\Delta}Cq analysis: (loci for target-gene or housekeeping-gene) and (the target gene is undigested or digested with endonuclease).
#' \cr
#' The quartet of the four Cq data, \code{housek0}, \code{target0} (these two are undigested), \code{housek1}, and \code{target1} (they are digested) should be prepared as four numeric vectors having the same length. If the \eqn{2 \times 2} combination is not complete, use another function, \code{\link{knownqpcr_unpaired}()}, where you can deal with situations, for example, only a part of the undigested sample was quantified.
#' \cr
#' One more variable, \code{trueY} is needed to run the estimation process; it a numeric vector having the same length with the four Cq data. It holds the exact allele-mixing ratio for each quartet (also see the code example). Optionally, you can fix the relative DNA concentration between the replicates with a parameter vector \code{A}.
#' \cr
#' The \code{\link{knownqpcr}()} function can also deal with general \eqn{\Delta\Delta}Cq analyses. In such cases, samples with any mixing ratios should be marked as `digested samples' i.e., either of  \code{housek1} or \code{target1}, depending on the loci to be amplified. The arguments of the corresponding undigested samples, \code{housek0} and \code{target0}, are not required. The parameter \code{baseChange} (the change rate of DNA contents before/after the endonuclease digestion) is then not included in the estimation result. You can also use \code{\link{knownqpcr_unpaired}()} to analyze general \eqn{\Delta\Delta}Cq data.
#' @param housek0,target0,housek1,target1 Measured Cq values. Numeric vectors having the same length as \code{A} and \code{trueY}. Any of the values should not be duplicated (any single Cq measure must not be recycled). If the Cq dataset has missing values, there are two ways. 1) Specify missing values explicitly as NA and use \code{\link{knownqpcr}()}, or 2) Use another function \code{\link{knownqpcr_unpaired}()}.
#' \cr
#' In RED-\eqn{\Delta\Delta}Cq method, \code{housek0} and \code{target0} corresponds to the intact test samples (not digested with endonuclease) amplified with the primer sets for housekeeping- and target-loci, respectively. In general \eqn{\Delta\Delta}Cq analyses, no samples are marked as \code{housek0} or \code{target0} by the user (only \code{housek1} and \code{target1} should be input).
#' @param trueY A numeric vector having the same length as the Cq data vectors. \code{trueY[i]} signifies the exact frequency of the mutant allele in the \emph{i}th sample. The values must be between 0 and 1. To improve the estimation accuracy, y == 1 (pure mutant solution) should be added to your experimental design.
#' @param A Optionally, you can specify relative DNA content between the samples, as a numeric vector having the same length as the Cq data vectors. It is the counterpart of the \code{N} argument in \code{\link{freqpcr}()}, but an element of \code{A} is not restricted to integer. Because the concentration as a whole is also adjusted with the parameter \code{meanDNA} (see Value section), this variable should be used exclusively to reflect the relative DNA contents between the samples solutions. Otherwise, it is better to keep its default setting (1 for all replicates).
#' @param XInit A named vector specifying the initial sizes of the auxiliary parameters in the optimization. Defined in the natural log scale; e.g. \code{zeroAmount = -5} corresponds to the residue rate of \code{exp(-5)} = 0.007. It is highly recommended to keep the default values.
#' @param method A string specifying the optimization algorithm used in \code{\link[stats]{optim}()}. The default is \code{BFGS}, which is plausible in most situation.
#' @param trace Non-negative integer. If positive, \code{\link[stats]{optim}()} outputs trace information. The default is 0 (no information).
#' @param report The frequency of reports if \code{trace} is positive. Defaults to every 10 iterations.
#' @param quiet Suppress the output of the function sent to stdout?.
#' @inheritParams freqpcr
#' @return A table containing the estimated values for the following parameters:
#' \enumerate{
#'     \item\code{meanDNA} is the template DNA concentration (of housekeeping gene, per unit volume of sample solution) compared to the threshold line of PCR.
#'     \item\code{targetScale} (\eqn{\delta_{T}}) is the relative template DNA amount of the target to the houskeeping loci.
#'     \item\code{baseChange} (\eqn{\delta_{B}}) is the change rate in the DNA amount after endonuclease digestion (in RED-\eqn{\Delta\Delta}Cq method). In general \eqn{\Delta\Delta}Cq analyses (neither \code{housek0} nor \code{target0} is input), this parameter is not returned. In both cases, \code{baseChange} is not used in \code{\link{freqpcr}()}.
#'     \item\code{sdMeasure} (\eqn{\sigma_{c}}) is the measurement error (standard deviation) at each Cq value.
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
#' knownqpcr(  housek0=d.cmp$housek0, target0=d.cmp$target0,
#'             housek1=d.cmp$housek1, target1=d.cmp$target1,
#'             trueY=d.cmp$trueY, A=d.cmp$A   )
#'
#' # In general DeltaDeltaCq analyses, the experimental design will not include
#' # dedicated control samples. The function then runs without 'housek0' and 'target0'.
#' knownqpcr(  housek1=d.cmp$housek1, target1=d.cmp$target1,
#'             trueY=d.cmp$trueY, A=d.cmp$A  )
#' @export
#' @family estimation procedures
knownqpcr <- function(  housek1, target1, housek0=NULL, target0=NULL, trueY, A=rep(1, length(trueY)),
                        XInit=c(meanDNA=-10, targetScale=0, baseChange=0, sdMeasure=1, zeroAmount=-5, EPCR=0),
                        method="BFGS", pvalue=0.05, trace=0, report=10, quiet=FALSE  ) {

    arg.len <- all(length(A), length(trueY), length(housek1), length(target1))
    if (!arg.len) {
        stop(paste("Error: the arguments A, trueY, housek1, and target1 must have the same length."))
    }
    # If all(length(housek0), length(target0), length(housek1), length(target1)), then the dataset is
    # treated as a "quartet," which is the case in RED-DeltaDeltaCq analyses.
    quartet <- TRUE
    if (is.null(housek0) | is.null(target0)) {
        # If either housek0 or target0 is absent, then the dataset is treated as "duo."
        quartet <- FALSE
        # The parameter 'baseChange' is not defined for duo. It is the case in general DeltaDeltaCq analyses.
        XInit <- XInit[names(XInit) != "baseChange"]
        cat("\nEither the arguments housek0 or target0 was not specified. 'baseChange' will not be estimated.\n")
    } else {
        if (!all(length(housek1), length(target1), length(housek0), length(target0))) {
            stop(paste( "Error: when you input the RED-DeltaDeltaCq dataset in the 'paired' format, \n",
                        "all of housek0, target0, housek1 and target1 must have the same length. \n",
                        "Fill the missing data with NA, or use 'knownqpcr_unpaired' function instead.", sep="" ))
        }
    }

    if (quartet) {
        Cq.long <- c(housek0, target0, housek1, target1)
        A.long <- rep(A, 4)
        Y.long <- rep(trueY, 4)
        Gene <- c(rep(0, length(housek0)), rep(1, length(target0)), rep(0, length(housek1)), rep(1, length(target1)))
        Digest <- c(rep(0, length(housek0) + length(target0)), rep(1, length(housek1) + length(target1)))
    } else {
        Cq.long <- c(housek1, target1)
        A.long <- rep(A, 2)
        Y.long <- rep(trueY, 2)
        Gene <- c(rep(0, length(housek1)), rep(1, length(target1)))
        Digest <- rep(1, length(housek1) + length(target1))
        # If there is a data of Y == 1, it is automatically re-assigned to Digest == 0 (control sample).
        # Note that the data structure is 'duo' where the data lengths of the test sample and control can differ.
        Digest[Y.long==1] <- 0
        if (any(Y.long==1)) {
            cat("There are test samples with trueY == 1.",
                "They are assigned internally to Digest == 0 (control sample).\n", sep=" ")
        }
    }
    A.long <- A.long[!is.na(Cq.long)]
    Y.long <- Y.long[!is.na(Cq.long)]
    Gene <- Gene[!is.na(Cq.long)]
    Digest <- Digest[!is.na(Cq.long)]
    Cq.long <- Cq.long[!is.na(Cq.long)] # trimming of Cq must comes last.

    if(quiet==FALSE) {
        cat("\n", "Dataset internally used to estimation:", "\n", sep="")
        data.matrix <- cbind(Y.long, Digest, Gene, A.long, Cq.long)
        colnames(data.matrix) <- c( "Allele_ratio", "Sample_(control=0,_test=1)", "Gene_(housek=0,_target=1)",
                                    "Rel_DNA_content", "Cq" )
        print(data.matrix, max=ncol(data.matrix)*200)
        cat("\n")
    }

    if (quartet) {
        Z <- optim( par=XInit, fn=.knownqpcr_loglike, gr=NULL, A.long, Y.long, Digest, Gene, Cq.long,
                    method=method, control=list(fnscale=-1, trace=trace, REPORT=report), hessian=TRUE )
        SE <- sqrt(-diag(solve(Z$hessian))) # When maximized, diag(inv._of_hessian) will be negative.
        qvalue <- -qnorm(p=pvalue/2, mean=0, sd=1) # if pvalue=0.05, it returns 2.5 percentile of N(0, 1)

        result <- matrix(rep(Z$par[1:6], each=3), ncol=3, byrow=TRUE)
        result[, 2] <- result[, 2]-qvalue*SE
        result[, 3] <- result[, 3]+qvalue*SE
        result <- exp(result)
        rownames(result) <- c(  "meanDNA (DNA content per unit)",
                                "targetScale (rel. target content)",
                                "baseChange (after digestion)", # baseChange is only estimated in quartet
                                "SD (Cq measurement error)",
                                "zeroAmount (S-target after digestion)",
                                "EPCR (PCR multiplification efficiency)"  )
    } else {
        Z <- optim( par=XInit, fn=.knownqpcr_loglike_duo, gr=NULL, A.long, Y.long, Digest, Gene, Cq.long,
                    method=method, control=list(fnscale=-1, trace=trace, REPORT=report), hessian=TRUE )
        SE <- sqrt(-diag(solve(Z$hessian))) # When maximized, diag(inv._of_hessian) will be negative.
        qvalue <- -qnorm(p=pvalue/2, mean=0, sd=1) # if pvalue=0.05, it returns 2.5 percentile of N(0, 1)

        result <- matrix(rep(Z$par[1:5], each=3), ncol=3, byrow=TRUE)
        result[, 2] <- result[, 2]-qvalue*SE
        result[, 3] <- result[, 3]+qvalue*SE
        result <- exp(result)
        rownames(result) <- c(  "meanDNA (DNA content per unit)",
                                "targetScale (rel. target content)",
                                "SD (Cq measurement error)",
                                "zeroAmount (S-target after digestion)",
                                "EPCR (PCR multiplification efficiency)"  )
    }
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



#' @title Estimate auxiliary parameters when the sample pairs are incomplete.
#'
#' @description A variant of \code{\link{knownqpcr}()} used when the qPCR analyses are only conducted for part of the sample sets (i.e., the combination of the four Cq measurements has missing values). To deal with the situation, the function accepts the observed Cq values concatenated into a single vector (the argument \code{Cq}) accompanied with the experimental conditions (the arguments \code{Digest} and \code{Gene}). Their exact allele mixing ratios are known as \code{trueY}.
#' @param Digest Numeric vector having the same length as \code{Gene}, \code{trueY}, and \code{Cq}. In the RED-\eqn{\Delta\Delta}Cq method, it specify the test sample is intact (= 0) or digested with endonuclease (= 1). In general \eqn{\Delta\Delta}Cq analyses, all elements must be marked 1, which trigger the exclusion of \code{baseChange} from the estimation process. NA is not allowed.
#' @param Gene Numeric vector that specify the amplified gene is housekeeping (= 0) or target (= 1). NA is not allowed.
#' @param trueY A numeric vector. \code{trueY[i]} signifies the exact frequency of the mutant allele in the \emph{i}th sample. The values must be between 0 and 1, and NA is not allowed. To improve the estimation accuracy, y == 1 (pure mutant solution) should be added to your experimental design.
#' @param Cq Measured Cq values. NAs are allowed, and then the data vectors are trimmed internally.
#' @param A Optionally, you can specify relative DNA content between the samples, as a numeric vector having the same length as \code{Digest}, \code{Gene}, \code{trueY}, and \code{Cq}. It is the counterpart of the \code{N} argument in \code{\link{freqpcr}()}, but an element of \code{A} is not restricted to integer. Because the concentration as a whole is also adjusted with the parameter \code{meanDNA} (see Value section), this variable should be used exclusively to reflect the relative DNA contents between the samples solutions. Otherwise, it is better to keep its default setting (1 for all replicates).
#' @param XInit A named vector specifying the initial sizes of the auxiliary parameters in the optimization. Defined in the natural log scale; e.g. \code{zeroAmount = -5} corresponds to the residue rate of \code{exp(-5)} = 0.007. It is highly recommended to keep the default values.
#' @param method A string specifying the optimization algorithm used in \code{\link[stats]{optim}()}. The default is \code{BFGS}, which is plausible in most situation.
#' @param trace Non-negative integer. If positive, \code{\link[stats]{optim}()} outputs trace information. The default is 0 (no information).
#' @param report The frequency of reports if \code{trace} is positive. Defaults to every 10 iterations.
#' @param quiet Suppress the output of the function sent to stdout?.
#' @inheritParams freqpcr
#' @return A table containing the estimated parameter values, which is same as \code{\link{knownqpcr}()}.
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
                                method="BFGS", pvalue=0.05, trace=0, report=10, quiet=FALSE ) {

    arg.len <- all(length(Digest), length(Gene), length(trueY), length(Cq))
    if (!arg.len) {
        stop(paste("Error: the arguments Digest, Gene, trueY, and Cq must have the same length."))
    }
    quartet <- TRUE
    if (min(Digest)>0) {
        # If all the samples are marked as "digested," then the dataset is treated as "duo."
        quartet <- FALSE
        # The parameter 'baseChange' is not defined for duo. It is the case in general DeltaDeltaCq analyses.
        XInit <- XInit[names(XInit) != "baseChange"]
        cat("\nAs all samples are marked `digested', `baseChange' will not be estimated.\n")
    }

    Digest <- Digest[!is.na(Cq)]
    Gene <- Gene[!is.na(Cq)]
    Y.long <- trueY[!is.na(Cq)]
    A.long <- A[!is.na(Cq)]
    Cq.long <- Cq[!is.na(Cq)] # trimming of Cq must comes last.

    # ここは必要なのか？
    if (!quartet) {
        # If there is a data of Y == 1, it is automatically re-assigned to Digest == 0 (control sample).
        # Note that the data structure is 'duo' where the data lengths of the test sample and control can differ.
        Digest[Y.long==1] <- 0
        if (any(Y.long==1)) {
            cat("There are test samples with trueY == 1.",
                "They are assigned internally to Digest == 0 (control sample).\n", sep=" ")
        }
    }

    if(quiet==FALSE) {
        cat("\n", "Dataset internally used to estimation:", "\n", sep="")
        data.matrix <- cbind(Y.long, Digest, Gene, A.long, Cq.long)
        colnames(data.matrix) <- c( "Allele_ratio", "Sample_(control=0,_test=1)", "Gene_(housek=0,_target=1)",
                                    "Rel_DNA_content", "Cq" )
        print(data.matrix, max=ncol(data.matrix)*200)
        cat("\n")
    }

    if (quartet) {
        Z <- optim( par=XInit, fn=.knownqpcr_loglike, gr=NULL, A.long, Y.long, Digest, Gene, Cq.long,
                    method=method, control=list(fnscale=-1, trace=trace, REPORT=report), hessian=TRUE )
        SE <- sqrt(-diag(solve(Z$hessian))) # When maximized, diag(inv._of_hessian) will be negative.
        qvalue <- -qnorm(p=pvalue/2, mean=0, sd=1) # if pvalue=0.05, it returns 2.5 percentile of N(0, 1)

        result <- matrix(rep(Z$par[1:6], each=3), ncol=3, byrow=TRUE)
        result[, 2] <- result[, 2]-qvalue*SE
        result[, 3] <- result[, 3]+qvalue*SE
        result <- exp(result)
        rownames(result) <- c(  "meanDNA (DNA content per unit)",
                                "targetScale (rel. target content)",
                                "baseChange (after digestion)", # baseChange is only estimated in quartet
                                "SD (Cq measurement error)",
                                "zeroAmount (S-target after digestion)",
                                "EPCR (PCR multiplification efficiency)"  )
    } else {
        Z <- optim( par=XInit, fn=.knownqpcr_loglike_duo, gr=NULL, A.long, Y.long, Digest, Gene, Cq.long,
                    method=method, control=list(fnscale=-1, trace=trace, REPORT=report), hessian=TRUE )
        SE <- sqrt(-diag(solve(Z$hessian))) # When maximized, diag(inv._of_hessian) will be negative.
        qvalue <- -qnorm(p=pvalue/2, mean=0, sd=1) # if pvalue=0.05, it returns 2.5 percentile of N(0, 1)

        result <- matrix(rep(Z$par[1:5], each=3), ncol=3, byrow=TRUE)
        result[, 2] <- result[, 2]-qvalue*SE
        result[, 3] <- result[, 3]+qvalue*SE
        result <- exp(result)
        rownames(result) <- c(  "meanDNA (DNA content per unit)",
                                "targetScale (rel. target content)",
                                "SD (Cq measurement error)",
                                "zeroAmount (S-target after digestion)",
                                "EPCR (PCR multiplification efficiency)"  )
    }
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
#' @description Internal function to return the log-likelihood getting the four Cq measurements under true allele frequency for each bulk sample is known. This is a variant for the 'duo' structure dataset and \code{baseChange} is not calculated.
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
    obs.housek1 <- -log(  xS+xR  )/log(1+EPCR) # as we set relative DNA content by A[i], baseChange is not needed.
    obs.target1 <- -log( (zeroAmount*xS+xR)*targetScale )/log(1+EPCR)
    obs.matrix <- cbind(obs.housek0, obs.target0, obs.housek1, obs.target1)

    for(j in c(1:length(A))) {
        mar <- dnorm(Cq[j], mean=obs.matrix[j, 2*Digest[j] + Gene[j] + 1], sd=sdMeasure, log=TRUE)
        ans <- ans + mar
    }
    return(ans)
}








## .knownqpcr_loglike_paired was deprecated as .knownqpcr_loglike was also adopted in paired data.

#' @title Log-likelihood of obtaining Cq values when exact allele content is known.
#'
#' @description The internal function to calculate the log-likelihood by obtaining the sets of four Cq measurements under true allele frequency for each bulk sample is known. As the largest difference from \code{\link{.freqpcr_loglike}()}, this function does not take account for the uncertainty of the individual DNA yield.
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
