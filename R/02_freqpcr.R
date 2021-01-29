# 02 freqpcr

setOldClass("proc_time") # Register the S3 class as formal class

#' S4 class that stores the estimation result of \code{\link{freqpcr}()}.
#'
#' @slot report A matrix of the simultaneous parameter estimation result. The rows represent the parameters \code{P}, \code{K}, \code{targetScale} (\eqn{\delta_{T}}), \code{sdMeasure} (\eqn{\sigma_{c}}), and \code{EPCR} (\eqn{\eta}).
#' @slot obj Returned value of the optimizer function \code{\link{nlm}()} as a list.
#' @slot cal.time Calculation time of \code{\link{nlm}()}, stored as a \code{proc_time} class object.
#' @export
CqFreq <- setClass("CqFreq",
    slots = c(
        report="matrix",
        obj="list",
        cal.time="proc_time"
    )
)



#' @title Estimate population allele frequency from the set of Cq measurements.

#' @description The function estimates the population allele frequency using the dataset of Cq values measured over \code{length(N)} bulk samples, each of which has a sample size of \code{N[i]} as the number of individuals included. \code{N[i]} can be 1, meaning that every individual is analyzed separately. For the \emph{i}th sample, the four Cq values were measured as \code{housek0[i]}, \code{target0[i]}, \code{housek1[i]}, and \code{target1[i]}.
#' \cr
#' The function can estimate up to five parameters simultaneously when the Cq sets are available for more than two (bulk) samples: \code{P}, \code{K}, \code{targetScale}, \code{sdMeasure}, and \code{EPCR}. However, \code{EPCR} is not reliable in many cases the user is encouraged to know the size beforehand by using \code{\link{knownqpcr}()} or \code{\link{knownqpcr_unpaired}()} (either function will return the same result). By default, \code{EPCR} is not estimated unless \code{EPCR = NULL} is specified explicitly.
#' \cr
#' Since version 0.3.2, the user can also use an experimental `continuous model' by specifying \code{A} instead of \code{N}. That is, each sample DNA is directly extracted from the environment and the sample allele ratio \code{y} follows \code{y ~ Beta(apk, a(1-p)k)} instead of \code{y ~ Beta(mk, (n-m)k), m ~ Binomial(n, p)}, where \code{p} and \code{k} are the population allele frequency and the gamma shape parameter of the individual DNA yield, respectively. The continuous model is likely when each sample directly comes from the environment e.g., water sampling in an eDNA analysis, or cultures in petri dishes.
#' @param N Sample sizes as a numeric vector. \code{N[i]} signifies the number of individuals (both for haploidy and diploidy) contained in the \emph{i}th bulk sample. If \code{N} is not applicable (= even not 1), feed \code{A} instead of \code{N} and then the estimation process runs with the `continuous model'.
#' @param A Used instead of \code{N} in the continuous model, which is a scale factor to control the relative DNA content between samples. \code{A[i]} can take any positive value. Considering the case you have arranged each sample by e.g. water filrtation or extraction from a culture in a petri dish, it is convenient to define the unit size of \code{A[i]} == 1.0 to be same as the vessel volume (e.g. 2.0 for two petri dishs, 0.5 for half bottle of water, etc.). When neither \code{N} nor \code{A} is specified by the user, the function automatically runs as the continuous model that follows \code{y ~ Beta(pk, (1-p)k)}, where \code{A[i] = 1} is assumed for all samples.
#' @param housek0 A numeric vector. In RED-\eqn{\Delta\Delta}Cq method, \code{housek0} is the Cq values of the test sample without the restriction enzyme digestion, which is amplified with the primer set for a housekeeping gene. In a general \eqn{\Delta\Delta}Cq analyses, \code{housek0} is defined for the control sample (typically, 100\% mutant) solution, which is also amplified with the primer set for the housekeeping gene.
#' \cr
#' The following four Cq arguments, \code{housek0}, \code{target0}, \code{housek1}, and \code{target1}, all must have the same data length. If either \code{N} or \code{A} argument is present, it also must be the same length as the four Cq vectors. If the Cq dataset has missing values, they must be filled with NA so that the length of the data does not differ between the arguments.
#' @param target0 A numeric vector. In RED-\eqn{\Delta\Delta}Cq method, \code{target0[i]} signifies the measured Cq value of the \emph{i}th bulk sample without the digestion, which is treated with the primer set that amplifies both alleles, wild-type (or S: susceptible) and mutant (or R: resistant to a pesticide), on the target locus. In general \eqn{\Delta\Delta}Cq analyses, \code{target0} is the Cq values of the pure-mutant control sample, which is amplified with a mutant-allele-specific primer set.
#' @param housek1 A numeric vector. The Cq values of the test sample measured on the housekeeping gene after the restriction enzyme digestion (in RED-\eqn{\Delta\Delta}Cq method), or the test sample amplified on the housekeeping gene (in general \eqn{\Delta\Delta}Cq analyses).
#' @param target1 A numeric vector. For each test sample with unknown allele-ratio, \code{target1[i]} is defined as the Cq value for the target locus amplified after the restriction enzyme digestion (in RED-\eqn{\Delta\Delta}Cq method), or the target locus amplified with the R-allele-specific primer set (in general \eqn{\Delta\Delta}Cq analyses).
#' @param P Scalar. Population allele frequency from which the test samples are derived. Default is \code{NULL} and to be estimated. If the parameter is known, it is given as a numeric between 0 and 1.
#' @param K Scalar. The gamma shape parameter of the individual DNA yield. Default is \code{NULL} and to be estimated. If known, given as a positive numeric.
#' @param targetScale (\eqn{\delta_{T}}) Scalar. The relative template DNA amount of the target locus to the houskeeping locus. If known, given as a positive numeric.
#' @param sdMeasure (\eqn{\sigma_{c}}) Scalar. The measurement error on each Cq value following Normal(0, \eqn{\sigma_{c}^2}). If known, given as a positive numeric.
#' @param EPCR (\eqn{\eta}) Scalar. Amplification efficiency per PCR cycle. If known, given as a positive numeric. When \code{EPCR = 1}, template DNA doubles every cycle (\code{EPCR + 1 = 2}).
#' @param XInit0 Optionally the initial value for the parameter optimization can be specified, but the user is highly recommended to keep the argument as is. Unlike \code{XInit} in \code{\link{knownqpcr}()}, the argument is not directly passed to the optimizer; each value is used only when the parameter is unknown.
#' @param zeroAmount (In RED-\eqn{\Delta\Delta}Cq method) residue rate of restriction enzyme digestion, or (in general \eqn{\Delta\Delta}Cq analyses) small portion of the off-target allele on the target locus of the test sample, which will be amplified in the PCR. It needs to be always specified by the user as a number between 0 and 1, usually near 0.
#' @param beta Whether to use the beta distribution to approximate the sample allele ratio instead of specifying individual gamma distribution for each of the allelic DNA amounts? Default is \code{TRUE}, which accelerates the calculation.
#' @param diploid The target organism diploidy? Default is \code{FALSE}, assuming haploidy. Current implementation assumes i.i.d. between the amounts of R and S chromosomes owned by a heterozygote individual, which is unlikely in many animals but necessary for the calculation in a realistic time.
#' @param pvalue The two-sided confidence interval is calculated at the last iteration based on the significance level. Default is 0.05, which returns the 95\% CI (2.5 to 97.5 percentile) based on the Hessian matrix.
#' @param gradtol,steptol,iterlim Control parameters passed to \code{\link[stats]{nlm}()}. \code{gradtol} and \code{steptol} are the positive scalars giving the tolerance to terminate the algorithm and the minimum allowable relative step length. \code{iterlim} specifies the maximum number of iterations to be performed before the program is terminated. Usually 30 iterations are enough.
#' @param maxtime A positive scalar to set the maximum calculation time in seconds to abort the optimizer (and return error). The total calculation time largely depends on \code{N[i]}, the number of individuals contained in each bulk sample.
#' @param print.level \code{print.level=1} (the default) shows the initial values of the parameters and likelihood as well as the output in the last iteration. \code{print.level=2} shows the parameter values and gradients in every step. \code{print.level=0} does not output any intermediate state to R console, simply returning the result summary.
#' @param ... Additional arguments passed to the function.
#' @return Object of the S4 class \linkS4class{CqFreq}. The slot \code{report} is a matrix and each row contains the estimated parameter value with 100*(1-pvalue)\% confidence interval. The following parameters are returned:
#' \enumerate{
#'     \item\code{P}, the population allele frequency from which the test samples are derived.
#'     \item\code{K}, the gamma shape parameter of the individual DNA yield.
#'     \item\code{targetScale} (\eqn{\delta_{T}}), the relative template DNA amount of the target to the houskeeping loci.
#'     \item\code{EPCR} (\eqn{\eta}), the amplification efficiency per PCR cycle.
#'     \item\code{sdMeasure} or "Cq measurement error" (\eqn{\sigma_{c}}).
#' }
#' @section Choise of the parameters to be estimated:
#' Estimation is conducted only for parameters specified, explicitly or implicitly, as \code{NULL}. If one feeds a value e.g. \code{K=1} or \code{sdMeasure=0.24}, it is then treated as fixed parameter.
#' \cr
#' \code{EPCR} is the amplification efficiency of PCR and usually not estimable from the experiments with unknown allele ratios. You should verify the size of \code{EPCR} and \code{zeroAmount} prior to the experiment with unknown samples. The functions \code{\link{knownqpcr}()} and \code{\link{knownqpcr_unpaired}()} provide the procedure to estimate the sizes of the experimental parameters using the DNA solutions of known allele mixing ratios.
#' \cr
#' For the unknown parameters, \code{XInit0} optionally specifies initial values for the optimization using \code{\link[stats]{nlm}()} though the use of internal default is highly recommended. The specification as a fixed parameter has higher priority than the specification in \code{XInit0}. Every user-specified parameter values, fixed or unknown, must be given in linear scale (e.g. between 0 and 1 for the allele frequency); they are transformed internally to log or logit.
#' @examples
#' # Dummy Cq dataset with six bulk samples, each of which comprises of eight haploids.
#' EPCR <- 0.95; zeroAmount <- 0.0016; # True values for the two must be known.
#' P <- 0.75
#' dmy_cq <- make_dummy(   rand.seed=1, P=P, K=2, ntrap=6, npertrap=8,
#'                         scaleDNA=1e-07, targetScale=1.5, baseChange=0.3,
#'                         EPCR=EPCR, zeroAmount=zeroAmount,
#'                         sdMeasure=0.3, diploid=FALSE   )
#' print(dmy_cq)
#' # Estimation with freqpcr, where P, K, targetScale, and baseChange are marked unknown.
#' result <- freqpcr( N=dmy_cq@N, housek0=dmy_cq@housek0, target0=dmy_cq@target0,
#'                    housek1=dmy_cq@housek1, target1=dmy_cq@target1,
#'                    EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=2 )
#' print(result)
#' # Estimation with freqpcr, where you have knowledge on the size of K.
#' result <- freqpcr( N=dmy_cq@N, housek0=dmy_cq@housek0, target0=dmy_cq@target0,
#'                    housek1=dmy_cq@housek1, target1=dmy_cq@target1,
#'                    K=2,
#'                    EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=1 )
#' # Estimation with freqpcr when the model is continuous.
#' result <- freqpcr( A=dmy_cq@N, housek0=dmy_cq@housek0, target0=dmy_cq@target0,
#'                    housek1=dmy_cq@housek1, target1=dmy_cq@target1,
#'                    EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=2 )
#' @export
#' @family estimation procedures
freqpcr <- function(N=NULL, A=NULL, housek0, target0, housek1, target1,
                    P=NULL, K=NULL, targetScale=NULL, sdMeasure=NULL, EPCR=0.99,
                    XInit0=c(P=NULL, K=NULL, targetScale=NULL, sdMeasure=NULL, EPCR=NULL),
                    zeroAmount=0.0016, beta=TRUE, diploid=FALSE,
                    pvalue=0.05, gradtol=1e-4, steptol=1e-9, iterlim=100, maxtime=600, print.level=1, ...) {
    if (is.null(EPCR) & length(target1)<3) {
        stop(paste("Error: when EPCR is unknown, please supply Cq set at least of length>=3"))
    }
    if (beta & is.null(K)) {
        warning(paste("Warning: if beta==TRUE, you should specify 'K' as a fixed parameter (or default 1.0)"))
    }
    if (is.null(EPCR)) {
        warning(paste("Warning: you should specify 'EPCR' as a fixed parameter"))
    }

    model.continuous <- FALSE
    if (any(is.na(N)) | is.null(N)) {
        model.continuous <- TRUE
        if (print.level>0) {
            cat("\nWarning: N is NULL or containing NA elements.",
                "The function runs with the continuous model.\n", sep=" ")
        }
    } else if (length(N) < length(target1)) {
        model.continuous <- TRUE
        if (print.level>0) {
            cat("\nWarning: N is shorter than length(target1).",
                "N is discarded and the continuous model, Beta(APK, A(1-P)K), is assumed.\n", sep=" ")
        }
    } else {
        if (print.level>0) {
            cat("\nThe model is discrete.\n")
        }
    }

    DCW <- target0-housek0
    DCD <- target1-housek1
    deldel <- DCD-DCW

    # a fixed-length vector that stores each parameter is target of estimation （TRUE） or fixed （FALSE）
    is.unknown <- c(    P=is.null(P), K=is.null(K), targetScale=is.null(targetScale),
                        sdMeasure=is.null(sdMeasure), EPCR=is.null(EPCR))
    # Full parameter vector is initialized with the default values. In linear scale.
    param0.full <- c(   P=1.99^-max(mean(deldel), 0.1), K=1.0, targetScale=exp(-mean(DCW)*log(1.99)),
                        sdMeasure=0.5, EPCR=0.99   )
    # Extract fixed parameters.
    para.fixed <- numeric(0)
    try({
        para.fixed <- c(P=P, K=K, targetScale=targetScale, sdMeasure=sdMeasure, EPCR=EPCR)
        # fixed part of the parameter vector is substituted with user-specified values.
        param0.full[names(is.unknown[!is.unknown])] <- para.fixed
    }, silent=TRUE) # try is used because it returns error when all parameters are unknown (substitution=NULL)
    # If the initial values for the unknown parameters are specified by user
    if (length(XInit0)>0) {
        nam.giv.init <- intersect(names(XInit0), names(is.unknown[is.unknown]))
        param0.full[nam.giv.init] <- XInit0[nam.giv.init]
    }

    # 'is.unknown' and 'para.fixed' stand for "which parameters are unknown" and "list of fixed parameters"
    # 'XInit' is input for the optimizer: only containning unknown parameters and scales are transformed.
    XInit <- log(param0.full[is.unknown]) #
    if (is.unknown["P"]==TRUE) {
        XInit["P"] <- qlogis(param0.full["P"]) # only 'P' is transformed in logit.
    }
    if(print.level>0) {
        cat("\nInitial values given by the user:\n")
        print(param0.full)
        cat("Which parameters are unknown?\n")
        print(is.unknown)
        cat("Fixed parameters:\n")
        print(para.fixed)
        cat("Initial value for optimization:\n")
        print(XInit)
        cat("\nGo ahead.\n")
    }

#    stop("Stop here.") # for debug

    # Optimization for the elements of XInit
    fscale <- 1 # minimization
    Param <- rep(NA, length(param0.full))
    setTimeLimit(elapsed=maxtime) # Functions to set elapsed time limits for the current session {}.
    ptime0 <- proc.time()
    success.nlm <- try({
        if (model.continuous==TRUE) {
            Z <- nlm(   f=.freqpcr_loglike_cont, p=XInit,
                        A=A, DCW=DCW, DCD=DCD, zeroAmount=zeroAmount, para.fixed=para.fixed,
                        beta=beta,
                        hessian=TRUE, fscale=fscale, print.level=print.level,
                        gradtol=gradtol, steptol=steptol, iterlim=iterlim   )
        } else {
            Z <- nlm(   f=.freqpcr_loglike, p=XInit,
                        N=N, DCW=DCW, DCD=DCD, zeroAmount=zeroAmount, para.fixed=para.fixed,
                        beta=beta, diploid=diploid,
                        hessian=TRUE, fscale=fscale, print.level=print.level,
                        gradtol=gradtol, steptol=steptol, iterlim=iterlim   )
        }
        Param[which(is.unknown)] <- Z$estimate
    })
    cal.time <- proc.time()-ptime0
    setTimeLimit(Inf)
    flush.console()

    # make the template of the result table
    result <- matrix(numeric(5*6), byrow=FALSE, ncol=6) # nrow == length(param0.full)
    rownames(result) <- c(  "P (R-allele frequency)",
                            "K (gamma shape parameter)",
                            "targetScale (relative amount of target locus)",
                            "Cq measurement error (SD)",
                            "EPCR (Duplication efficiency of PCR)")
    colnames(result) <- c(  "Estimate", "Fixed", "(scaled)", "(scaled.SE)",
                            paste(deparse(100*   pvalue/2) , "%", sep=""),
                            paste(deparse(100*(1-pvalue/2)), "%", sep="") )
    qvalue <- qnorm(p=pvalue/2, mean=0, sd=1, lower.tail=FALSE) # to be 2.5 percentile of N(0, 1) if pvalue=0.05

    # if calculation of Z failed:
    if (class(success.nlm) == "try-error") {
        cat("Calculation was terminated:\n")
        print(paste("Maximum elapsed time was set", maxtime, "seconds:", "elapsed", cal.time[3], "seconds.", sep=" "))
        cat("\n")
        if (print.level > 0) {
            title <- paste( "Maximum-likelihood estimates with the two-sided ",
                            deparse(100*(1-pvalue)), "% CIs", sep="")
            cat(title, "\n", sep=""); print(result); cat("\n");
        }
        return( CqFreq(report=result, obj=list(iterations=NA), cal.time=cal.time) )
    } else {
        if(print.level > 0) {
            cat("Optimization ends in\n")
            print(cal.time)
            cat("\n")
            print(Z)
        }
        Est <- exp(Param) # Optimization was conducted in log scale.
        Est[1] <- plogis(Param[1]) # Only P was estimated in logit scale.
        Est[which(!is.unknown)] <- para.fixed # Overwrite the fixed parameters

        # Standard errors are calculated from Hessian matrix: sometimes fail to have the inverse
        SE <- rep(NA, length(param0.full))
        SEraw <- numeric(length(XInit))
        success.solve <- try({
            SEraw <- sqrt(diag(solve(fscale*Z$hessian))) # When minimized, diag() will be positive.
            SE[which(is.unknown)] <- SEraw
        }, silent=TRUE)

        result[, 1] <- Est
        result[, 2] <- !is.unknown
        result[, 3] <- Param
        result[, 4] <- SE
        result[, 5] <- Param-qvalue*SE
        result[, 6] <- Param+qvalue*SE
        result[1, 5:6] <- plogis(result[1, 5:6])
        result[2:nrow(result), 5:6] <- exp(result[2:nrow(result), 5:6])

        if (print.level > 0) {
            title <- paste( "Maximum-likelihood estimates with the two-sided ",
                            deparse(100*(1-pvalue)), "% CIs", sep="")
            cat(title, "\n", sep=""); print(result); cat("\n");
        }
        return( CqFreq(report=result, obj=Z, cal.time=cal.time) )
    }
}


#' @title Log-likelihood of obtaining Cq values under given parameter set.
#'
#' @description The internal function is called from the optimizer (\code{\link[stats]{nlm}()}) running in \code{\link{freqpcr}()}. It defines the log-likelihood by obtaining the two \eqn{\Delta}Cq values (differences in the four Cq measurements) provided that the allele mixing ratio for each bulk sample is given together with other parameters. This function was vectorized over multiple bulk samples.
#' @param X Numeric vector that stores the parameter values to be optimized via \code{\link{nlm}()}: \code{P} in logit scale and \code{K}, \code{targetScale}, \code{sdMeasure}, and \code{EPCR} in log scale.
#' @param N Sample sizes as a numeric vector. \code{N[i]} signifies the number of individuals (both for haploidy and diploidy) contained in the \emph{i}th bulk sample.
#' @param DCW,DCD Numeric vectors having the same length as \code{N}. They store the measured values of the two \eqn{\Delta}Cq: \code{DCW (= target0 - housek0)} and \code{DCD (= target1 - housek1)}. They can contain NA (simply ignored in the calculation).
#' @param para.fixed Named numeric vector that stores the fixed parameters inherited from \code{\link{freqpcr}()}, if specified. By default (\code{NULL}), all the parameters (\code{P}, \code{K}, \code{targetScale}, \code{sdMeasure}, and \code{EPCR}) are unknown. Unlike \code{X}, each element value is set in linear scale.
#' @inheritParams freqpcr
#' @return Scalar of the log likelihood.
.freqpcr_loglike <- function(X, N, DCW, DCD, zeroAmount, para.fixed=NULL, beta=TRUE, diploid=FALSE) {
    if (is.null(para.fixed)) {
        para.fixed <- numeric(0)
    }
    names(X) <- setdiff(c("P", "K", "targetScale", "sdMeasure", "EPCR"), names(para.fixed))
    P           <- ifelse(is.na(X["P"]),            para.fixed["P"],            plogis(X["P"])          )
    K           <- ifelse(is.na(X["K"]),            para.fixed["K"],            exp(X["K"])             )
    targetScale <- ifelse(is.na(X["targetScale"]),  para.fixed["targetScale"],  exp(X["targetScale"])   )
    sdMeasure   <- ifelse(is.na(X["sdMeasure"]),    para.fixed["sdMeasure"],    exp(X["sdMeasure"])     )
    EPCR        <- ifelse(is.na(X["EPCR"]),         para.fixed["EPCR"],         exp(X["EPCR"])     )

    ploidy <- ifelse(diploid, 2, 1) # a diploid individual is treated as "two independent haploids" (tricky!)
    xsm <- 2.0 # xsm is the accumulation of normal error. 2 for the delta Cq values (DCW, DCD).
#    cat("¥n")
#    print(c(P, K, targetScale, sdMeasure, EPCR=EPCR, ploidy=ploidy))

    ans <- c(LogLike=0.0) # Define log likelihood. From Eq. 11, DCW is not affected by m.
    LL.DCW <- dnorm(DCW, mean=-log(targetScale)/log(1.0+EPCR), sd=sqrt(xsm)*sdMeasure, log=TRUE)
    LL.DCW <- LL.DCW[!is.na(DCW)] # if DCW contained NA values, these results should be removed.

    # DCD is affected by (zXS+XR)/(XS+XR)=z+(1-z)y, where y=XR/(XS+XR)

    # Case1: 0 out of N[h] are R.
    ans0 <- dbinom(0, size=ploidy*N, prob=P, log=FALSE) *
            dnorm(DCD, mean=-log(zeroAmount*targetScale)/log(1.0+EPCR), sd=sqrt(xsm)*sdMeasure, log=FALSE)

    # Case3: N[h] out of N[h] are R. 2*N in the case of diploidy.
    ansN <- dbinom(ploidy*N, size=ploidy*N, prob=P, log=FALSE) *
            dnorm(DCD, mean=-log(targetScale)/log(1.0+EPCR), sd=sqrt(xsm)*sdMeasure, log=FALSE)

    # Case2: m out of N[h] are R. Only considered when N[h] >= 2.
    # When N[h] == 1, m only takes 0 and 1. These cases are coverd with ans0 and ansN.
    ansm <- N*0
    if (diploid==FALSE) {
        # Case of haploidy
        for (h in 1:length(N)) {
            if (N[h]>1) {
                m <- sequence(N[h]-1) # m:= 1, 2, 3, ..., N[h]-1
                n <- rep(N[h], times=N[h]-1) # n:= N[h], N[h], ...
#                print(rbind(m, n, rep(DCD[h], times=N[h]-1)))
                if (beta==TRUE) {
                    ansm.long <-dbinom(m, size=n, prob=P, log=FALSE) *
                                .integrate_beta(rep(DCD[h], times=N[h]-1), K*m, K*(n-m),
                                                zeroAmount, targetScale, sdMeasure, xsm, EPCR,
                                                cubmethod="hcubature", relTol=1e-1, absTol=1e-8, maxEval=10^6)
                } else {
                    ansm.long <-dbinom(m, size=n, prob=P, log=FALSE) *
                                .integrate_gamma(rep(DCD[h], times=N[h]-1), K*m, K*(n-m),
                                                zeroAmount, targetScale, sdMeasure, xsm, EPCR,
                                                cubmethod="hcubature", relTol=1e-1, absTol=1e-8, maxEval=10^6)
                }
#                print(ansm.long)
                ansm[h] <- sum(ansm.long)
            }
        }
    } else {
        # Case of diploidy
        # Assuming i.i.d. between the amounts of R and S chromosomes owned by heterozygotes.
        # We can use \code{\link{.integrate_gamma}()} or \code{\link{.integrate_beta}()} with the assumption:
        # shapeR=K*(colSums(combi)+m1-m0), shapeS=K*(colSums(combi)-m1+m0), obs=rep(DCD[h], times=length(m1))
        for (h in 1:length(N)) {
            if (N[h]>1) {
                m1a <- rep(c(0, sequence(N[h]-1)), times=N[h]) # c(0, sequence(N[h]-1): 0, 1, 2, 3, 4, 5 when N[h]=6
                m0a <- rep(c(0, sequence(N[h]-1)), each=N[h])
                nha <- m1a+m0a # "mha" may exceed the sample size, N[h]; such the cases should be removed.
                m1 <- m1a[nha<=N[h]] # number of RR individuals (not the chromosome copies) in the bulk sample
                m0 <- m0a[nha<=N[h]] # number of SS individuals
                mh <-N[h]-m1-m0 # number of hetero individuals
                combi <- rbind(m1, mh, m0) # RR, SR, SS
#                print(combi)
                if (beta==TRUE) {
                    ansm.long <-apply(combi, 2, function(x) dmultinom(x, prob=c(P^2, 2*P*(1-P), (1-P)^2))) *
                                .integrate_beta(rep(DCD[h], times=length(m1)), K*(N[h]+m1-m0), K*(N[h]-m1+m0),
                                                zeroAmount, targetScale, sdMeasure, xsm, EPCR,
                                                cubmethod="hcubature", relTol=1e-1, absTol=1e-8, maxEval=10^6)
                } else {
                    ansm.long <-apply(combi, 2, function(x) dmultinom(x, prob=c(P^2, 2*P*(1-P), (1-P)^2))) *
                                .integrate_gamma(rep(DCD[h], times=length(m1)), K*(N[h]+m1-m0), K*(N[h]-m1+m0),
                                                zeroAmount, targetScale, sdMeasure, xsm, EPCR,
                                                cubmethod="hcubature", relTol=1e-1, absTol=1e-8, maxEval=10^6)
                }
                ansm[h] <- sum(ansm.long)
            } else {
                # when N[h] == 1 and heterozygote: same DNA quantities for xR and xS (i.i.d is not assumed).
                # DCD[h] when xR==xS -> (1+z)*X/2X = (1+zeroAmount)/2
                ansm[h] <-  dmultinom(c(0, 1, 0), prob=c(P^2, 2*P*(1-P), (1-P)^2)) *
                            dnorm(  DCD[h], mean=-(log(targetScale)+log((1+zeroAmount)/2))/log(1.0+EPCR),
                                    sd=sqrt(xsm)*sdMeasure  )
            }
        }
    }
#    print(rbind(ans0, ansm, ansN))
    flush.console()
    LL.DCD <- log(ans0+ansm+ansN)
    LL.DCD <- LL.DCD[!is.na(DCD)] # if DCD contained NA values, these results should be removed.
    return(-sum(c(LL.DCW, LL.DCD)))
}


#' @title Log-likelihood when sample allele ratio is continuous.
#'
#' @description Called from \code{\link{freqpcr}()} instead of \code{\link{.freqpcr_loglike}()} when the model is `continuous'. This function assumes that each sample does not consists of \code{n} individual organisms with certain genotypes, but the result of a direct sampling from the continuous population having the allele ratio around \code{p:(1-p)}. Each sample allele ratio is considered to follow \code{Beta(pk, (1-p)k)}.
#' @param X Numeric vector that stores the parameter values to be optimized via \code{\link{nlm}()}: \code{P} in logit scale and \code{K}, \code{targetScale}, \code{sdMeasure}, and \code{EPCR} in log scale.
#' @param A Relative DNA content between the samples in the continuous model, as a numeric vector.
#' @param DCW,DCD Numeric vectors. They store the measured values of the two \eqn{\Delta}Cq: \code{DCW (= target0 - housek0)} and \code{DCD (= target1 - housek1)}.
#' @param para.fixed Named numeric vector that stores the fixed parameters inherited from \code{\link{freqpcr}()}, if specified. By default (\code{NULL}), all the parameters (\code{P}, \code{K}, \code{targetScale}, \code{sdMeasure}, and \code{EPCR}) are unknown. Unlike \code{X}, each element value is set in linear scale.
#' @inheritParams freqpcr
#' @return Scalar of the log likelihood.
.freqpcr_loglike_cont <- function(X, A, DCW, DCD, zeroAmount, para.fixed=NULL, beta=TRUE) {
    if (is.null(para.fixed)) {
        para.fixed <- numeric(0)
    }
    names(X) <- setdiff(c("P", "K", "targetScale", "sdMeasure", "EPCR"), names(para.fixed))
    P           <- ifelse(is.na(X["P"]),            para.fixed["P"],            plogis(X["P"])          )
    K           <- ifelse(is.na(X["K"]),            para.fixed["K"],            exp(X["K"])             )
    targetScale <- ifelse(is.na(X["targetScale"]),  para.fixed["targetScale"],  exp(X["targetScale"])   )
    sdMeasure   <- ifelse(is.na(X["sdMeasure"]),    para.fixed["sdMeasure"],    exp(X["sdMeasure"])     )
    EPCR        <- ifelse(is.na(X["EPCR"]),         para.fixed["EPCR"],         exp(X["EPCR"])     )

    xsm <- 2.0 # xsm is the accumulation of normal error. 2 for the delta Cq values (DCW, DCD).

    ans <- c(LogLike=0.0)
    LL.DCW <- dnorm(DCW, mean=-log(targetScale)/log(1.0+EPCR), sd=sqrt(xsm)*sdMeasure, log=TRUE)
    LL.DCW <- LL.DCW[!is.na(DCW)]

    # assuming X_R ~ Gamma(pk, theta), x_S ~ Gamma((1-p)k, theta),
    # and y_R ~ Beta(shape1=pk, shape2=(1-p)k) under the regional allele frequency p.

    ansm <- DCD*0
    for (h in 1:length(DCD)) {
        if (beta==TRUE) {
            ansm[h] <-.integrate_beta(  DCD[h], A*K*P, A*K*(1-P),
                                        zeroAmount, targetScale, sdMeasure, xsm, EPCR,
                                        cubmethod="hcubature", relTol=1e-1, absTol=1e-8, maxEval=10^6  )
        } else {
            ansm[h] <-.integrate_gamma( DCD[h], A*K*P, A*K*(1-P),
                                        zeroAmount, targetScale, sdMeasure, xsm, EPCR,
                                        cubmethod="hcubature", relTol=1e-1, absTol=1e-8, maxEval=10^6 )
        }
    }
#    print(ansm)
    flush.console()
    LL.DCD <- log(ansm)
    LL.DCD <- LL.DCD[!is.na(DCD)] # if DCD contained NA values, these results should be removed.
    return(-sum(c(LL.DCW, LL.DCD)))
}
