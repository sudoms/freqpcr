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

#' @description The function estimates the population allele frequency using the dataset of Cq values measured over \code{length(N)} bulk samples, each of which has the sample size of \code{N[i]} as the number of individuals included. \code{N[i]} can be 1, meaning every individual is analyzed separately. For the \emph{i}th (bulk) sample, the four Cq values have been measured as \code{target0[i]}, \code{target1[i]}, \code{housek0[i]}, and \code{housek1[i]}, respectively.
#' \cr
#' The function can estimate up to five parameters simultaneously when the Cq measures are available for more than two bulk samples: \code{P}, \code{K}, \code{targetScale}, \code{sdMeasure}, and \code{EPCR}, though the last one is not reliable in many cases and users are encouraged to know the size for \code{EPCR} beforehand by using \code{\link{knownqpcr}()} or \code{\link{knownqpcr_unpaired}()}. By default, parameters other than \code{EPCR} are estimated.
#' @param N Sample sizes as a numeric vector. \code{N[i]} signifies the number of individuals (both for haploidy and diploidy) contained in the \emph{i}th bulk sample.
#' @param housek0 A numeric vector having the same length as \code{N}. In RED-\eqn{\Delta\Delta}Cq method, \code{housek0} is the Cq values of the test sample without the restriction enzyme digestion, which is amplified with the primer set for a housekeeping gene. In a general \eqn{\Delta\Delta}Cq analyses, \code{housek0} is defined for the control sample (100\% mutant) solution, which is also amplified with the primer set for the housekeeping gene.
#' @param target0 A numeric vector having the same length as \code{N}. In RED-\eqn{\Delta\Delta}Cq method, \code{target0[i]} signifies the measured Cq value of the \emph{i}th bulk sample without the digestion, which is treated with the primer set that amplifies both alleles, wild-type (or S: susceptible) and mutant (or R: resistant to a pesticide), on the target locus. In general \eqn{\Delta\Delta}Cq analyses, \code{target0} is the Cq values of the control sample, which is defined as a sample solution that contains 100\% R, amplified with an R-specific primer set.
#' @param housek1 A numeric vector having the same length as \code{N}. The Cq values of the test sample measured on the housekeeping gene after the restriction enzyme digestion (in RED-\eqn{\Delta\Delta}Cq method), or the test sample amplified on the housekeeping gene (in general \eqn{\Delta\Delta}Cq analyses).
#' @param target1 A numeric vector having the same length as \code{N}. For each test sample with unknown allele-ratio, \code{target1[i]} is defined as the Cq value for the target locus amplified with the universal primer set after the restriction enzyme digestion (in RED-\eqn{\Delta\Delta}Cq method), or the target locus amplified with the R-allele-specific primer set (in general \eqn{\Delta\Delta}Cq analyses).
#' @param P Population allele frequency from which the test samples are retrieved. Default is \code{NULL} and to be estimated. If the parameter is known, it should be given as a numeric between 0 and 1.
#' @param K The gamma shape parameter of the individual DNA yield. Default is \code{NULL} and to be estimated. If known, given as a positive numeric.
#' @param targetScale (\eqn{\delta_{T}}) The relative template DNA amount of the target locus to the houskeeping locus. If known, given as a positive numeric.
#' @param sdMeasure (\eqn{\sigma_{c}}) The measurement error on each Cq value following Normal(0, \eqn{\sigma_{c}^2}). If known, given as a positive numeric.
#' @param EPCR (\eqn{\eta}) Amplification efficiency per PCR cycle. If known, given as a positive numeric. When \code{EPCR = 1}, template DNA doubles every cycle (\code{EPCR + 1 = 2}).
#' @param XInit0 Optionally the initial value for the parameter optimization can be specified, but the user is highly recommended to keep the argument as is. Unlike XInit in \code{\link{knownqpcr}()}, the argument is not directly fed to the optimizer. Each value is used only when the parameter is unknown.
#' @param zeroAmount (In RED-\eqn{\Delta\Delta}Cq method) residual rate of restriction enzyme digestion, or (in general \eqn{\Delta\Delta}Cq analyses) small portion of the off-target allele amplified in the PCR process. It needs to be always specified by the user as a number between 0 and 1, usually near 0.
#' @param beta Whether to use the beta distribution for the approximation of allele frequency instead of the two gamma distributions for the DNA quantities of the two alleles? Default is \code{TRUE}, which considerably accelerates the calculation.
#' @param diploid Is the target organism diploidy? Default is \code{FALSE}, assuming haploidy. Current implementation on the diploids assumes i.i.d. between the amounts of R and S chromosomes owned by a heterozygote, which is unlikely in many animals but necessary for the calculation in a realistic time.
#' @param pvalue The two-sided confidence interval is calculated at the last iteration based on the significance level. Default is 0.05, which returns the 95\% CI (2.5 to 97.5 percentile) based on the Hessian matrix.
#' @param gradtol,steptol,iterlim Arguments passed to \code{\link[stats]{nlm}()}. \code{gradtol} and \code{steptol} are the positive scalars giving the tolerance to terminate the algorithm and the minimum allowable relative step length. \code{iterlim} specifies the maximum number of iterations to be performed before the program is terminated. Usually 30 iterations are enough.
#' @param maxtime a positive scalar to set the maximum calculation time in seconds to abort the optimizer (and return error). The total calculation time largely depends on \code{N[i]}, the number of individuals contained in each bulk sample.
#' @param print.level \code{print.level=1} (the default) shows the initial values of the parameters and likelihood as well as the output in the last iteration. \code{print.level=2} shows the parameter values and gradients in every step. \code{print.level=0} does not output any intermediate state to R console, simply returning the result summary.
#' @param ... Additional arguments passed to the function.
#' @return Object of the S4 class \linkS4class{CqFreq}. The slot \code{report} is a matrix and each row contains the estimated parameter value with 100*(1-pvalue)\% confidence interval. The following parameters are returned:
#' \enumerate{
#'     \item\code{P} Population allele frequency from which the test samples are retrieved.
#'     \item\code{K} The gamma shape parameter of the individual DNA yield.
#'     \item\code{targetScale} (\eqn{\delta_{T}}), the relative template DNA amount of the target to the houskeeping loci.
#'     \item\code{EPCR} (\eqn{\eta}), the amplification efficiency per PCR cycle.
#'     \item\code{sdMeasure} or "Cq measurement error" (\eqn{\sigma_{c}}).
#' }
#' @section Choise of the parameters to be estimated:
#' Estimation is conducted only for parameters specified, explicitly or implicitly, as \code{NULL}. If one feeds a value e.g. \code{K=1} or \code{sdMeasure=0.24}, it is then treated as fixed parameter.
#' \cr
#' \code{EPCR} is the amplification efficiency of PCR and usually not estimable from the experiments with unknown allele ratios. You should verify the size of \code{EPCR} and \code{zeroAmount} (the residue rate of the restriction enzyme digestion in the RED-\eqn{\Delta\Delta}Cq method) prior to the experiment with unknown samples. The functions \code{\link{knownqpcr}()} and \code{\link{knownqpcr_unpaired}()} provide the procedure to estimate the sizes of the experimental parameters using the DNA solutions of known allele mixing ratios.
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
#' result <- freqpcr( dmy_cq@N,
#'                    dmy_cq@housek0, dmy_cq@target0, dmy_cq@housek1, dmy_cq@target1,
#'                    EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=2 )
#' print(result)
#' # Estimation with freqpcr, where you have knowledge on the size of baseChange.
#' result <- freqpcr( dmy_cq@N,
#'                    dmy_cq@housek0, dmy_cq@target0, dmy_cq@housek1, dmy_cq@target1,
#'                    baseChange=0.3,
#'                    EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=1 )
#' @export
#' @family estimation procedures
freqpcr <- function(N, housek0, target0, housek1, target1,
                    P=NULL, K=NULL, targetScale=NULL, sdMeasure=NULL, EPCR=0.99,
                    XInit0=c(P=NULL, K=NULL, targetScale=NULL, sdMeasure=NULL, EPCR=NULL),
                    zeroAmount=0.0016, beta=TRUE, diploid=FALSE,
                    pvalue=0.05, gradtol=1e-4, steptol=1e-9, iterlim=100, maxtime=600, print.level=1, ...) {
    if (length(N)<1) {
        stop(paste("Error: please supply N of length>=1"))
    }
    if (is.null(EPCR) & length(N)<3) {
        stop(paste("Error: when EPCR is unknown, please supply N at least of length>=3"))
    }
    if (beta & is.null(K)) {
        warning(paste("Warning: if beta==TRUE, you should specify 'K' as a fixed parameter (or default 1.0)"))
    }
    if (is.null(EPCR)) {
        warning(paste("Warning: you should specify 'EPCR' as a fixed parameter"))
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
        cat("\nGiven initial values by the user:\n")
        print(param0.full)
        cat("Parameters are unknown?\n")
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
        Z <- nlm(   f=.freqpcr_loglike, p=XInit,
                    N, DCW, DCD, zeroAmount=zeroAmount, para.fixed=para.fixed,
                    beta=beta, diploid=diploid,
                    hessian=TRUE, fscale=fscale, print.level=print.level,
                    gradtol=gradtol, steptol=steptol, iterlim=iterlim )
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


#' @title Log-likelihood of getting Cq values under given parameter set.
#'
#' @description Internal function called from the optimizer (\code{\link[stats]{nlm}()}) running in \code{\link{freqpcr}()}. It defines the log-likelihood getting the two \eqn{\Delta}Cq values (differences in the four Cq measurements) provided the allele mixing ratio for each bulk sample is given together with other parameters. The function is vectorized over multiple (bulk) samples.
#' @param X Numeric vector that stores the parameter values to be optimized via \code{\link{nlm}()}, \code{P} in logit scale and \code{K}, \code{targetScale}, \code{sdMeasure}, and \code{EPCR} in log scale.
#' @param DCW,DCD Numeric vectors having the same length as \code{N}. They store the measured values of the two \eqn{\Delta}Cq: \code{DCW (= target0 - housek0)} and \code{DCD (= target1 - housek1)}.
#' @param para.fixed Named numeric vector that stores the fixed parameters in \code{\link{freqpcr}()}, if specified by the user. By default (\code{NULL}), all the parameters are considered unknown. Unlike \code{X}, each element should be the raw parameter value in linear scale.
#' @inheritParams freqpcr
#' @return A scalar of the log likelihood.
.freqpcr_loglike <- function(  X, N, DCW, DCD, zeroAmount, para.fixed=NULL,
                                beta=TRUE, diploid=FALSE  ) {
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
    return(-sum(LL.DCW+log(ans0+ansm+ansN)))
}
