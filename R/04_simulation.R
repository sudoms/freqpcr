# 04 simulation

#' S4 class that stores the dummy Cq data for performance test.
#'
#' A dummy Cq dataset suitable for the test of the package, typically the output of \code{\link{makeCqList}()}.
#' @slot N Sample sizes as a numeric vector. \code{N[i]} signifies the number of individuals (both for haploidy and diploidy) contained in the \emph{i}th bulk sample.
#' @slot m As for haploidy, \code{m} is a matrix with 2 rows and \code{ntrap} columns. It is the matrix with 3 rows and \code{ntrap} columns for diploidy, for which each column means the segregation ratio of the R (mutant) and S (wild type) individuals in each bulk sample. \code{m[1, i]} stores the number of R (or R homozygous) individuals while \code{N[i] = sum(m[, i])} specifies the total number in the bulk sample. As for haploidy, \code{m[2, i]} stands for the number of S individuals. \code{m[2, i]} and \code{m[3, i]} in diploidy stand for the number of heterozygotes and S homozygotes, respectively.
#' @slot xR,xS Numeric vector of the same length with N. \code{xR[i]} stores the amount of the template DNA for the R allele in the \emph{i}th bulk sample.
#' @slot housek0,target0,housek1,target1 Numeric vectors of the same lengths with N. Store the Cq values.
#' @slot DCW \eqn{\Delta}Cq value measured on the bulk sample without endonuclease digestion, \code{DCW} is defined as (\code{target0-housek0}).
#' @slot DCD \eqn{\Delta}Cq value measured on the bulk sample after endonuclease digestion, \code{DCD} is defined as (\code{target1-housek1}).
#' @slot deldel \eqn{\Delta\Delta}Cq value, \code{deldel} is defined as (\code{DCD - DCW}).
#' @slot RFreqMeasure A classical index of the allele frequency calculated for each bulk sample, which is defined as \code{(1.0+EPCR)^(-deldel)}.
#' @slot ObsP Defined as \code{min(RFreqMeasure, 1)} because \code{RFreqMeasure} can exceed 1 by definition.
#' @slot rand.seed The seed of the random-number generator (RNG) fed to the current R session to generate dummy \code{m}, \code{xR} and \code{xS} data. It is usually specified by user.
#' @export
CqList <- setClass("CqList",
    slots = c(
        N="numeric",
        m="matrix",
        xR="numeric", xS="numeric",
        housek0="numeric", target0="numeric", housek1="numeric", target1="numeric",
        DCW="numeric", DCD="numeric",
        deldel="numeric",
        RFreqMeasure="numeric",
        ObsP="numeric",
        rand.seed="numeric"
    )
)




#' Generate dummy DNA dataset ready for allele-frequency estimation.
#'
#' The function generates a dummy dataset of the Cq values. You can directly feed the output of this function to the first argument of \code{\link{simCqList}()}.
#' @param rand.seed Seed for the R built-in random-number-generator used in the data generation.
#' @param P Population allele frequency from which the test samples are generated. It is given as a numeric between 0 and 1.
#' @param K The gamma shape parameter of the individual DNA yield, given as a positive numeric.
#' @param ntrap,npertrap Scalar specifying the number of bulk samples (\code{ntrap}) or the numbers of individuals (both for haploidy and diploidy) contained in each bulk sample (\code{npertrap}). Currently limited to the case all bulk samples have the same sample size e.g. (4+4+4) when \code{ntrap = 3} and \code{npertrap = 4} hold.
#' @param scaleDNA Small positive scalar that specifies the scale parameter of the gamma distribution appriximating the yield of DNA from (per-haploid) individual. The yield of \code{2*scaleDNA} is expected from a diploid. The quantity is determined as the relative amount (in linear scale) to the termination threshold of the realtime PCR cycle.
#' @param targetScale (\eqn{\delta_{T}}) The relative template DNA amount of the target locus to the houskeeping locus, given as a positive numeric.
#' @param baseChange (\eqn{\delta_{B}}) The change rate in the template DNA quantities after the restriction enzyme digestion (in the RED-delta delta Cq method), given as a positive numeric. This parameter is not used in \code{\link{freqpcr}()}, but essential for \code{\link{knownqpcr}()} and \code{\link{makeCqList}()}.
#' @param EPCR (\eqn{\eta}) Amplification efficiency per PCR cycle, given as a positive numeric. When \code{EPCR = 1}, template DNA doubles every cycle (\code{EPCR + 1 = 2}).
#' @param sdMeasure (\eqn{\sigma_{c}}) The measurement error on each Cq value following Normal(0, \eqn{\sigma_{c}^2}), given as a positive number.
#' @inheritParams freqpcr
#' @return Object of the S4 class \linkS4class{CqList}, storing the dummy experiment data of Cq-based qPCR analysis. Note that a \linkS4class{CqList} object in no way contains original information on \code{P}, \code{K}, \code{targetScale}, \code{sdMeasure}, and \code{EPCR}.
#' @examples
#' P <- 0.25
#' # Calculate the segregation ratios for six bulk samples (1000 individuals for each).
#' rbinom(n=6, size=1000, prob=P)
#' rmultinom(n=6, size=1000, prob=c(P, 1-P))
#' rmultinom(6, size=1000, prob=c(P^2, 2*P*(1-P), (1-P)^2))
#'
#' # Make a dummy Cq dataset of six bulk DNA samples,
#' # each of which comprises eight haploid individuals.
#' dmy_cq <- makeCqList( rand.seed=1, P=0.75, K=2, ntrap=6, npertrap=8, scaleDNA=1e-07,
#'                       targetScale=1.5, baseChange=0.3, EPCR=0.95, zeroAmount=1e-3,
#'                       sdMeasure=0.3, diploid=FALSE )
#' print(dmy_cq)
#' @export
makeCqList <- function( rand.seed, P, K, ntrap, npertrap, scaleDNA=(1/K)*1e-06,
                        targetScale, baseChange, EPCR, zeroAmount, sdMeasure, diploid=FALSE ) {
    set.seed(rand.seed)
    N <- rep(npertrap, ntrap)
    if (!diploid) {
        m <- rmultinom(n=ntrap, size=npertrap, prob=c(P, 1-P)) # R and S in each trap
        xR <- rgamma(n=ntrap, shape=K*m[1, ], scale=scaleDNA) # first row of m: R
        xS <- rgamma(n=ntrap, shape=K*m[2, ], scale=scaleDNA) # second row of m: S
    } else {
        m <- rmultinom(n=ntrap, size=npertrap, prob=c(P^2, 2*P*(1-P), (1-P)^2)) # RR, RS, and SS in each trap
        mx <- matrix(rgamma(n=length(c(m)), shape=c(m)*K, scale=scaleDNA), byrow=FALSE, nrow=nrow(m))
        xS <- 2*mx[3, ]+mx[2, ]
        xR <- 2*mx[1, ]+mx[2, ]
    }
    housek0 <- rnorm(n=ntrap, mean=-log(xR+xS)/log(1.0+EPCR), sd=sdMeasure)
    target0 <- rnorm(n=ntrap, mean=-(log(targetScale) + log(xR+xS))/log(1.0+EPCR), sd=sdMeasure)
    housek1 <- rnorm(n=ntrap, mean=-(log(baseChange) + log(xR+xS))/log(1.0+EPCR), sd=sdMeasure)
    target1 <- rnorm(n=ntrap, mean=-(log(targetScale) + log(baseChange) + log(xR+zeroAmount*xS))/log(1.0+EPCR),
                        sd=sdMeasure)
    set.seed(NULL)
    DCW <- target0-housek0
    DCD <- target1-housek1
    deldel <- DCD - DCW
    RFreqMeasure <- (1.0+EPCR)^(-deldel)
    ObsP <- ifelse(RFreqMeasure>1, 1, RFreqMeasure)
    return( CqList( N=N, m=m,
                    xR=xR, xS=xS, housek0=housek0, target0=target0, housek1=housek1, target1=target1,
                    DCW=DCW, DCD=DCD, deldel=deldel,
                    RFreqMeasure=RFreqMeasure, ObsP=ObsP, rand.seed=rand.seed ) )
}



#P <- 0.25
#K <- 100
#scaleDNA <- 10^-2
#(m <- rmultinom(6, size=100000, prob=c(P^2, 2*P*(1-P), (1-P)^2)))
#(mx <- matrix(rgamma(n=length(c(m)), shape=c(m)*K, scale=scaleDNA), byrow=FALSE, nrow=nrow(m)))




#' Simulate parameter estimation based on user-generated dummy data.
#'
#' Wrapper of \code{\link{freqpcr}()} function suitable for the test of the function using randomly generated datasets i.e., the output of \code{\link{makeCqList}()}.
#' @param CqList Object belonging to the \linkS4class{CqList} class, typically the output from \code{\link{makeCqList}()}. Having the slots \code{N}, \code{target0}, \code{target1}, \code{housek0}, and \code{housek1}, all of which are numeric vectors of the same length.
#' @param P,K,targetScale,sdMeasure If NULL, the parameter is considered unknown and estimated internally via \code{\link{freqpcr}()}. If the value is specified, it is passed to \code{\link{freqpcr}()} as a fixed parameter. \code{EPCR} and \code{zeroAmount} are always treated as fixed parameters, for which values must be supplied.
#' @param beta,diploid,maxtime,print.level Parameters which are passed directly to \code{\link{freqpcr}()}.
#' @param aux Additional information to be displayed. The default is \code{NULL}. If user input some values, it is echoed to stdout together with the contents of the argument \code{CqList}. This option is convenient when you want to record the original dummy dataset and the corresponding result sequentially e.g. using \code{capture.output()}.
#' @inheritParams freqpcr
#' @return Object of the S4 class \linkS4class{CqFreq}.
#' @examples
#' K <- 2 # You already know the size of K in this case.
#' EPCR <- 0.97 # The sizes of EPCR and zeroAmount must always be supplied.
#' zeroAmount <- 1.6e-03
#' is.diploid <- FALSE
#'
#' # First, make a dummy Cq dataset of six bulk DNA samples,
#' # each of which comprises eight haploid individuals.
#' dmy_cq <- makeCqList( rand.seed=1, P=0.75, K=K, ntrap=6, npertrap=8, scaleDNA=1e-07,
#'                       targetScale=1.5, baseChange=0.3, EPCR=EPCR,
#'                       zeroAmount=zeroAmount, sdMeasure=0.3, diploid=is.diploid )
#'
#' # Then, estimate the population allele frequency from the dummy Cq dataset.
#' simCqList( CqList=dmy_cq, EPCR=EPCR, zeroAmount=zeroAmount, K=K,
#'             beta=TRUE, diploid=is.diploid, maxtime=60, print.level=2, aux="test" )
#'
#' # If the maximum calculation time was too short to converge, nlm() returns error.
#' # simCqList() avoids error and returns a matrix filled with zeros.
#' simCqList( CqList=dmy_cq, EPCR=EPCR, zeroAmount=zeroAmount, K=K,
#'             beta=FALSE, diploid=is.diploid, maxtime=0.01, print.level=2 )
#' @export
#' @family estimation procedures
simCqList <- function( CqList, EPCR, zeroAmount,
                        P=NULL, K=NULL, targetScale=NULL, sdMeasure=NULL,
                        beta, diploid, maxtime, print.level, aux=NULL, ... ) {
    if (!is.null(aux)) {
        cat("\n\n")
        cat("--------------------------------------------------------------------------------")
        cat("\n")
        print(c(PID=Sys.getpid()), quote=FALSE)
        print(aux, quote=FALSE)
        cat("\n")
        print(CqList)
        cat("\n")
        print(Sys.time(), quote=FALSE)
    }
    e1 <- try( {
        result.list <- freqpcr( N=CqList@N, # number of individuals (not chromosome sets)
                                target0=CqList@target0, target1=CqList@target1,
                                housek0=CqList@housek0, housek1=CqList@housek1,
                                P=P, K=K, targetScale=targetScale, sdMeasure=sdMeasure, EPCR=EPCR,
                                zeroAmount=zeroAmount, beta=beta, diploid=diploid,
                                # pvalue=0.05, gradtol=1e-4, steptol=1e-7, iterlim=100,
                                maxtime=maxtime,
                                print.level=print.level, ... )
    } )
    if (class(e1)=="try-error") {
        pvalue <- 0.05
        result <- matrix(numeric(5*6), byrow=FALSE, ncol=6) # nrow == length(param0.full)
        rownames(result) <- c(  "P (R-allele frequency)",
                                "K (gamma shape parameter)",
                                "targetScale (relative amount of target locus)",
                                "Cq measurement error (SD)",
                                "EPCR (Duplication efficiency of PCR)")
        colnames(result) <- c(  "Estimate", "Fixed", "(scaled)", "(scaled.SE)",
                                paste(deparse(100*   pvalue/2) , "%", sep=""),
                                paste(deparse(100*(1-pvalue/2)), "%", sep="") )
        ptime0 <- proc.time()
        cal.time <- proc.time()-ptime0
        result.list <- CqFreq(report=result, obj=list(iterations=NA), cal.time=cal.time)
    }
    print(Sys.time(), quote=FALSE)
    rm(e1)
    gc(); gc();
    return(result.list)
}