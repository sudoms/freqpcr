# 01 integrator

#' @title Integration of likelihood function based on the beta assumption.

#' @description Internal function to integrate the likelihood getting the \eqn{\Delta}Cq value (the argument \code{del}) over the entire range of the allele ratio (0 to 1). Vectorized for multiple bulk samples. Having the same arguments with .integrate_gamma.
#' @param del Numeric vector of the observed \eqn{\Delta}Cq values.
#' @param SHR The gamma shape parameters for the mutant (R) portion of the bulk samples. Should be the same vector length as del. Each element of SHR is defined as K*(the assumed number of R allele in the bulk sample: 1, 2, 3, ..., n-1).
#' @param SHS The gamma shape parameters for the wild (S) portion of the bulk samples. Should be the same length as del. Each element of SHS is defined as K*(the assumed number of S allele in the bulk sample).
#' @param xsm Specify the accumulation of the standard deviation of the Cq measuring errors when the *-Cq values are fed as difference. For \eqn{\Delta}Cq values, sdMeasure times two. For \eqn{\Delta\Delta}Cq, sdMeasure times four. Default is two and used in most cases.
#' @param cubmethod Cubature method passed to the integrator function. See the section "Methods for cubintegrate".
#' @param relTol The maximum tolerance passed to the cubature method. Though the default of cubature::cubintegrate function is 1e-5, the accuracy is reduced here to acceralate the integration.
#' @param absTol The absolute tolerance passed to the cubature method. The default is 1e-8, which is less accurate than the default of cubintegrate function (1e-12) but considered enough for the estimation.
#' @param maxEval Maximum number of function evaluations needed. The default is 10^6, which is same as the cubintegrate default.
#' @inheritParams freqpcr
#' @return A numeric vector of marginal likelihoods having the same length as \code{del}.
#' @family integrators
#' @section Methods for cubintegrate: The following methods are available for \code{\link[cubature]{cubintegrate}()}: cubmethod = c("hcubature", "pcubature", "cuhre", "divonne", "suave", "vegas").
#' \cr
#' \code{hcubature} is considerably fast, but less accurate with larger \code{reltol}. \code{cuhre} is moderately fast and the most accurate in most \code{relTol} range. If you can wait longer, \code{hcubature} with \code{relTol = 1e-4} or \code{cuhre} with \code{relTol = 1e-1} is recommended.
#' At \code{reltol = 1e-1}, \code{hcubature} is three times faster than \code{cuhre}, but the log-likelihood fluctuate by 0.1 (by 0.001 in \code{cuhre}).
#' The speed and accuracy of \code{hcubature} at \code{reltol = 1e-5} is comparable with \code{cuhre} at \code{reltol = 1e-1}.
#' \cr
#' \code{pcubature} and \code{divonne} frequently returns \code{NaN} and are not recommended.
#' \cr
#' \code{suave} and \code{vegas} are Monte Carlo integration and slow. They are certainly accurate even at \code{reltol = 1}, but interior to \code{cuhre} with same \code{reltol}.
#' @keywords internal
.integrate_beta <- function (   del, SHR, SHS, zeroAmount, targetScale, sdMeasure, xsm=2, EPCR,
                                cubmethod="hcubature", relTol=1e-1, absTol=1e-8, maxEval=10^6) {
    fa <- function(x, SHR, SHS, del, zeroAmount, targetScale, sdMeasure, xsm, EPCR) {
        dbeta(x, shape1=SHR, shape2=SHS) *
        dnorm(  del, 
                mean=-(log(targetScale)+log((zeroAmount+(1-zeroAmount)*x)))/log(1.0+EPCR), 
                sd=sqrt(xsm)*sdMeasure  )
    }
    return( cubintegrate(  f=fa, lower=0, upper=1, fDim=length(del),
                        method=cubmethod, relTol=relTol, absTol=absTol, maxEval=maxEval, nVec=1L,
                        SHR=SHR, SHS=SHS, del=del, zeroAmount=zeroAmount, targetScale=targetScale,
                        sdMeasure=sdMeasure, xsm=xsm, EPCR=EPCR)$integral )
}


#' @title Double integration of likelihood over the two DNA quantities obeying the gamma distributions.

#' @description Internal function to integrate the likelihood getting \eqn{\Delta}Cq value (the argument \code{del}) over the entire range of the DNA quantities of the two alleles, 0 <= x_S < Inf and 0 <= x_R < Inf. Vectorized for multiple bulk samples. It shares the arguments with .integrate_beta.
#' @inheritParams .integrate_beta
#' @return A numeric vector of marginal likelihoods having the same length as \code{del}.
#' @family integrators
#' @keywords internal
.integrate_gamma <- function (  del, SHR, SHS, zeroAmount, targetScale, sdMeasure, xsm=2, EPCR,
                                cubmethod="hcubature", relTol=1e-1, absTol=1e-8, maxEval=10^6) {
    # double integration of R and S DNA quantities  (x_R, x_S) over {0 to Inf, 0 to Inf}
    # transformed in a polar coordinates system (r, theta) over {0 to 1, 0 to pi/2}
    # xsm: For delta Cq, sdMeasure times two. For delta-delta Cq, sdMeasure times four.
    fa <- function(x=c(0, 0), SHR, SHS, del, zeroAmount, targetScale, sdMeasure, xsm, EPCR) {
        xS <- (x[1]/(1-x[1]))*cos(x[2]) # x[1] := r, x[2] := theta
        xR <- (x[1]/(1-x[1]))*sin(x[2])
        return(
            dgamma(xR, shape=SHR, scale=1) *
            dgamma(xS, shape=SHS, scale=1) *
            dnorm(  del,
                    mean=-(log(targetScale)+log(zeroAmount*xS+xR)-log(xS+xR))/log(1.0+EPCR),
                    sd=sqrt(xsm)*sdMeasure  ) *
            (x[1]/((1-x[1])^3))
        ) # when m=0 or m=n, SHR or SHS are always 0 and dgamma returns 0.
    }
    return( cubintegrate( f=fa, lower=c(0, 0), upper=c(1, pi/2), fDim=length(del),
                        method=cubmethod, relTol=relTol, absTol=absTol, maxEval=maxEval, nVec=1L,
                        SHR=SHR, SHS=SHS, del=del, zeroAmount=zeroAmount, targetScale=targetScale,
                        sdMeasure=sdMeasure, xsm=xsm, EPCR=EPCR)$integral )
}
