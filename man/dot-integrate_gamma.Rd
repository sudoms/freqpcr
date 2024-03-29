% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_integrator.R
\name{.integrate_gamma}
\alias{.integrate_gamma}
\title{Double integration of likelihood over the two DNA quantities obeying the gamma distributions.}
\usage{
.integrate_gamma(
  del,
  SHR,
  SHS,
  zeroAmount,
  targetScale,
  sdMeasure,
  xsm = 2,
  EPCR,
  cubmethod = "hcubature",
  relTol = 0.1,
  absTol = 1e-08,
  maxEval = 10^6
)
}
\arguments{
\item{del}{Numeric vector of the observed \eqn{\Delta}Cq values.}

\item{SHR}{The gamma shape parameters for the mutant (R) portion of the bulk samples. Should be the same vector length as del. Each element of SHR is defined as K*(the assumed number of R allele in the bulk sample: 1, 2, 3, ..., n-1).}

\item{SHS}{The gamma shape parameters for the wild (S) portion of the bulk samples. Should be the same length as del. Each element of SHS is defined as K*(the assumed number of S allele in the bulk sample).}

\item{zeroAmount}{(In RED-\eqn{\Delta\Delta}Cq method) residue rate of restriction enzyme digestion, or (in general \eqn{\Delta\Delta}Cq analyses) small portion of the off-target allele on the target locus of the test sample, which will be amplified in the PCR. It needs to be always specified by the user as a number between 0 and 1, usually near 0.}

\item{targetScale}{(\eqn{\delta_{T}}) Scalar. The relative template DNA amount of the target locus to the houskeeping locus. If known, given as a positive numeric.}

\item{sdMeasure}{(\eqn{\sigma_{c}}) Scalar. The measurement error (standard deviation) on each Cq value following Normal(0, \eqn{\sigma_{c}^2}). If known, given as a positive numeric.}

\item{xsm}{Specify the accumulation of the standard deviation of the Cq measuring errors when the *-Cq values are fed as difference. For \eqn{\Delta}Cq values, sdMeasure times two. For \eqn{\Delta\Delta}Cq, sdMeasure times four. Default is two and used in most cases.}

\item{EPCR}{(\eqn{\eta}) Scalar. Amplification efficiency per PCR cycle. If known, given as a positive numeric. When \code{EPCR = 1}, template DNA doubles every cycle (\code{EPCR + 1 = 2}).}

\item{cubmethod}{Cubature method passed to the integrator function. See the section "Methods for cubintegrate".}

\item{relTol}{The maximum tolerance passed to the cubature method. Though the default of cubature::cubintegrate function is 1e-5, the accuracy is reduced here to acceralate the integration.}

\item{absTol}{The absolute tolerance passed to the cubature method. The default is 1e-8, which is less accurate than the default of cubintegrate function (1e-12) but considered enough for the estimation.}

\item{maxEval}{Maximum number of function evaluations needed. The default is 10^6, which is same as the cubintegrate default.}
}
\value{
A numeric vector of marginal likelihoods having the same length as \code{del}.
}
\description{
Internal function to integrate the likelihood getting \eqn{\Delta}Cq value (the argument \code{del}) over the entire range of the DNA quantities of the two alleles, 0 <= x_S < Inf and 0 <= x_R < Inf. Vectorized for multiple bulk samples. It shares the arguments with .integrate_beta.
}
\seealso{
Other integrators: 
\code{\link{.integrate_beta}()}
}
\concept{integrators}
\keyword{internal}
