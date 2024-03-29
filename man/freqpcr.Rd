% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_freqpcr.R
\name{freqpcr}
\alias{freqpcr}
\title{Estimate population allele frequency from the set of Cq measurements.}
\usage{
freqpcr(
  N,
  A,
  housek0,
  target0,
  housek1,
  target1,
  P = NULL,
  K = NULL,
  targetScale = NULL,
  sdMeasure = NULL,
  EPCR = 0.99,
  XInit0 = c(P = NULL, K = NULL, targetScale = NULL, sdMeasure = NULL, EPCR = NULL),
  zeroAmount = NULL,
  beta = TRUE,
  diploid = FALSE,
  pvalue = 0.05,
  gradtol = 1e-04,
  steptol = 1e-09,
  iterlim = 100,
  maxtime = 600,
  print.level = 1,
  ...
)
}
\arguments{
\item{N}{Sample sizes as a numeric vector. \code{N[i]} signifies the number of individuals (both for haploidy and diploidy) contained in the \emph{i}th bulk sample. \code{N} must not contain a missing value (\code{NA}). If \code{N} is not applicable (= even not 1), feed \code{A} instead of \code{N} and then the estimation process runs with the `continuous model'.}

\item{A}{Use instead of \code{N} in the continuous model. This is a scale factor to control the relative DNA content between samples. \code{A[i]} can take any positive value, but must not be \code{NA}. Considering the case you have arranged each sample by e.g. water filrtation or extraction from a culture in a petri dish, it is convenient to define the unit size of \code{A[i]} == 1.0 to be same as the vessel volume (e.g. 2.0 for two petri dishs, 0.5 for half bottle of water, etc.). When neither \code{N} nor \code{A} is specified by the user, the function stops. If both \code{N} and \code{A} are specified, only \code{N} is evaluated.}

\item{housek0}{A numeric vector. In RED-\eqn{\Delta\Delta}Cq method, \code{housek0} is the Cq values of the test sample without the restriction enzyme digestion, which is amplified with the primer set for a housekeeping gene. In general \eqn{\Delta\Delta}Cq analyses, \code{housek0} is defined for the control sample (typically, 100\% mutant) solution, which is also amplified with the primer set for the housekeeping gene. 
\cr
Since \code{v0.4.0}, you can run \code{\link{freqpcr}()} without specifying \code{housek0} and \code{target0} (a \eqn{\Delta}Cq method). As this setting halves the effective data points, it is recommended to fix other parameters, especially \code{targetScale}. 
\cr
The four Cq arguments, \code{housek0}, \code{target0}, \code{housek1}, and \code{target1}, all must have the same data length. They also must be the same length as \code{N} or \code{A}. If the Cq dataset has missing values, they must be filled with NA so that the length of the data vectors will not differ.}

\item{target0}{In RED-\eqn{\Delta\Delta}Cq method, \code{target0[i]} signifies the measured Cq value of the \emph{i}th bulk sample without the digestion, for which both alleles, wild-type (S: susceptible) and mutant (R: resistant to a pesticide), on the target locus are amplified. In general \eqn{\Delta\Delta}Cq analyses, \code{target0} is the Cq values of the pure-S control sample, which is amplified with a R-allele-specific primer set.}

\item{housek1}{The Cq values of the test sample measured on the housekeeping gene after the restriction enzyme digestion (in RED-\eqn{\Delta\Delta}Cq method), or the test sample amplified on the housekeeping gene (in general \eqn{\Delta\Delta}Cq analyses).}

\item{target1}{For each test sample with unknown allele-ratio, \code{target1[i]} is defined as the Cq value for the target locus amplified after the restriction enzyme digestion (in RED-\eqn{\Delta\Delta}Cq method), or the target locus amplified with the R-allele-specific primer set (in general \eqn{\Delta\Delta}Cq analyses).}

\item{P}{Scalar. Population allele frequency from which the test samples are derived. Default is \code{NULL} and to be estimated. If the parameter is known, it is given as a numeric between 0 and 1.}

\item{K}{Scalar. The gamma shape parameter of the individual DNA yield. Default is \code{NULL} and to be estimated. If known, given as a positive numeric.}

\item{targetScale}{(\eqn{\delta_{T}}) Scalar. The relative template DNA amount of the target locus to the houskeeping locus. If known, given as a positive numeric.}

\item{sdMeasure}{(\eqn{\sigma_{c}}) Scalar. The measurement error (standard deviation) on each Cq value following Normal(0, \eqn{\sigma_{c}^2}). If known, given as a positive numeric.}

\item{EPCR}{(\eqn{\eta}) Scalar. Amplification efficiency per PCR cycle. If known, given as a positive numeric. When \code{EPCR = 1}, template DNA doubles every cycle (\code{EPCR + 1 = 2}).}

\item{XInit0}{Optionally the initial value for the parameter optimization can be specified, but it is strongly recommended to keep the argument as is. Unlike \code{XInit} in \code{\link{knownqpcr}()}, the argument is not directly passed to the optimizer; used only when each parameter is set unknown (the parameter is absent or specified as NULL).}

\item{zeroAmount}{(In RED-\eqn{\Delta\Delta}Cq method) residue rate of restriction enzyme digestion, or (in general \eqn{\Delta\Delta}Cq analyses) small portion of the off-target allele on the target locus of the test sample, which will be amplified in the PCR. It needs to be always specified by the user as a number between 0 and 1, usually near 0.}

\item{beta}{Whether to use the beta distribution to approximate the sample allele ratio instead of specifying individual gamma distribution for each of the allelic DNA amounts? Default is \code{TRUE}, which accelerates the calculation.}

\item{diploid}{Is the target organism diploidy? Default is \code{FALSE}, assuming haploidy. Current implementation of diploidy assumes i.i.d. between the amounts of R and S chromosomes owned by a heterozygote individual, which is unlikely in many animals but necessary for the calculation in a realistic time.}

\item{pvalue}{The two-sided confidence interval is calculated at the last iteration at given significance level. Default is 0.05, which returns the 95\% Wald's CI (2.5 to 97.5 percentile) based on the Hessian matrix.}

\item{gradtol, steptol, iterlim}{Control parameters passed to \code{\link[stats]{nlm}()}. \code{gradtol} and \code{steptol} are the positive scalars giving the tolerance to terminate the algorithm and the minimum allowable relative step length. \code{iterlim} specifies the maximum number of iterations to be performed before the program is terminated (and evaluated at the last iteration). Usually 30 iterations are enough.}

\item{maxtime}{A positive scalar to set the maximum calculation time in seconds to abort the optimizer (and return error). The total calculation time largely depends on \code{N[i]}, the number of individuals contained in each bulk sample.}

\item{print.level}{\code{print.level=1} (the default) shows the initial values of the parameters and likelihood as well as the output in the last iteration. \code{print.level=2} shows the parameter values and gradients in every step. \code{print.level=0} does not output any intermediate state to R console, simply returning the result summary.}

\item{...}{Additional arguments passed to the function.}
}
\value{
Object of the S4 class \linkS4class{CqFreq}. The slot \code{report} is a matrix and each row contains the estimated parameter value with 100*(1-pvalue)\% confidence interval. The following parameters are returned:
\enumerate{
    \item\code{P}, the population allele frequency from which the test samples are derived.
    \item\code{K}, the gamma shape parameter of the individual DNA yield.
    \item\code{targetScale} (\eqn{\delta_{T}}), the relative template DNA amount of the target to the houskeeping loci.
    \item\code{EPCR} (\eqn{\eta}), the amplification efficiency per PCR cycle.
    \item\code{sdMeasure} or "Cq measurement error" (\eqn{\sigma_{c}}).
}
}
\description{
The function estimates the population allele frequency using the dataset of Cq values measured over \code{length(N)} bulk samples, each of which has a sample size of \code{N[i]} as the number of individuals included. \code{N[i]} can be 1, meaning that every individual is analyzed separately. For the \emph{i}th sample, the four Cq values were measured as \code{housek0[i]}, \code{target0[i]}, \code{housek1[i]}, and \code{target1[i]}. The function can estimate up to five parameters simultaneously when the Cq sets are available for more than two (bulk) samples: \code{P}, \code{K}, \code{targetScale}, \code{sdMeasure}, and \code{EPCR}.
\cr
Since v0.3.2, user can also use an experimental `continuous model' by specifying \code{A} instead of \code{N}. That is, each sample DNA is directly extracted from the environment and the sample allele ratio \code{y} follows \code{y ~ Beta(apk, a(1-p)k)} instead of \code{y ~ Beta(mk, (n-m)k), m ~ Binomial(n, p)}, where \code{p} and \code{k} are the population allele frequency and the gamma shape parameter of the individual DNA yield, respectively. Each element of \code{A}, \code{a} is a scaling factor of relative DNA contents between the samples. The continuous model is likely when each sample directly comes from the environment e.g., water sampling in an eDNA analysis or cell culture in a petri dish.
\cr
Since v0.4.0, \code{\link{freqpcr}()} also works without specifying \code{housek0} and \code{target0} i.e., it can estimate population allele frequency from \eqn{\Delta}Cq values instead of \eqn{\Delta\Delta}Cq. In this setting, the sizes of \code{targetScale} and \code{sdMeasure} should be fixed in addition to \code{EPCR} and \code{zeroAmount}.
}
\section{Choise of the parameters to be estimated}{

Estimation is conducted only for parameters where the values are not specified or specified explicitly as \code{NULL}. If one feeds a value e.g. \code{K=1} or \code{sdMeasure=0.24}, it is then treated as fixed parameter. Notwithstanding, \code{EPCR} is estimated only when \code{EPCR = NULL} is specified explicitly.
\cr
You must verify the size of \code{EPCR} and \code{zeroAmount} beforehand because they are not estimable from the experiments with unknown allele ratios. Although \code{targetScale} and \code{sdMeasure} are estimable within \code{\link{freqpcr}()}, it is better to feed the known values, especially when you have only a few bulk samples (length(N) <= 3). Fixing \code{targetScale} and \code{sdMeasure} is also strongly recommended when \code{housek0} and \code{target0} are absent (\eqn{\Delta}Cq method). The functions \code{\link{knownqpcr}()} or \code{\link{knownqpcr_unpaired}()}, depending on your data format, provide the procedure to estimate the sizes of the experimental parameters using the DNA solutions of known allele mixing ratios. 
\cr
For the unknown parameters, \code{XInit0} optionally specifies initial values for the optimization using \code{\link[stats]{nlm}()} though the use of internal default is strongly recommended. The specification as a fixed parameter has higher priority than the specification in \code{XInit0}. Every user-specified parameter values, fixed or unknown, must be given in linear scale (e.g. between 0 and 1 for the allele frequency); they are transformed internally to log or logit.
}

\examples{
\donttest{
# Dummy Cq dataset with six bulk samples, each of which comprises of eight haploids.
EPCR <- 0.95; zeroAmount <- 0.0016; # True values for the two must be known.
P <- 0.75
dmy_cq <- make_dummy(   rand.seed=1, P=P, K=2, ntrap=6, npertrap=8,
                        scaleDNA=1e-07, targetScale=1.5, baseChange=0.3,
                        EPCR=EPCR, zeroAmount=zeroAmount,
                        sdMeasure=0.3, diploid=FALSE   )
print(dmy_cq)

# Estimation with freqpcr, where P, K, targetScale, and baseChange are marked unknown.
result <- freqpcr( N=dmy_cq@N, housek0=dmy_cq@housek0, target0=dmy_cq@target0,
                   housek1=dmy_cq@housek1, target1=dmy_cq@target1,
                   EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=0 )
print(result)

# Estimation with freqpcr, assuming diploidy.
result <- freqpcr( N=dmy_cq@N, housek0=dmy_cq@housek0, target0=dmy_cq@target0,
                   housek1=dmy_cq@housek1, target1=dmy_cq@target1,
                   EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, diploid=TRUE )

# Estimation where you have knowledge on the size of K.
result <- freqpcr( N=dmy_cq@N, housek0=dmy_cq@housek0, target0=dmy_cq@target0,
                   housek1=dmy_cq@housek1, target1=dmy_cq@target1,
                   K=2,
                   EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=2 )
# (>= v0.3.2)
# Provided the model is continuous (specify A instead of N).
result <- freqpcr( A=dmy_cq@N, housek0=dmy_cq@housek0, target0=dmy_cq@target0,
                   housek1=dmy_cq@housek1, target1=dmy_cq@target1,
                   K=2, EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=1 )
# (>= v0.4.0)
# If the dataset lacks control samples (housek0 and target0 are absent).
# Fixing the sizes of targetScale and sdMeasure is strongly recommended.
result <- freqpcr( N=dmy_cq@N, housek1=dmy_cq@housek1, target1=dmy_cq@target1,
                   K=2, EPCR=EPCR, zeroAmount=zeroAmount,
                   targetScale=1.5, sdMeasure=0.3, beta=TRUE, print.level=1 )
}
}
\seealso{
Other estimation procedures: 
\code{\link{knownqpcr_unpaired}()},
\code{\link{knownqpcr}()},
\code{\link{sim_dummy}()}
}
\concept{estimation procedures}
