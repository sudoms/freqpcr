% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_simulation.R
\name{sim_dummy}
\alias{sim_dummy}
\title{Simulate freqpcr estimation based on user-generated dummy data.}
\usage{
sim_dummy(
  CqList,
  EPCR,
  zeroAmount,
  P = NULL,
  K = NULL,
  targetScale = NULL,
  sdMeasure = NULL,
  beta,
  diploid,
  maxtime,
  print.level,
  aux = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{CqList}{Object belonging to the \linkS4class{CqList} class, typically the output from \code{\link{make_dummy}()}. Having the slots \code{N}, \code{housek0}, \code{target0}, \code{housek1}, and \code{target1}, all of which are numeric vectors of the same length.}

\item{EPCR}{(\eqn{\eta}) Amplification efficiency per PCR cycle, given as a positive numeric. When \code{EPCR = 1}, template DNA doubles every cycle (\code{EPCR + 1 = 2}).}

\item{zeroAmount}{A numeric between 0 and 1, usually near 0, giving the residue rate of restriction enzyme digestion in RED-\eqn{\Delta\Delta}Cq method.}

\item{P, K, targetScale, sdMeasure}{If NULL (default), the parameter is considered unknown and estimated via \code{\link{freqpcr}()}. If a value is specified, it is passed to \code{\link{freqpcr}()} as a fixed parameter. On the contrary, \code{EPCR} and \code{zeroAmount} are always treated as fixed parameters, for which values must be supplied.}

\item{beta, diploid, maxtime, print.level}{Configuration parameters which are passed directly to \code{\link{freqpcr}()}.}

\item{aux}{Additional information to be displayed on the console. The default is \code{NULL}. If some value is input by the user, it is echoed to stdout together with the contents of the argument \code{CqList}. This option is convenient when you want to record the original dummy dataset and the corresponding result sequentially e.g. using \code{capture.output()}.}

\item{verbose}{Prints more information e.g. system time. Default is \code{TRUE}.}

\item{...}{Additional arguments passed to \code{\link{freqpcr}()}.}
}
\value{
Object of the S4 class \linkS4class{CqFreq}, which is same as \code{\link{freqpcr}()}.
}
\description{
Wrapper of \code{\link{freqpcr}()} suitable for the performance test using a randomly-generated data object.
}
\examples{
# Prepare the parameter values.
K <- 2 # You already know the size of K in this case.
EPCR <- 0.97 # The sizes of EPCR and zeroAmount must always be supplied.
zeroAmount <- 1.6e-03
is.diploid <- FALSE

# First, make a dummy Cq dataset with six bulk DNA samples,
# each of which comprises of eight haploid individuals.
dmy_cq <- make_dummy( rand.seed=1, P=0.75, K=K, ntrap=6, npertrap=8, scaleDNA=1e-07,
                      targetScale=1.5, baseChange=0.3, EPCR=EPCR,
                      zeroAmount=zeroAmount, sdMeasure=0.3, diploid=is.diploid )

# Estimate the population allele frequency on the dummy dataset, presupposing K = 2.
sim_dummy(  CqList=dmy_cq, EPCR=EPCR, zeroAmount=zeroAmount,
            K=K,
            beta=TRUE, diploid=is.diploid, maxtime=60, print.level=2, aux="test" )

# If the maximum calculation time was too short to converge, nlm() returns error.
# Then sim_dummy() returns a matrix filled with zeros.
sim_dummy(  CqList=dmy_cq, EPCR=EPCR, zeroAmount=zeroAmount,
            beta=TRUE, diploid=is.diploid, maxtime=0.01, print.level=2 )
}
\seealso{
Other estimation procedures: 
\code{\link{freqpcr}()},
\code{\link{knownqpcr_unpaired}()},
\code{\link{knownqpcr}()}
}
\concept{estimation procedures}
