# freqpcr
Interval estimation of population allele frequency based on quantitative PCR double-Delta Cq measures from bulk samples

# Features
Interval estimation of the population allele frequency from qPCR analysis based on the restriction enzyme digestion (RED)-DeltaDeltaCq method (Osakabe et al. 2017), as well as general DeltaDeltaCq analysis.
Compatible with the Cq measurement of DNA extracted from multiple individuals at once, so called ``group-testing'', this model assumes that the quantity of DNA extracted from an individual organism follows a gamma distribution.
Therefore, the point estimate is robust regarding the uncertainty of the DNA yield.

# Requirement
* R (>= ~~2.14~~3.6)
* cubature (https://cran.r-project.org/package=cubature)

# Installation
```
library(remotes)
install_github("sudoms/freqpcr")

library(freqpcr)
packageVersion("freqpcr")
```
If there are errors (converted from warning), it might be the case the dependent package 'cubature' has been built on a newer version of R (https://github.com/r-lib/remotes/issues/403).
```
** byte-compile and prepare package for lazy loading
Error: (converted from warning) package 'cubature' was built under R version 3.6.3
```
Then, set the following environment variable: R_REMOTES_NO_ERRORS_FROM_WARNINGS="true" and run install_github() again.
```
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
install_github("sudoms/freqpcr")

library(freqpcr)
packageVersion("freqpcr")
```

# Usage
```
library(freqpcr)
```

## First, define the parameters
```
P <- 0.15
K <- 4
sdMeasure <- 0.2

# These auxiliary parameters must be known
EPCR <- 0.97
zeroAmount <- 0.0016
```
P:
* The population allele frequency to be estimated

K:
* The gamma shape parameter of the individual DNA yield

sdMeasure:
* The measurement error on each Cq (Ct) value, following Normal(0, sdMeasure)

EPCR:
* Amplification efficiency per PCR cycle, given as a positive numeric.
* When EPCR = 1, template DNA doubles every cycle (EPCR + 1 = 2).

zeroAmount:
* (In RED-DeltaDelta Cq method) residual rate of restriction enzyme digestion.
* (in general DeltaDelta Cq analyses) small portion of the off-target allele amplified in the PCR process.
* It needs to be always specified by the user as a number between 0 and 1, usually near 0.

## Make a dummy Cq dataset
With four bulk samples, each of which comprises eight haploid individuals
```
dmy_cq <- make_dummy(   rand.seed=71, P=P, K=K, ntrap=4, npertrap=8,
                        scaleDNA=1e-06, targetScale=1.2, baseChange=0.2,
                        EPCR=EPCR, zeroAmount=zeroAmount,
                        sdMeasure=sdMeasure, diploid=FALSE  )
print(dmy_cq)
```

```
> print(dmy_cq)
An object of class "CqList"
Slot "N":
[1] 8 8 8 8

Slot "m":
     [,1] [,2] [,3] [,4]
[1,]    1    1    1    0
[2,]    7    7    7    8

Slot "xR":
[1] 6.654827e-07 1.540772e-06 7.998461e-07 0.000000e+00

Slot "xS":
[1] 5.541084e-06 7.909241e-06 8.052976e-06 1.004673e-05

Slot "housek0":
[1] 17.57120 16.90864 17.14117 16.88547

Slot "target0":
[1] 17.65665 16.84725 16.58306 16.90106

Slot "housek1":
[1] 19.90382 19.08542 19.59406 19.42896

Slot "target1":
[1] 23.13257 21.59854 22.97203 28.34278

Slot "DCW":
[1]  0.08544953 -0.06139074 -0.55811003  0.01559466

Slot "DCD":
[1] 3.228745 2.513126 3.377968 8.913814

Slot "deldel":
[1] 3.143296 2.574517 3.936078 8.898219

Slot "RFreqMeasure":
[1] 0.11868765 0.17453873 0.06933585 0.00239759

Slot "ObsP":
[1] 0.11868765 0.17453873 0.06933585 0.00239759

Slot "rand.seed":
[1] 71
```
## Estimation with freqpcr
* When P, K, targetScale, and baseChange are marked unknown.
```
result <- freqpcr(  N=dmy_cq@N, housek0=dmy_cq@housek0, target0=dmy_cq@target0,
                    housek1=dmy_cq@housek1, target1=dmy_cq@target1,
                    EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=2  )
print(result)
```

```
> print(result)
An object of class "CqFreq"
Slot "report":
                                                 Estimate Fixed  (scaled) (scaled.SE)       2.5%       97.5%
P (R-allele frequency)                         0.09773189     0 -2.222684  0.60914920 0.03178086   0.2633220
K (gamma shape parameter)                     20.92728472     0  3.041054  1.83522514 0.57354357 763.5884568
targetScale (relative amount of target locus)  1.11922896     0  0.112640  0.08911953 0.93985370   1.3328388
Cq measurement error (SD)                      0.20973065     0 -1.561931  0.32845069 0.11017528   0.3992451
EPCR (Duplication efficiency of PCR)           0.97000000     1        NA          NA         NA          NA

Slot "obj":
$minimum
[1] 6.094915

$estimate
[1] -2.222684  3.041054  0.112640 -1.561931

$gradient
[1] -3.401052e-05 -8.275831e-05 -5.170264e-05  8.878366e-05

$hessian
            [,1]        [,2]       [,3]       [,4]
[1,]  2.71023746  0.05094368   1.168535 -0.1766755
[2,]  0.05094368  0.37630056   2.045198 -0.6539723
[3,]  1.16853487  2.04519835 147.389557  6.4578192
[4,] -0.17667547 -0.65397225   6.457819 11.1504638

$code
[1] 1

$iterations
[1] 12


Slot "cal.time":
   user  system elapsed
   0.70    0.14    0.85
```
* When you have knowledge on the exact size of K (= set as a fixed parameter).
```
result <- freqpcr(  N=dmy_cq@N, housek0=dmy_cq@housek0, target0=dmy_cq@target0,
                    housek1=dmy_cq@housek1, target1=dmy_cq@target1,
                    K=4,
                    EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=1  )
```

# Note
* Download PDF help of the latest release https://github.com/sudoms/freqpcr/releases/latest
* Preprint is available from bioRxiv https://doi.org/10.1101/2021.01.19.427228
* Simulation dataset: https://figshare.com/account/home#/collections/5258027


# citation
```
citation("freqpcr")
```


# Author
Masaaki Sudo (https://orcid.org/0000-0001-9834-9857)
National Agriculture and Food Research Organization (NARO), Japan


# License
GNU GPL (>= 3)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# History
* v0.3.2 (2021.01.26) Fixed: R version ~~2.14~~ -> 3.6 and higher
* v0.3.1 (2021.01.20) First published (first draft of the preprint)
