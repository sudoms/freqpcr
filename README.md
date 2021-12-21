# freqpcr
Interval estimation of population allele frequency based on quantitative PCR ΔΔCq measures from bulk samples

* [README in Japanese（日本語: equivalent to this document + vignette）](https://github.com/sudoms/freqpcr/blob/master/README.jp.md)

# Features
Interval estimation of the population allele frequency from qPCR analysis based on the restriction enzyme digestion (RED)-ΔΔCq method (Osakabe et al. 2017), as well as general  ΔΔCq analysis.
Compatible with the Cq measurement of DNA extracted from multiple individuals at once, so called ``group-testing'', this model assumes that the quantity of DNA extracted from an individual organism follows a gamma distribution.
Therefore, the point estimate is robust regarding the uncertainty of the DNA yield.

See the package vignette [Introduction to freqpcr] and [published article (Sudo and Osakabe 2021 Molecular Ecology Resources)](https://doi.org/10.1111/1755-0998.13554) for the model structure.

Prerequisites are:

* Want to know the percentage of the specific allele at the population level (population allele frequency) relative to other allele(s) on the same locus.

* For the population in question, multiple DNA solutions (bulk samples) have been obtained, each of which consists of multiple individuals extracted at once.

* Quantitative PCR based on the real-time PCR analysis has been performed on each bulk sample and a set of four Cq values is available.

* The ΔΔCq value calculated from the Cq quartet represents the allele frequency of the bulk sample (sample allele frequency).

The main function `freqpcr()` takes as input data the number of individuals included in each bulk sample (`N`) and the set of Cq values (`housek0`, `target0`, `housek1`, `target1`) measured for each bulk sample (since it is a ΔΔCq method, there are usually four Cq values).

Prior to the estimation with samples with unknown allele ratios, auxiliary experimental parameters e.g. amplification efficiency during real-time PCR are required. Functions for estimating experimental parameters (`knownqpcr()` and `knownqpcr_unpaired()` functions) are also provided, using multiple DNA solutions with known allele mixing ratios.

> Since v0.4.0, the `freqpcr()` function has been extended to accept datasets lacking control samples (`housek0` and `target0`) i.e., a ΔCq method provided the sizes of the auxiliary experimental parameters are known.


# Requirement
* R (>= 3.6)
* cubature (https://cran.r-project.org/package=cubature)

# Installation

```R
library(remotes)
install_github("sudoms/freqpcr")

library(freqpcr)
packageVersion("freqpcr")
```

If there are errors (converted from warning), it might be the case the dependent package 'cubature' has been built on a newer version of R (https://github.com/r-lib/remotes/issues/403).

```R
** byte-compile and prepare package for lazy loading
Error: (converted from warning) package 'cubature' was built under R version 3.6.3
```

Then, set the following environment variable: `R_REMOTES_NO_ERRORS_FROM_WARNINGS="true"` on you R session and run install_github() again.

```R
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
install_github("sudoms/freqpcr")
```

# Usage

```R
library(freqpcr)
```

## First, define the parameters for dummy data generation

Dummy Cq data is generated for demo use. As for real qPCR experiments, you first need to estimate those auxiliary parameters using samples with known multiple allele-mixing ratios and `knownqpcr()` or `knownqpcr_unpaired()` functions. See the vignette for detail.

```R
P <- 0.15
K <- 4

# These auxiliary parameters must be pre-defined 
EPCR <- 0.97
zeroAmount <- 0.0016
scaleDNA <- 1e-06 # used in make_dummy() but not in freqpcr()
targetScale <- 1.2 
baseChange <- 0.2# used in make_dummy() but not in freqpcr()
sdMeasure <- 0.2
```

P:
* The population allele frequency to be estimated

K:
* The gamma shape parameter of the individual DNA yield

EPCR:
* Amplification efficiency per PCR cycle, given as a positive numeric.
* When EPCR = 1, template DNA doubles every cycle (EPCR + 1 = 2).

zeroAmount:
* (In RED-ΔΔCq method) residual rate of restriction enzyme digestion.
* (in general ΔΔCq and ΔCq analyses) small portion of the off-target allele amplified in the PCR process.
* It needs to be always specified by the user as a number between 0 and 1, usually near 0.

scaleDNA:
* Relative amount of the template DNA (of the housekeeping gene) compared to the threshold of realtime PCR. 
* If DNA increases twofold in a PCR cycle (`EPCR = 1`), `scaleDNA` = 10^-6 means Cq value is approximately 20.

targetScale:
* Relative DNA amount of the target locus compared with the housekeeping gene.

baseChange:
* Relative amount of the template DNA (of the housekeeping gene) before and after the restriction enzyme treatment in the RED-ΔΔCq method (as an involuntary change due to experimental manipulation).
* The variable does not exist in other common ΔΔCq and ΔCq methods.

sdMeasure:
* The measurement error on each Cq (Ct) value, following Normal(0, sdMeasure)

## Make a dummy Cq dataset

With four bulk samples, each of which comprises eight haploid individuals

```R
dmy_cq <- make_dummy(   rand.seed=71, P=P, K=K, ntrap=4, npertrap=8,
                        scaleDNA=scaleDNA, targetScale=targetScale, 
                        baseChange=targetScale, EPCR=EPCR, 
                        zeroAmount=zeroAmount,
                        sdMeasure=sdMeasure, diploid=FALSE  )
print(dmy_cq)
```

```R
> print(dmy_cq)
An object of class "CqList"
Slot "N":
[1] 8 8 8 8

Slot "m":
     [,1] [,2] [,3] [,4]
[1,]    1    1    1    0
[2,]    7    7    7    8

Slot "xR":
[1] 2.661931e-06 6.163087e-06 3.199384e-06 0.000000e+00

Slot "xS":
[1] 2.216434e-05 3.163697e-05 3.221190e-05 4.018694e-05

Slot "housek0":
[1] 15.52662 14.86406 15.09659 14.84089

Slot "target0":
[1] 15.61207 14.80267 14.53848 14.85648

Slot "housek1":
[1] 17.85924 17.04084 17.54948 17.38438

Slot "target1":
[1] 21.08799 19.55396 20.92745 26.29820

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

* When `P`, `K`, `sdMeasure`, and `targetScale` are marked unknown.

```R
result <- freqpcr(  N=dmy_cq@N, housek0=dmy_cq@housek0, target0=dmy_cq@target0,
                    housek1=dmy_cq@housek1, target1=dmy_cq@target1,
                    EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=2  )
print(result)
```

```R
> print(result)
An object of class "CqFreq"
Slot "report":
                                                 Estimate Fixed  (scaled) (scaled.SE)       2.5%       97.5%
P (R-allele frequency)                         0.09773189     0 -2.222684  0.60914923 0.03178086   0.2633220
K (gamma shape parameter)                     20.92728471     0  3.041054  1.83522515 0.57354355 763.5884712
targetScale (relative amount of target locus)  1.11922896     0  0.112640  0.08911953 0.93985370   1.3328388
Cq measurement error (SD)                      0.20973065     0 -1.561931  0.32845070 0.11017528   0.3992451
EPCR (Duplication efficiency of PCR)           0.97000000     1        NA          NA         NA          NA

Slot "obj":
$minimum
[1] 6.094915

$estimate
[1] -2.222684  3.041054  0.112640 -1.561931

$gradient
[1] -3.400973e-05 -8.275889e-05 -5.170087e-05  8.878422e-05

$hessian
            [,1]        [,2]       [,3]       [,4]
[1,]  2.71023719  0.05094362   1.168535 -0.1766756
[2,]  0.05094362  0.37630056   2.045198 -0.6539723
[3,]  1.16853469  2.04519835 147.389558  6.4578190
[4,] -0.17667556 -0.65397228   6.457819 11.1504636

$code
[1] 1

$iterations
[1] 12


Slot "cal.time":
   user  system elapsed
   0.72    0.15    0.88
```
* When you have knowledge on the exact size of K (= set as a fixed parameter).

```R
result <- freqpcr(  N=dmy_cq@N, housek0=dmy_cq@housek0, target0=dmy_cq@target0,
                    housek1=dmy_cq@housek1, target1=dmy_cq@target1,
                    K=4,
                    EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=1  )
```


# External links

* GitHub repository https://github.com/sudoms/freqpcr

* PDF help of the latest version https://github.com/sudoms/freqpcr/releases/latest

* Paper https://doi.org/10.1111/1755-0998.13554 
<br>
Sudo, M., & Osakabe, M. (2021). freqpcr: Estimation of population allele frequency using qPCR ΔΔCq measures from bulk samples. Molecular Ecology Resources, 00, 1–14. Creative Commons Attribution 4.0 International (CC BY 4.0)

* Author website https://sudori.info/english.html


# citation

```R
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

* NEWS.md
