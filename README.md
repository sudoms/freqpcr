# freqpcr
Interval estimation of population allele frequency based on quantitative PCR double-Delta Cq measures from bulk samples

<!--# DEMO
"hoge"の魅力が直感的に伝えわるデモ動画や図解を載せる-->

# Features
Interval estimation of the population allele frequency from qPCR analysis based on the restriction enzyme digestion (RED)-DeltaDeltaCq method (Osakabe et al. 2017), as well as general DeltaDeltaCq analysis. 
Compatible with the Cq measurement of DNA extracted from multiple individuals at once, so called ``group-testing'', this model assumes that the quantity of DNA extracted from an individual organism follows a gamma distribution. 
Therefore, the point estimate is robust regarding the uncertainty of the DNA yield.

# Requirement
* R (>= 2.14)
* cubature

# Installation
```
library(remotes)
install_github("sudoms/freqpcr")
```

# Usage
```
library(freqpcr)
```

## First, define the parameters.
```
# The population allele frequency to be estimated
P <- 0.75
# The gamma shape parameter of the individual DNA yield (to be estimated)
K <- 2
# The measurement error (SD) on each Cq (Ct) value following Normal(0, SD).
sdMeasure <- 0.3

# These auxiliary parameters must be known
EPCR <- 0.95
zeroAmount <- 0.0016
```

EPCR: 
* Amplification efficiency per PCR cycle, given as a positive numeric.
* When EPCR = 1, template DNA doubles every cycle (EPCR + 1 = 2).

zeroAmount:
* (In RED-DeltaDelta Cq method) residual rate of restriction enzyme digestion.
* (in general DeltaDelta Cq analyses) small portion of the off-target allele amplified in the PCR process.
* It needs to be always specified by the user as a number between 0 and 1, usually near 0.

## Make a dummy Cq dataset with six bulk samples, each of which comprises of eight haploid individuals.
```
dmy_cq <- make_dummy(   rand.seed=1, P=P, K=K, ntrap=6, npertrap=8,
                        scaleDNA=1e-07, targetScale=1.5, baseChange=0.3,
                        EPCR=EPCR, zeroAmount=zeroAmount,
                        sdMeasure=sdMeasure, diploid=FALSE  )
print(dmy_cq)
```
## Estimation with freqpcr
* When P, K, targetScale, and baseChange are marked unknown.
```
result <- freqpcr(  dmy_cq@N,
                    dmy_cq@housek0, dmy_cq@target0, dmy_cq@housek1, dmy_cq@target1,
                    EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=2  )
print(result)
```
* When you have knowledge on the size of baseChange (= set as a fixed parameter).
```
result <- freqpcr(  dmy_cq@N,
                    dmy_cq@housek0, dmy_cq@target0, dmy_cq@housek1, dmy_cq@target1,
                    baseChange=0.3,
                    EPCR=EPCR, zeroAmount=zeroAmount, beta=TRUE, print.level=1  )
```

# Note
* Preprint is soon available from bioRxiv
* Simulation dataset: https://figshare.com/account/home#/collections/5258027

# Author
Masaaki Sudo (https://orcid.org/0000-0001-9834-9857)
National Agriculture and Food Research Organization (NARO), Japan

# License
GNU GPL (>= 3)
