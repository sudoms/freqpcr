
# freqpcr 0.4.0

* First submitted version to CRAN

* The paper (Sudo and Osakabe 2021) was published online (https://doi.org/10.1111/1755-0998.13554 )

* Vignette has been added

## New features

* The `freqpcr()` function has been extended to accept datasets lacking control samples (`housek0` and `target0`)

# freqpcr 0.3.5 

## Bug fixes
* (2021.02.17) `Lapack routine dgesv: system is exactly singular` in `knownqpcr()` when Cq contained NAs

# freqpcr 0.3.4 

* (2021.02.16) Released version of 0.3.3 (the code did not change)

# freqpcr 0.3.3

## New features

* (2021.02.13) Continuous distribution in sample allele ratio was supported. `knownqpcr()` function now officially accepts missing Cq observations (NA)

# freqpcr 0.3.2 

* (2021.01.29) Fixed: R version 2.14 -> 3.6 and higher

# freqpcr 0.3.1

* (2021.01.20) First published (first draft of the preprint)
