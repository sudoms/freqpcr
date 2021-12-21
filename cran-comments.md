## Test environments

* local: R version 3.6.1 x86_64-w64-mingw32/x64 (64-bit)
* win-builder: R version 4.1.2 x86_64-w64-mingw32 (64-bit)
* rhub: Fedora Linux, R-devel, clang, gfortran
* rhub: Ubuntu Linux 20.04.1 LTS, R-release, GCC

## R CMD check results

0 errors | 0 warnings | 2 notes 

```
Possibly mis-spelled words in DESCRIPTION:
  Cq (10:241)
  DeltaDeltaCq (4:26, 10:136, 10:198)
  Osakabe (10:157)
  al (10:168)
  et (10:165)
  qPCR (10:74)
```
> inst/WORDLIST was added.

```
Found the following (possibly) invalid URLs:
  URL: https://github.com/sudoms/freqpcr/README.jp.md
    From: inst/doc/freqpcr-intro.html
    Status: 404
    Message: Not Found
```
> The link is going to work after v0.4.0 is released and open to the public on github.com

## revdepcheck results

There are currently no downstream dependencies for this package