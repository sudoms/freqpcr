## Test environments

* local: R version 3.6.1 x86_64-w64-mingw32/x64 (64-bit)
* win-builder: R version 4.1.2 x86_64-w64-mingw32 (64-bit)
* rhub: Fedora Linux, R-devel, clang, gfortran
* rhub: Ubuntu Linux 20.04.1 LTS, R-release, GCC

## R CMD check results

0 errors | 0 warnings | 1 note 

```
New submission

Possibly mis-spelled words in DESCRIPTION:
  Cq (10:241)
  DeltaDeltaCq (4:26, 10:136, 10:198)
  Osakabe (10:157)
  al (10:168)
  et (10:165)
  qPCR (10:74)
```
> inst/WORDLIST was added.

## Revisions

> Title in DESCRIPTION was reduced to 65 characters.

> doi was added to the description.

> \dontrun{} in the example of `freqpcr` function was changed to \donttest{} as it takes 6 sec.

> Message output was re-designed and the 'quiet' argument in `knownqpcr` function was renamed as 'verbose'. This also affects the example code in the vignette (line 381).

## revdepcheck results

There are currently no downstream dependencies for this package