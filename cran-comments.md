## Test environments
* local R installation (Windows), R 4.0.4
* local R installation (macOS), R 4.0.5
* ubuntu 20.04 (on GitHub Actions), (devel and release)
* windows-latest (on GitHub Actions), (release)
* macOS-latest (on GitHub Actions), (release)

## R CMD check results

0 errors | 0 warnings | 1 notes

    installed size is  7.1Mb
    sub-directories of 1Mb or more:
      enumpart   2.0Mb
      libs       2.5Mb

## Reverse Dependencies
There are no reverse dependencies to check.

## Additional Notes:
* This is to fix the CRAN check issues with version 3.0.1. In particular, the OGR error in building vignettes is fixed by downgrading the PROJ string that was included in the data with the package. The tests are now properly backwards compatible with changes in matrices. The Solaris compiler issues were resovled by clarifying some C++ types that were ambiguous.

* On some builds, we are getting a NOTE that the installed directory is 6 or more Mbs.