## Test environments
* local R installation (Windows), R 4.0.4
* local R installation (macOS), R 4.0.5
* ubuntu 20.04 (on GitHub Actions), (devel and release)
* windows-latest (on GitHub Actions), (release)
* macOS-latest (on GitHub Actions), (release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse Dependencies
There are no reverse dependencies to check.

## Additional Notes:
* This addresses the Solaris build issue by downgrading the PROJ/gdal within included data by several versions. RHub cannot support the full set of necessary packages for checking, so we have thinned out the packages to test and readded them as necessary for submission. Additionally, we have checked on a very old version of Ubuntu with identical PROJ/gdal setting to what Solaris has on CRAN. To the best of our knowledge, this should work, but if not, is there any way to ignore Solaris builds as long as it has no warnings or checks on the standard Linux/Windows/macOS builds?
* Fixes off by one error for clang-UBSAN and passed check with change on rhub::check_with_sanitizers()