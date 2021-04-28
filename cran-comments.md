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
* Spaces were removed from the doi and arXiv links in the description to make links clickable.
* Fixes Rd files with missing return flags in Roxygen documentation
* dontrun is changed to donttest for cases where examples can be reliably run but will take significant time. In many cases, the dontrun was dropped completely when the example completes in under 5 seconds. In the case of the enumpart family of functions, we keep to dontrun, as these can take significant time to run and will fail if system software requirements declared in the DESCRIPTION file are not met.
* enumpart examples are shifted to save to tempdir(), rather than the home directory to comply with CRAN policy.
* The change to the random seed kind in test setup was removed.