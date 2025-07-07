## Test environments

* local R installation (macOS), R 4.5.0
* local R installation (Windows), R 4.5.0
* ubuntu-latest (on GitHub Actions), (release)
* ubuntu-latest (on GitHub Actions), (old release 1)
* macOS-latest (on GitHub Actions), (release)
* Windows-latest (on GitHub Actions), (release)
* Windows-latest (on Winbuilder), (devel and release)


## R CMD check results

0 errors | 0 warnings | 1 note

This requires a change in maintainer email from `christopherkenny@fas.harvard.edu` to `ctkenny@proton.me`.
An email was sent to CRAN maintainers on 2025-07-07 referring to this change.

## Reverse Dependencies

There are no reverse dependencies to check.

## Additional Notes

* Fixes documentation links to external packages, such as correcting links from `tibble` to `tibble::tibble`.
* Corrects safe access to internet resources for vignettes.
