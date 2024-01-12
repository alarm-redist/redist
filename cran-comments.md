## Test environments
* local R installation (macOS), R 4.3.2
* local R installation (Windows), R 4.3.2
* ubuntu 20.04 (on GitHub Actions), (devel)
* ubuntu 22.04 (on GitHub Actions), (release)
* ubuntu 22.04 (on GitHub Actions), (old release)
* macOS-latest (on GitHub Actions), (release)
* Windows (on GitHub Actions), (release)
* Windows (on Winbuilder), (devel and release)


## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse Dependencies

There are no reverse dependencies to check.

## Additional Notes

* Fixes itemize braces Rd note for `redist.calc.frontier.size.Rd`, `redist.crsg.Rd`, and `redist.rsg.Rd`.
* Fixes "memory not mapped" error in `persily` function, which was due to be removed from the package in this release.
