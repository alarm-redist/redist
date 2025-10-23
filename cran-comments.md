## Test environments

* local R installation (Windows 11), R 4.5.1
* local R installation (macOS 11.4), R 4.5.1
* ubuntu-latest (on GitHub Actions), (oldrel-1, devel, and release)
* windows-latest (on GitHub Actions), (release)
* macOS-latest (on GitHub Actions), (release)
* Windows (on Winbuilder), (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

New maintainer:
  Christopher T. Kenny <ctkenny@proton.me>
Old maintainer(s):
  Christopher T. Kenny <christopherkenny@fas.harvard.edu>

## Reverse Dependencies

We checked 4 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

## Additional Notes

* This includes updating the maintainer email to ctkenny@proton.me. Please see the email from christopherkenny@fas.harvard.edu on 2025-08-29 for confirmation.
* Fixes CRAN issue regarding pareto score tests to avoid error on 1-dimensional frontier.
