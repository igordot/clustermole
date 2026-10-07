## Test environments

* local R installation (macOS): R 4.6.1
* ubuntu-latest (on GitHub Actions), R 4.4, R-release, R-devel
* macOS (on GitHub Actions): R-release
* Windows (on GitHub Actions): R-release
* win-builder: R-devel

## R CMD check results

There were no ERRORs or WARNINGs.

## Resubmission

This is a resubmission. In this version I have:

* Moved GSVA, GSEABase, and singscore to Suggests to reduce the package load time.
* Skipped some slow tests on CRAN.
