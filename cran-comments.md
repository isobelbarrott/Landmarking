## Update

This is an update from version 1.0.2 to 1.1.0. The main changes are:

* Prediction performance is now assessed with `riskRegression::Score()` (time-dependent AUC, Brier score and index of prediction accuracy, with standard errors), so package `pec` is no longer used.
* The Brier score, which was underestimated in earlier versions, has been corrected.
* Cross-validation for the LME landmark model is now fully nested.
* A smooth calibration curve and calibration metrics are available (`get_calibration()`, `plot(type = "smooth")`).
* `R CMD check` now runs the test suite; slow tests are skipped on CRAN.

See NEWS.md for details.

## Test environments

* local Windows 11 install (R release)
* GitHub Actions: windows-latest (R release), macos-latest (R release), ubuntu-latest (R devel, R release and R oldrel-1)
* win-builder (R devel)
* mac builder (R release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no reverse dependencies.
