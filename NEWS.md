# Landmarking (development version)

## Changes that affect results

* Performance is now assessed with `riskRegression::Score()`. `prediction_error` contains the time-dependent AUC, the Brier score, their standard errors (from the influence function) and the index of prediction accuracy (IPA, a scaled Brier score). The C-index has been removed, as it is not a proper measure for evaluating predicted risks at a fixed time, and package `pec` is no longer used. Code that reads `prediction_error$c_index` needs updating.

* The Brier score was too low in earlier versions. It was calculated on follow-up censored at the horizon, so individuals who were event-free at the horizon were treated as censored before it. For the example data at landmark age 60, the Brier score is now 0.138 rather than 0.111 for the LOCF model, and 0.114 rather than 0.088 for the LME model (with the same predictions). The C-index was not affected.

* Cross-validation for LME models (`fit_LME_landmark()` with `k` or `cross_validation_df`) is now fully nested. For each fold, the LME model and the survival model are both fitted to the other folds only, and time is standardised using the other folds' measurements. Previously, the survival model for a fold was trained on values predicted by LME models that had been fitted with that fold's data (including its measurements after the landmark age when `include_data_after_x_L = TRUE`), which could make cross-validated performance optimistic. For the example data at landmark age 60 the difference is small (AUC 0.778 before and 0.782 now, with the corrected metrics). `model_LME_standardise_time` now has one element per fold. LOCF models, and LME models fitted without cross-validation, give the same results as before.

* `b` is deprecated: standard errors no longer use the bootstrap, and supplying `b` gives a warning.

## New features

* `plot()` has a new option `type = "smooth"`, which draws a smooth calibration curve estimated with pseudo-observations and gives the integrated calibration index (ICI), E50 and E90. The new function `get_calibration()` returns the curve and these metrics. The binned plot is still the default.

* `responses_LME` may contain missing values, for example when each response is measured at different times.

* `predict()` gives an error if `x_hor` is beyond the horizon the model was fitted with, or if `cv_fold` is not given for a model fitted with cross-validation. Fitted models store their horizon in `x_hor`.

* The documentation of `include_data_after_x_L` explains the trade-off of using measurements after the landmark age. The default is unchanged.

## Bug fixes and other changes

* Fix `fit_LME_landmark()` and `predict()` failing when `individual_id` is not `"id"` and random slopes are used in the survival submodel.

* Accept integer values for `x_L`, `x_hor`, `k`, `b` and `n` (e.g. `x_L = 60:61`).

* Fitted landmark models now store their resolved arguments in `args`. `predict()` and `plot()` use these instead of re-evaluating the original call, which failed when arguments were passed as variables. Models fitted with earlier versions still work.

* Fix `predict()` for LME models when `newdata` has unused factor levels in the `individual_id` column.

* A single `cross_validation_df` data frame can now be used with several landmark ages, and may be in long format.

* Missing cross-validation folds are now detected at every landmark age, not only the last.

* Add tests for argument handling, cross-validation folds, the `standard_cox` and `fine_gray` submodels, and `predict()` (including for a single new individual).

* Replace the RDS snapshot tests with numerical checks of predictions and prediction error against reference values.

* `R CMD check` now runs the test suite (adds `tests/testthat.R`). Slow tests are skipped on CRAN.

* Fix `fit_LME_landmark()` failing when a cross-validation fold has no individuals in a category of one of the `predictors_LME`. Character predictors are now converted to factors, so that every fold, and `predict()`, uses the same levels.

* `fit_LME_landmark()` now warns when individuals in the risk set are excluded from the survival submodel because they have no observation of one of the `responses_LME` on or before the landmark age. Their measurements are still used to fit the LME model. Results are unchanged.

* `fit_LOCF_landmark()` and `fit_LME_landmark()` now give a message when no cross-validation is used, as `prediction_error` is then the in-sample (apparent) performance.

* `fit_LOCF_landmark()`, `fit_LME_landmark()` and `predict()` now check that `x_hor` is greater than `x_L`.

* `add_cv_number()` no longer changes the global random number generator state. Fold assignments are unchanged, but bootstrap standard errors (argument `b`) now depend on the session's random seed: call `set.seed()` beforehand for reproducible standard errors.

* `plot()` now gives an informative error when there are too many tied predictions to form `n` groups, and its check for event types with no events now works.

* Internal fixes to `mixoutsamp()`: the checks of `newdata` now include the random-effects variables and look at the right columns.

* Documentation: bootstrap standard errors hold the fitted models fixed; `return_ids_with_LOCF()` includes measurements at `x_L`; the random slopes are per standard deviation of time when `standardise_time = TRUE`.

* Put the `NEWS.md` sections in order.

# Landmarking 1.0.2

* Change argument names in fit_LME_longitudinal() for clarity

* Fix bug to allow one longitudinal outcome in fit_LME_longitudinal(), in addition to two or more longitudinal outcomes

# Landmarking 1.0.1

* Fix bug that returned prediction of survival, rather than prediction of event, for the coxph survival model.

* Fix bug for non-competing risks scenario in get_model_assessment()

* Change the slopes in the survival model to be fixed slopes plus random slopes in function fit_LME_longitudinal()

* Add feature to use all the data or only data up to and including the landmark age in the validation dataset of function fit_LME_longitudinal()

# Landmarking 1.0.0

* Added a `NEWS.md` file to track changes to the package.
