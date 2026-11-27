# Landmarking (development version)

* Fix `fit_LME_landmark()` and `predict()` failing when `individual_id` is not `"id"` and random slopes are used in the survival submodel.

* Accept integer values for `x_L`, `x_hor`, `k`, `b` and `n` (e.g. `x_L = 60:61`).

* Fitted landmark models now store their resolved arguments in `args`. `predict()` and `plot()` use these instead of re-evaluating the original call, which failed when arguments were passed as variables. Models fitted with earlier versions still work.

* Fix `predict()` for LME models when `newdata` has unused factor levels in the `individual_id` column.

* A single `cross_validation_df` data frame can now be used with several landmark ages, and may be in long format.

* Missing cross-validation folds are now detected at every landmark age, not only the last.

* Add tests for argument handling, cross-validation folds, the `standard_cox` and `fine_gray` submodels, and `predict()` (including for a single new individual).

* Replace the RDS snapshot tests with numerical checks of predictions and prediction error against reference values.

# Landmarking 1.0.2

# Landmarking 1.0.0

* Added a `NEWS.md` file to track changes to the package.

# Landmarking 1.0.1

* Fix bug that returned prediction of survival, rather than prediction of event, for the coxph survival model.

* Fix bug for non-competing risks scenario in get_model_assessment()

* Change the slopes in the survival model to be fixed slopes plus random slopes in function fit_LME_longitudinal()

* Add feature to use all the data or only data up to and including the landmark age in the validation dataset of function fit_LME_longitudinal()

# Landmarking 1.0.2

* Change argument names in fit_LME_longitudinal() for clarity

* Fix bug to allow one longitudinal outcome in fit_LME_longitudinal(), in addition to two or more longitudinal outcomes
